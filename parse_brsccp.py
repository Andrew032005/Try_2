import re
import json
import argparse
import ast
from pathlib import Path

def parse_brsccp_log(file_path):
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Error: File {file_path} not found.")

    log_content = path.read_text(encoding='utf-8')

    # Initialize data structure
    data = {
        "metadata": {},
        "warnings": {
            "pdb_discontinuities": [],
            "large_shifts": []
        },
        "optimization": {
            "best_trajectory": {},
            "improved_iterations": []
        },
        "energetics": {
            "complex_1": {},
            "complex_2": {}
        },
        "results": {}
    }

    # Regular Expressions
    re_pdb_warning = re.compile(r"WARNING\s+PDBConstructionWarning:\s*(.*)")
    re_shift_warning = re.compile(r"WARNING\s+(\w+['\"]?)\s+@\s+(\w+)\((\d+)\)\s+has shift\s+([\d.]+)")
    re_best_traj = re.compile(r"Best optimization trajectory:\s*(\w+);\s*E=([-.\d]+)\s*Eh")
    re_improved_iter = re.compile(r"Better structure was found at (\d+) iteration!")
    re_energy = re.compile(r"Energy of (complex|ligand|protein) (1|2) (?:after optimization|in dot):\s*([-.\d]+)")
    re_ddg = re.compile(r"Relative binding energy, ddG.*:\s*([-.\d]+)\s+kcal/mol")
    re_activity = re.compile(r"Activity ratio,.*:\s*([-.\d]+)")
    re_runtime = re.compile(r"TOTAL RUN TIME:\s*(.*)")

    # Metadata extraction
    # Find the block starting with "Input parameters" and ending with a long dash line
    # We use a more robust regex to find the section
    input_params_section_match = re.search(r"-+\s*Input parameters\s*-+\n(.*?)\n\s*-{10,}", log_content, re.DOTALL | re.IGNORECASE)
    if input_params_section_match:
        params_text = input_params_section_match.group(1)

        current_key = None
        current_value = []

        for line in params_text.splitlines():
            if not line.strip(): continue

            # Check if line contains '='
            match = re.match(r"^\s*(\w+)\s*=\s*(.*)$", line)
            if match:
                if current_key:
                    save_metadata(data["metadata"], current_key, " ".join(current_value))
                current_key = match.group(1)
                current_value = [match.group(2).strip()]
            elif current_key:
                current_value.append(line.strip())

        if current_key:
            save_metadata(data["metadata"], current_key, " ".join(current_value))

    lines = log_content.splitlines()
    for line in lines:
        # PDB Warnings
        pdb_match = re_pdb_warning.search(line)
        if pdb_match:
            data["warnings"]["pdb_discontinuities"].append(pdb_match.group(1).strip())
            continue

        # Shift Warnings
        shift_match = re_shift_warning.search(line)
        if shift_match:
            data["warnings"]["large_shifts"].append({
                "atom": shift_match.group(1),
                "residue": shift_match.group(2),
                "res_id": int(shift_match.group(3)),
                "shift_value": float(shift_match.group(4))
            })
            continue

        # Optimization
        traj_match = re_best_traj.search(line)
        if traj_match:
            data["optimization"]["best_trajectory"] = {
                "name": traj_match.group(1),
                "energy_Eh": float(traj_match.group(2))
            }
            continue

        iter_match = re_improved_iter.search(line)
        if iter_match:
            data["optimization"]["improved_iterations"].append(int(iter_match.group(1)))
            continue

        # Energetics
        energy_match = re_energy.search(line)
        if energy_match:
            etype = energy_match.group(1) # complex, ligand, protein
            eidx = energy_match.group(2)  # 1, 2
            eval_val = float(energy_match.group(3))

            target_key = f"complex_{eidx}"

            # Map type to subkey
            if "in dot" in line:
                subkey = f"{etype}_energy_dot"
            else:
                subkey = f"{etype}_energy_opt"

            data["energetics"][target_key][subkey] = eval_val
            continue

        # Results
        ddg_match = re_ddg.search(line)
        if ddg_match:
            data["results"]["ddG_kcal_mol"] = float(ddg_match.group(1))
            continue

        act_match = re_activity.search(line)
        if act_match:
            data["results"]["activity_ratio"] = float(act_match.group(1))
            continue

        runtime_match = re_runtime.search(line)
        if runtime_match:
            data["results"]["total_run_time"] = runtime_match.group(1).strip()
            continue

    return data

def save_metadata(metadata_dict, key, val):
    try:
        if (val.startswith("'") and val.endswith("'")) or (val.startswith('"') and val.endswith('"')):
            inner = val[1:-1]
            metadata_dict[key] = inner
        elif val.lower() == 'true':
            metadata_dict[key] = True
        elif val.lower() == 'false':
            metadata_dict[key] = False
        else:
            try:
                metadata_dict[key] = ast.literal_eval(val)
            except (ValueError, SyntaxError):
                metadata_dict[key] = val
    except Exception:
        metadata_dict[key] = val

def main():
    parser = argparse.ArgumentParser(description="Parse BRSCCP log file into structured JSON.")
    parser.add_argument("log_file", help="Path to the BRSCCP log file")
    parser.add_argument("-o", "--output", help="Path to the output JSON file", default="result.json")

    args = parser.parse_args()

    parsed_data = parse_brsccp_log(args.log_file)

    if parsed_data:
        with open(args.output, "w", encoding='utf-8') as f:
            json.dump(parsed_data, f, indent=4)

        print(f"Successfully parsed {args.log_file} -> {args.output}")

        # Short report
        results = parsed_data.get("results", {})
        metadata = parsed_data.get("metadata", {})

        print("\n--- Summary ---")
        print(f"Replacement: {metadata.get('replacement_type', 'N/A')}")
        print(f"ddG: {results.get('ddG_kcal_mol', 'N/A')} kcal/mol")
        print(f"Activity ratio: {results.get('activity_ratio', 'N/A')}")
        print(f"Total Run Time: {results.get('total_run_time', 'N/A')}")
        print(f"Large shifts found: {len(parsed_data['warnings']['large_shifts'])}")

if __name__ == "__main__":
    main()
