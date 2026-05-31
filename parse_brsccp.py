import re
import json
import argparse
import ast
from pathlib import Path

# Regular Expressions
RE_PDB_WARNING = re.compile(r"WARNING\s+PDBConstructionWarning:\s*(.*)")
RE_SHIFT_WARNING = re.compile(r"WARNING\s+(\w+['\"]?)\s+@\s+(\w+)\((\d+)\)\s+has shift\s+([\d.]+)")
RE_BEST_TRAJ = re.compile(r"Best optimization trajectory:\s*(\w+);\s*E=([-.\d]+)\s*Eh")
RE_IMPROVED_ITER = re.compile(r"Better structure was found at (\d+) iteration!")
RE_ENERGY = re.compile(r"Energy of (complex|ligand|protein) (1|2) (?:after optimization|in dot):\s*([-.\d]+)")
RE_DDG = re.compile(r"Relative binding energy, ddG.*:\s*([-.\d]+)\s+kcal/mol")
RE_ACTIVITY = re.compile(r"Activity ratio,.*:\s*([-.\d]+)")
RE_RUNTIME = re.compile(r"TOTAL RUN TIME:\s*(.*)")

def parse_metadata(content, data_dict):
    """Extract metadata from Input parameters section."""
    start_marker = "--- Input parameters ---"
    if start_marker in content:
        start_idx = content.find(start_marker) + len(start_marker)
        # Find the end of the section marked by a line of dashes
        end_match = re.search(r"\n\s*-{3,}", content[start_idx:])
        if end_match:
            end_idx = start_idx + end_match.start()
        else:
            end_idx = len(content)

        params_text = content[start_idx:end_idx]

        current_key = None
        current_value = []

        for line in params_text.splitlines():
            if not line.strip(): continue
            # Match "key = value" optionally preceded by whitespace
            match = re.search(r"^\s*(\w+)\s*=\s*(.*)$", line)
            if match:
                if current_key:
                    data_dict[current_key] = convert_type(" ".join(current_value))
                current_key = match.group(1)
                current_value = [match.group(2).strip()]
            elif current_key:
                current_value.append(line.strip())

        if current_key:
            data_dict[current_key] = convert_type(" ".join(current_value))

def convert_type(val):
    """Helper to convert string values to Python types."""
    try:
        val = val.strip()
        if len(val) >= 2 and ((val.startswith("'") and val.endswith("'")) or (val.startswith('"') and val.endswith('"'))):
            return val[1:-1]
        elif val.lower() == 'true':
            return True
        elif val.lower() == 'false':
            return False
        else:
            try:
                return ast.literal_eval(val)
            except (ValueError, SyntaxError):
                return val
    except Exception:
        return val

def parse_line(line, data):
    """Parse a single line and update the data dictionary."""
    # PDB Warnings
    pdb_match = RE_PDB_WARNING.search(line)
    if pdb_match:
        data["warnings"]["pdb_discontinuities"].append(pdb_match.group(1).strip())
        return

    # Shift Warnings
    shift_match = RE_SHIFT_WARNING.search(line)
    if shift_match:
        data["warnings"]["large_shifts"].append({
            "atom": shift_match.group(1),
            "residue": shift_match.group(2),
            "res_id": int(shift_match.group(3)),
            "shift_value": float(shift_match.group(4))
        })
        return

    # Optimization
    traj_match = RE_BEST_TRAJ.search(line)
    if traj_match:
        data["optimization"]["best_trajectory"] = {
            "name": traj_match.group(1),
            "energy_Eh": float(traj_match.group(2))
        }
        return

    iter_match = RE_IMPROVED_ITER.search(line)
    if iter_match:
        data["optimization"]["improved_iterations"].append(int(iter_match.group(1)))
        return

    # Energetics
    energy_match = RE_ENERGY.search(line)
    if energy_match:
        etype = energy_match.group(1) # complex, ligand, protein
        eidx = energy_match.group(2)  # 1, 2
        eval_val = float(energy_match.group(3))

        target_key = f"complex_{eidx}"
        subkey = f"{etype}_energy_dot" if "in dot" in line else f"{etype}_energy_opt"
        data["energetics"][target_key][subkey] = eval_val
        return

    # Results
    ddg_match = RE_DDG.search(line)
    if ddg_match:
        data["results"]["ddG_kcal_mol"] = float(ddg_match.group(1))
        return

    act_match = RE_ACTIVITY.search(line)
    if act_match:
        data["results"]["activity_ratio"] = float(act_match.group(1))
        return

    runtime_match = RE_RUNTIME.search(line)
    if runtime_match:
        data["results"]["total_run_time"] = runtime_match.group(1).strip()
        return

def parse_brsccp_log(file_path):
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Error: File {file_path} not found.")

    log_content = path.read_text(encoding='utf-8')

    data = {
        "metadata": {},
        "warnings": {"pdb_discontinuities": [], "large_shifts": []},
        "optimization": {"best_trajectory": {}, "improved_iterations": []},
        "energetics": {"complex_1": {}, "complex_2": {}},
        "results": {}
    }

    parse_metadata(log_content, data["metadata"])

    for line in log_content.splitlines():
        parse_line(line, data)

    return data

def main():
    parser = argparse.ArgumentParser(description="Parse BRSCCP log file into structured JSON.")
    parser.add_argument("log_file", help="Path to the BRSCCP log file")
    parser.add_argument("-o", "--output", help="Path to the output JSON file", default="result.json")

    args = parser.parse_args()

    try:
        parsed_data = parse_brsccp_log(args.log_file)
        with open(args.output, "w", encoding='utf-8') as f:
            json.dump(parsed_data, f, indent=4)

        print(f"Successfully parsed {args.log_file} -> {args.output}")

        results = parsed_data.get("results", {})
        metadata = parsed_data.get("metadata", {})

        print("\n--- Summary ---")
        print(f"Replacement: {metadata.get('replacement_type', 'N/A')}")
        print(f"ddG: {results.get('ddG_kcal_mol', 'N/A')} kcal/mol")
        print(f"Activity ratio: {results.get('activity_ratio', 'N/A')}")
        print(f"Total Run Time: {results.get('total_run_time', 'N/A')}")
        print(f"Large shifts found: {len(parsed_data['warnings']['large_shifts'])}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
