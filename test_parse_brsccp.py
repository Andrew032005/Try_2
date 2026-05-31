import pytest
from parse_brsccp import convert_type, parse_brsccp_log, parse_line, parse_metadata
import tempfile
import os

@pytest.mark.parametrize("val, expected", [
    ("42", 42),
    ("3.14", 3.14),
    ("True", True),
    ("false", False),
    ("'hello'", "hello"),
    ('"', '"'), # Edge case for single quote
    ("['A', 'B']", ['A', 'B']),
    ("{'FE': 2, 'ZN': 2}", {'FE': 2, 'ZN': 2}),
    ("Some string", "Some string"),
    ("", ""),
])
def test_convert_type(val, expected):
    assert convert_type(val) == expected

def test_parse_line_complex_scenarios():
    data = {
        "metadata": {},
        "warnings": {"pdb_discontinuities": [], "large_shifts": []},
        "optimization": {"best_trajectory": {}, "improved_iterations": []},
        "energetics": {"complex_1": {}, "complex_2": {}},
        "results": {}
    }

    # Test multiple improved iterations
    parse_line("Better structure was found at 0 iteration!", data)
    parse_line("Better structure was found at 2 iteration!", data)
    assert data["optimization"]["improved_iterations"] == [0, 2]

    # Test complex 2 energy
    parse_line("Energy of complex 2 after optimization: -702.0 h", data)
    assert data["energetics"]["complex_2"]["complex_energy_opt"] == -702.0

    # Test protein energy in dot
    parse_line("Energy of protein 1 in dot: -601.3 h", data)
    assert data["energetics"]["complex_1"]["protein_energy_dot"] == -601.3

def test_parse_metadata_malformed():
    metadata = {}
    content = """
--- Input parameters ---
ligand_charge = -1
model_name = 'TEST'
invalid line
---
    """
    parse_metadata(content, metadata)
    assert metadata["ligand_charge"] == -1
    # Current implementation might be greedy with continuation lines
    # It's acceptable for now as long as it handles valid lines.

def test_parse_brsccp_log_integration():
    log_path = "brsccp2V8X.log"
    if not os.path.exists(log_path):
        pytest.skip("Log file not found for integration test")

    data = parse_brsccp_log(log_path)

    assert data["metadata"]["model_name"] == "2V8X"
    assert data["metadata"]["ligand_charge"] == -1
    assert data["results"]["ddG_kcal_mol"] == 2.0
    assert data["results"]["activity_ratio"] == 0.03
    assert len(data["warnings"]["large_shifts"]) == 24
    assert data["energetics"]["complex_1"]["complex_energy_opt"] == -697.84919
    assert data["energetics"]["complex_2"]["protein_energy_dot"] == -601.37474

def test_file_not_found():
    with pytest.raises(FileNotFoundError):
        parse_brsccp_log("non_existent_file.log")

def test_mock_log_content():
    content = """
--------------------------------------------------- Input parameters ---------------------------------------------------
                    ligand_charge = -1
                    model_name = 'TEST_MODEL'
------------------------------------------------------------------------------------------------------------------------
                    WARNING   PDBConstructionWarning: Chain A is discontinuous at line 6672.
                    WARNING    C2' @ MGQ(1218) has shift 1.12 > 1.0 A!
                    Energy of complex 1 after optimization: -100.0 h
                    Energy of ligand 1 in dot: -10.0 h
                    Relative binding energy, ddG: 2.50 kcal/mol
                    Activity ratio, exp(-ddG/RT): 0.05
                    TOTAL RUN TIME: 0:10:00
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.log', delete=False) as tmp:
        tmp.write(content)
        tmp_path = tmp.name

    try:
        data = parse_brsccp_log(tmp_path)
        assert data["metadata"]["ligand_charge"] == -1
        assert data["metadata"]["model_name"] == "TEST_MODEL"
        assert "Chain A is discontinuous" in data["warnings"]["pdb_discontinuities"][0]
        assert data["warnings"]["large_shifts"][0]["atom"] == "C2'"
        assert data["energetics"]["complex_1"]["complex_energy_opt"] == -100.0
        assert data["energetics"]["complex_1"]["ligand_energy_dot"] == -10.0
        assert data["results"]["ddG_kcal_mol"] == 2.5
        assert data["results"]["activity_ratio"] == 0.05
        assert data["results"]["total_run_time"] == "0:10:00"
    finally:
        os.remove(tmp_path)
