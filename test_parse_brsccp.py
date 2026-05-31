import pytest
from parse_brsccp import save_metadata, parse_brsccp_log
import tempfile
import os

def test_save_metadata_types():
    metadata = {}

    # Int
    save_metadata(metadata, "test_int", "42")
    assert metadata["test_int"] == 42

    # Float
    save_metadata(metadata, "test_float", "3.14")
    assert metadata["test_float"] == 3.14

    # Bool
    save_metadata(metadata, "test_bool_true", "True")
    assert metadata["test_bool_true"] is True
    save_metadata(metadata, "test_bool_false", "false")
    assert metadata["test_bool_false"] is False

    # String
    save_metadata(metadata, "test_str", "'hello'")
    assert metadata["test_str"] == "hello"

    # List
    save_metadata(metadata, "test_list", "['A', 'B']")
    assert metadata["test_list"] == ['A', 'B']

    # Dict
    save_metadata(metadata, "test_dict", "{'FE': 2, 'ZN': 2}")
    assert metadata["test_dict"] == {'FE': 2, 'ZN': 2}

    # Multi-line list (joined by space)
    save_metadata(metadata, "test_multiline", "['CL', 'PO4', 'SO4']")
    assert metadata["test_multiline"] == ['CL', 'PO4', 'SO4']

def test_parse_brsccp_log_integration():
    # Since we have the actual log file, we can test against it
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

def test_missing_sections():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.log', delete=False) as tmp:
        tmp.write("Some random text\nRelative binding energy, ddG: 1.5 kcal/mol\n")
        tmp_path = tmp.name

    try:
        data = parse_brsccp_log(tmp_path)
        assert data["results"]["ddG_kcal_mol"] == 1.5
        assert data["metadata"] == {}
        assert data["warnings"]["large_shifts"] == []
    finally:
        os.remove(tmp_path)

def test_complex_idx_switch():
    content = """
                    Energy of complex 1 after optimization: -100.0 h
                    Energy of complex 2 after optimization: -200.0 h
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.log', delete=False) as tmp:
        tmp.write(content)
        tmp_path = tmp.name

    try:
        data = parse_brsccp_log(tmp_path)
        assert data["energetics"]["complex_1"]["complex_energy_opt"] == -100.0
        assert data["energetics"]["complex_2"]["complex_energy_opt"] == -200.0
    finally:
        os.remove(tmp_path)

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
