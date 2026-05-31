# Testing Strategy for BRSCCP Log Parser

This document outlines the testing strategy for the `parse_brsccp.py` script.

## 1. Unit Testing
Unit tests focus on individual functions to ensure they behave correctly in isolation.
- **`convert_type(val)`**: Verified with various inputs including integers, floats, booleans, lists, dictionaries, and quoted strings.
- **`parse_line(line, data)`**: Verified by passing individual log lines and checking if the `data` dictionary is updated correctly for warnings, energetics, and results.
- **`parse_metadata(content, data_dict)`**: Verified with well-formed and malformed "Input parameters" blocks.

## 2. Integration Testing
Integration tests ensure that the entire parsing pipeline works together from file reading to the final structured output.
- **Mock Log Content**: A complete mock log is parsed using `tempfile` to verify that all sections are correctly extracted in a single pass.
- **Real Log Integration**: The script is run against the provided `brsccp2V8X.log` to verify its performance on real-world data.

## 3. Edge Case Testing
- **Missing Files**: Verifies that `FileNotFoundError` is raised.
- **Missing Sections**: Verifies that the script doesn't crash if certain log sections (like metadata or energetics) are missing.
- **Malformed Data**: Verifies that the script handles unexpected characters or missing values gracefully.

## 4. How to Run Tests
Tests are implemented using `pytest`. To run them, use:

```bash
pytest test_parse_brsccp.py
```

For more detailed output:

```bash
pytest -v test_parse_brsccp.py
```
