# GTF File Comparison Tools

This directory contains tools for comparing and testing GTF (Gene Transfer Format) file modifications. The tools are designed to work together to ensure consistency in GTF file processing and to detect and report differences between GTF files.

## Components

### Scripts
- `compare_gtfs.py` - Main comparison tool for analyzing differences between two GTF files
- `test_gtf_comparison.py` - Unit tests for GTF comparison functionality
- `modify_gtf.py` - Script to modify GTF files (referenced in tests)

### Required Files
- `test_data/test1.gtf` - Test GTF file
- `Biotypes.tsv` - File containing allowed biotypes

## Features

The comparison tool analyzes:
- Structural differences in the first 8 GTF fields
- Attribute differences in the 9th field, including:
  - Reordered attributes
  - Extra or missing attributes
  - Different attribute values
- Gene-level differences
- Mitochondrial gene comparisons

## Usage

### Running GTF Comparison

```bash
python compare_gtfs.py <gtf1> <gtf2> --output-prefix <prefix>
```

Example:
```bash
python compare_gtfs.py test_data/test1.gtf modified_output.gtf --output-prefix comparison
```

This will generate three output files:
- `<prefix>_structural_diff.txt` - Differences in GTF structure
- `<prefix>_attribute_diff.txt` - Detailed attribute differences
- `<prefix>_gene_diff.txt` - Gene-level comparison results

### Running Tests

```bash
python -m unittest test_gtf_comparison.py -v
```

## GitHub Actions Integration

The repository includes GitHub Actions workflows that automatically:
1. Run GTF comparison tests
2. Generate comparison reports
3. Upload test artifacts

### Workflow Files
- `.github/workflows/gtf_tests.yml` - Main test workflow configuration

## Output Reports

### Structural Differences Report
Contains information about:
- Total row counts
- Row-by-row field differences
- Sample differences for each field

### Attribute Differences Report
Shows:
- Summary of attribute differences
- Detailed attribute comparisons
- Extra attributes in each file
- Value differences for common attributes

### Gene Differences Report
Includes:
- Total gene counts
- Unique genes in each file
- Mitochondrial gene analysis

## Requirements

- Python 3.x
- pandas
- Standard Python libraries (argparse, os, collections)

## Installation

No special installation required. Just ensure you have the required Python packages:

```bash
pip install pandas
```

## Directory Structure

```
build-indices/
├── test_data/
│   ├── test1.gtf
│   └── reference_outputs/
├── test_output/
│   └── comparison_files/
├── compare_gtfs.py
├── test_gtf_comparison.py
├── modify_gtf.py
└── Biotypes.tsv
```

## Contributing

When modifying these scripts:
1. Ensure all tests pass
2. Update test cases for new functionality
3. Maintain compatibility with GitHub Actions workflow
4. Update documentation as needed

## Error Handling

The scripts include comprehensive error handling for:
- Missing input files
- Malformed GTF content
- Directory creation/access issues
- Attribute parsing errors

## Output Examples

Example attribute difference report:
```
Attribute Differences (9th field):

Attribute Key           Only in GTF1    Only in GTF2    Different Values
-------------------------------------------------------------------------
gene_id                 0              0              5
gene_type              2              0              3
transcript_id          1              0              0

Detailed Attribute Differences:
Row 1:
  Different values:
    gene_id: GTF1="ENSG01", GTF2="ENSG01.1"
```

## Note

This comparison tool is sensitive to:
- GTF format variations
- Attribute ordering
- Whitespace differences
- Version numbers in IDs

Make sure your input files follow standard GTF formatting for best results.
