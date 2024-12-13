# Build_indices

## Quick reference

Copy and paste to pull this image:

#### `docker pull us.gcr.io/broad-gotc-prod/build-indices:1.0.0-2.7.10a-1663605340`

- __What is this image:__ This image is a Debian-based custom image with STAR installed and pre-configured along with python scripts to build indices.
- __What is STAR:__ Spliced Transcripts Alignment to a Reference (STAR) is a fast RNA-seq read mapper, with support for splice-junction and fusion read detection. STAR aligns reads by finding the Maximal Mappable Prefix (MMP) hits between reads (or read pairs) and the genome, using a Suffix Array index, [more info here](https://github.com/alexdobin/STAR).
- __How to see tool version used in image:__ Please see below.

## Versioning

Build_indices uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/build-indices:<image-version>-<star-version>-<unix-timestamp>`

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/build-indices:1.0.0-2.7.10a-1663605340
$ docker inspect us.gcr.io/broad-gotc-prod/build-indices:1.0.0-2.7.10a-1663605340
```

## Usage

### Build_indices Docker Container

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/build-indices:1.0.0-2.7.10a-1663605340 \
    build-indices bash
```

Then you can exec into the container and use STAR or any of the scripts accordingly. Alternatively, you can run one-off commands by passing the command as a docker run parameter.

## GTF Comparison Tools

This repository includes tools for comparing and testing GTF (Gene Transfer Format) file modifications. These tools ensure consistency in GTF processing and provide detailed comparison reports.

### Components

#### Scripts
- `compare_gtfs.py` - Analyzes differences between two GTF files
- `test_gtf_comparison.py` - Unit tests for GTF comparison functionality
- `modify_gtf.py` - Script to modify GTF files

#### Required Files
- `test_data/test1.gtf` - Test GTF file
- `Biotypes.tsv` - File containing allowed biotypes

### Features

The comparison tool analyzes:
- Structural differences in GTF fields
- Attribute differences, including:
  - Reordered attributes
  - Extra or missing attributes
  - Different attribute values
- Gene-level differences
- Mitochondrial gene comparisons

### Running GTF Comparison

```bash
python compare_gtfs.py <gtf1> <gtf2> --output-prefix <prefix>
```

Example:
```bash
python compare_gtfs.py test_data/test1.gtf modified_output.gtf --output-prefix comparison
```

### Testing

Run the test suite:
```bash
python -m unittest test_gtf_comparison.py -v
```

### GitHub Actions Integration

Automated testing is configured via GitHub Actions:
- Runs comparison tests
- Generates reports
- Uploads test artifacts

Configuration file: `.github/workflows/gtf_tests.yml`

### Output Reports

1. Structural Differences (`<prefix>_structural_diff.txt`):
   - Row counts
   - Field differences
   - Sample comparisons

2. Attribute Differences (`<prefix>_attribute_diff.txt`):
   - Attribute summaries
   - Detailed comparisons
   - Value differences

3. Gene Differences (`<prefix>_gene_diff.txt`):
   - Gene counts
   - Unique gene lists
   - MT gene analysis

### Requirements

- Python 3.x
- pandas
- Standard Python libraries

Install dependencies:
```bash
pip install pandas
```

### Directory Structure

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

### Error Handling

The tools include comprehensive error handling for:
- Missing files
- Malformed GTF content
- Directory issues
- Attribute parsing errors

### Contributing

When modifying these tools:
1. Ensure all tests pass
2. Update test cases for new features
3. Maintain Docker compatibility
4. Update documentation
5. Follow GitHub Actions workflow requirements

## Notes

- GTF comparison is sensitive to format variations
- Docker container provides consistent environment
- All scripts are accessible within the container
- Use reference files for reliable testing