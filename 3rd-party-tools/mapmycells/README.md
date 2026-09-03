# MapMyCells

## Quick reference

```bash
docker pull us.gcr.io/broad-gotc-prod/mapmycells:1.0.0-1.7.4-1725283419
```

- **What is this image:** Python 3.11 based image wrapping the MapMyCells `cell_type_mapper` library.
- **What is `cell_type_mapper`:** A toolkit from the Allen Institute for mapping single-cell transcriptomics data to hierarchical reference taxonomies.

## Synopsis

This image provides the `cell_type_mapper` command-line interfaces (`map_to_on_the_fly_markers` and `from_specified_markers`). It will be used by WARP to classify cell types in single-cell and single-nucleus RNA-seq datasets (e.g., Human MTG and Mouse WMB) by comparing query expression data against precomputed reference atlases.

## Image contents

- Base image: `python:3.11-slim`
- `cell_type_mapper` `1.7.4`
- Baked-in reference-atlas assets under `/opt/mapmycells/data/` (sourced from Allen Institute's
  public MapMyCells S3 bucket at build time, see the Dockerfile), so the WDL doesn't need to fetch
  them from external storage at runtime:
  - `precomputed_stats_ABC_revision_230821.h5` + `mouse_markers_230821.json` — Whole Mouse Brain / ABC atlas taxonomy
  - `precomputed_stats.20231120.sea_ad.MTG.h5` — Human MTG / SEA-AD taxonomy (markers are found at runtime for this one)
  - No gene-symbol-to-Ensembl-ID mapping db is baked in yet — building it needs a 15GB NCBI
    taxonomy download (`bkbit download-ncbi-taxonomy`), which doesn't fit a routine image build

## Versioning

This image follows the WARP-Tools tag convention
`us.gcr.io/broad-gotc-prod/mapmycells:<image-version>-<tool-version>-<unix-timestamp>`.
Past versions are tracked in [docker_versions.tsv](docker_versions.tsv).

## Inputs

| Input | Type / format | Required | Description |
| --- | --- | --- | --- |
| `--query_path` | `.h5ad` | yes | AnnData file containing raw transcriptomics counts |
| `--precomputed_stats.path` | `.h5` | yes | HDF5 file containing the precomputed hierarchical taxonomy stats |
| `--gene_mapping.db_path` | `.db` | no | SQLite database for translating gene symbols to Ensembl IDs |
| `--query_markers.serialized_lookup` | `.json` | no | JSON file containing predefined marker genes |

## Outputs

| Output | Format | Description |
| --- | --- | --- |
| `output.csv` | `CSV` | Cell type predictions mapping query cells to reference taxonomy |
| `output.json` | `JSON` | Extended summary results |

## Usage

Run an interactive shell in the image:

```bash
docker run --rm -it us.gcr.io/broad-gotc-prod/mapmycells:1.0.0-1.7.4-1725283419 bash
```

Run the `from_specified_markers` workflow:
```bash
docker run --rm -v $PWD:/data us.gcr.io/broad-gotc-prod/mapmycells:1.0.0-1.7.4-1725283419 \
    python -m cell_type_mapper.cli.from_specified_markers \
    --query_path /data/query.h5ad \
    --extended_result_path /data/output.json \
    --csv_result_path /data/output.csv \
    --precomputed_stats.path /data/stats.h5 \
    --query_markers.serialized_lookup /data/markers.json \
    --type_assignment.normalization raw
```

## Task descriptions

- **Pipeline:** `MapMyCells.wdl` (In development)
- **Task(s):** `RunMapMyCells` — Invokes `python -m cell_type_mapper.cli.*` to classify cells.

## Troubleshooting

To run the image independently of a WDL for testing, start a shell explicitly
(`docker run -it --rm <image-url> bash`). See
[BUILDING.md → Troubleshooting and running standalone](../../BUILDING.md#troubleshooting-and-running-standalone).
