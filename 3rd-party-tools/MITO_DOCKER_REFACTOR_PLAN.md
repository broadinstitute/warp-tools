# Mitochondria Docker Refactor Plan

## Goal

Consolidate the mtSwirl Python code from `warp/all_of_us/mitochondria/mtSwirl_refactor/` into
`warp-tools`, reduce from 3 Dockerfiles to 2, and ensure both images are self-contained and
descriptively named.

---

## Current State (3 Dockerfiles)

| # | Location | Base image | Purpose |
|---|----------|------------|---------|
| 1 | `warp-tools/3rd-party-tools/aou_mitochondrial_annotate_coverage/Dockerfile` | `hailgenetics/hail:0.2.127-py3.11` | Hail + gnomad + mtSwirl (git cloned). Used only as base for Image 3. Otherwise unused. |
| 2 | `warp/all_of_us/mitochondria/mtSwirl_refactor/coverage_db/Dockerfile` | `python:3.12-slim` | Lightweight — numpy/pandas/h5py/gcsfs only. Builds the coverage DB (HDF5). Standalone. |
| 3 | `warp/all_of_us/mitochondria/mtSwirl_refactor/Terra/Dockerfile` | Image 1 (above) | Image 1 + gcloud CLI + h5py/gcsfs + local Python COPYd in. Combines VCFs + homref from covdb. |

**Problem:** Image 1 exists only to serve as a base for Image 3. Images 2 and 3 both live
outside warp-tools. Python source is split across warp and warp-tools.

---

## Target State (2 Dockerfiles, all in warp-tools)

### Image A — `aou_mito_coverage_db`
- **Replaces:** Image 2 (`coverage_db/Dockerfile`)
- **Base:** `python:3.12-slim`
- **Purpose:** Builds the HDF5 coverage database. Standalone, no Hail dependency.
- **Scripts needed:** `coverage_db/build_coverage_db.py` only (no cross-package imports)

### Image B — `aou_mito_hail_processing`
- **Replaces:** Image 1 + Image 3 (combined into one)
- **Base:** `hailgenetics/hail:0.2.127-py3.11` (directly — no intermediate base image)
- **Purpose:** All Hail-based mtSwirl steps: sharding, merging, combining VCFs, finalizing MT, adding annotations
- **Scripts needed:** All root-level shared modules + Terra/ scripts (see below)

---

## Script Ownership (no shared scripts between images)

### Image A scripts (`aou_mito_coverage_db/scripts/`)
```
scripts/
└── coverage_db/
    ├── __init__.py
    ├── build_coverage_db.py
    ├── smoke_test_build_coverage_db.py
    └── smoke_test_build_coverage_db_skip_summary.py
```

### Image B scripts (`aou_mito_hail_processing/scripts/`)
```
scripts/
├── __init__.py
├── add_annotations.py
├── covdb_utils.py
├── merging_constants.py
├── merging_utils.py
└── Terra/
    ├── __init__.py
    ├── build_vcf_shard_mt.py
    ├── combine_vcfs_and_homref_from_covdb.py
    ├── finalize_mt_with_covdb.py
    ├── merge_mt_shards.py
    ├── shard_mt_by_samples.py
    └── union_mt_shards.py
```

Source for all scripts: `warp/all_of_us/mitochondria/mtSwirl_refactor/`

---

## New warp-tools Directory Structure

```
warp-tools/3rd-party-tools/
├── aou_mito_coverage_db/
│   ├── Dockerfile                          # adapted from coverage_db/Dockerfile
│   ├── scripts/
│   │   └── coverage_db/
│   │       ├── __init__.py
│   │       ├── build_coverage_db.py
│   │       ├── smoke_test_build_coverage_db.py
│   │       └── smoke_test_build_coverage_db_skip_summary.py
│   └── aou_mito_coverage_db.changelog.md
│
├── aou_mito_hail_processing/
│   ├── Dockerfile                          # merged from Image 1 + Image 3
│   ├── scripts/
│   │   ├── __init__.py
│   │   ├── add_annotations.py
│   │   ├── covdb_utils.py
│   │   ├── merging_constants.py
│   │   ├── merging_utils.py
│   │   └── Terra/
│   │       ├── __init__.py
│   │       ├── build_vcf_shard_mt.py
│   │       ├── combine_vcfs_and_homref_from_covdb.py
│   │       ├── finalize_mt_with_covdb.py
│   │       ├── merge_mt_shards.py
│   │       ├── shard_mt_by_samples.py
│   │       └── union_mt_shards.py
│   └── aou_mito_hail_processing.changelog.md
│
└── aou_mitochondrial_annotate_coverage/    # DELETE — replaced by aou_mito_hail_processing
    └── Dockerfile
```

---

## Dockerfile Changes

### `aou_mito_coverage_db/Dockerfile`
- Copy `coverage_db/Dockerfile` as-is
- Update `COPY` path to reflect new build context:
  - Old: `COPY ./ /opt/mtSwirl/generate_mtdna_call_mt/`
  - New: `COPY scripts/ /opt/mtSwirl/generate_mtdna_call_mt/`
- Update build instructions comment to reflect new location

### `aou_mito_hail_processing/Dockerfile`
Merge Image 1 and Image 3 into a single Dockerfile:
1. Start `FROM hailgenetics/hail:0.2.127-py3.11` (Image 1's base)
2. Keep all of Image 1's `RUN` layers (pip installs, git clones for gnomad_qc + gnomad-mitochondria)
3. Remove the `git clone https://github.com/jamesemery/mtSwirl.git` step — replaced by `COPY`
4. Add Image 3's additions:
   - gcloud CLI install block
   - `pip install h5py gcsfs`
5. Add `COPY scripts/ /opt/mtSwirl/generate_mtdna_call_mt/`
6. Set `PYTHONPATH=/opt/mtSwirl` and `ENTRYPOINT ["python3"]`

---

## WDL / Task Updates Required

Any WDL tasks in the `warp` repo that reference the old image names will need updating:
- `us.gcr.io/broad-gotc-prod/aou-mitochondrial-annotate-coverage:*` → `us.gcr.io/broad-gotc-prod/aou-mito-hail-processing:*`
- Search: `grep -rn "aou-mitochondrial-annotate-coverage\|aou_mitochondrial_annotate_coverage" warp/`

---

## Steps to Execute

1. Create `warp-tools/3rd-party-tools/aou_mito_coverage_db/` and copy scripts
2. Create `warp-tools/3rd-party-tools/aou_mito_hail_processing/` and copy scripts
3. Write `aou_mito_coverage_db/Dockerfile` (adapted from `coverage_db/Dockerfile`)
4. Write `aou_mito_hail_processing/Dockerfile` (merged from Image 1 + Image 3)
5. Create changelog files for both new directories
6. Delete `warp-tools/3rd-party-tools/aou_mitochondrial_annotate_coverage/`
7. Update WDL task docker references in warp repo
8. Build and smoke-test both images
9. Remove `warp/all_of_us/mitochondria/mtSwirl_refactor/coverage_db/Dockerfile`
   and `warp/all_of_us/mitochondria/mtSwirl_refactor/Terra/Dockerfile` (source of truth is now warp-tools)
