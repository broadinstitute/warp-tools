# WARP-Tools

WARP-Tools hosts the Docker images — and the in-house scripts and tools they contain — used by the [WARP](https://github.com/broadinstitute/warp) cloud-optimized pipelines. Each image is an execution environment for one or more processing steps in a WARP WDL workflow.

This README is a map of the repository. For build instructions and the Dockerfile style guide see **[BUILDING.md](BUILDING.md)**; for agent/automation guidance see **[AGENTS.md](AGENTS.md)**.

## Repository structure

| Path | What it is |
| --- | --- |
| `tools/` | The in-house **`warp-tools`** image — a single image bundling Broad-authored tools (`fastqpreprocessing` and `TagSort` in C++, plus helper `scripts/`). It has its own `Dockerfile`, `docker_build.sh`, `docker_versions.tsv`, and `warp-tools.changelog.md` at the top level. |
| `3rd-party-tools/` | One subdirectory per third-party image (~57). Each contains a `Dockerfile`, a `docker_build.sh`, a `docker_versions.tsv`, a `README.md`, and any baked-in scripts/assets. |
| `.github/workflows/` | CI. One `build-<tool>.yml` per image (builds & pushes to GCR), plus `trivy.yml` for vulnerability scanning. |

Every image directory carries its own `README.md` describing the tool and its usage — see the [tool README template](TOOL_README_TEMPLATE.md) for the standard structure.

## Where images live

All images are published to **Google Container Registry (GCR)**, the single source of truth:

```text
us.gcr.io/broad-gotc-prod/<image>
```

There is no Quay or other mirror.

## Versioning strategy

Images use **semantic tags** (never rolling tags like `latest`/`master`):

```text
us.gcr.io/broad-gotc-prod/<image>:<image-version>-<tool-version>-<unix-timestamp>
```

- **`image-version`** — `major.minor.patch` of the image; bumped for image-level changes (base OS, system packages, baked-in scripts).
- **`tool-version`** — version of the wrapped tool (e.g. `samtools` `1.19`).
- **`unix-timestamp`** — keeps every tag unique and avoids Cromwell image-cache collisions.

Each image's `docker_versions.tsv` records every published tag; the **last line is the version currently used by WARP**. WARP pipelines pin images by **immutable `@sha256` digest**, so a pipeline always runs an exact image — see [How images are consumed by WARP](BUILDING.md#how-images-are-consumed-by-warp-digest-pinning).

## Building images

Images are built by CI: open a pull request that touches an image's directory and its `build-<tool>.yml` workflow builds and pushes it to GCR (tagged by branch). A `docker_build.sh` next to each Dockerfile supports local builds. Full details, plus the Dockerfile style guide and how to add a new image, are in **[BUILDING.md](BUILDING.md)**.

## Citing WARP and WARP-Tools

When citing WARP and WARP-Tools, please use the following:

Degatano, K.; Awdeh, A.; Dingman, W.; Grant, G.; Khajouei, F.; Kiernan, E.; Konwar, K.; Mathews, K.; Palis, K.; Petrillo, N.; Van der Auwera, G.; Wang, C.; Way, J.; Pipelines, W. WDL Analysis Research Pipelines: Cloud-Optimized Workflows for Biological Data Processing and Reproducible Analysis. Preprints 2024, 2024012131. <https://doi.org/10.20944/preprints202401.2131.v1>
