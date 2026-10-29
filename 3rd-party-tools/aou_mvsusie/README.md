# cis mvSuSiE image

One Docker image for every task in `CisMvSusie.wdl`.

## Sources and provenance

| Component | Source | How it is included |
| --- | --- | --- |
| WARP container packaging | [broadinstitute/warp-tools](https://github.com/broadinstitute/warp-tools/tree/develop/3rd-party-tools/aou_mvsusie) | This Dockerfile and README |
| Analysis scripts and original container/environment design | [AoU-Multiomics-Analysis/mvSusieR](https://github.com/AoU-Multiomics-Analysis/mvSusieR) | Scripts are vendored in [scripts/](scripts/); the container adapts the original micromamba-based design to R installations without micromamba |
| `mvsusieR` R package | [stephenslab/mvsusieR](https://github.com/stephenslab/mvsusieR) | Installed from GitHub using `MVSUSIER_REF` (default: `HEAD`) |
| `susieR` R package | [stephenslab/susieR](https://github.com/stephenslab/susieR) | Installed from GitHub using `SUSIER_REF` (default: `HEAD`), or from CRAN when `SUSIER_VERSION` is set |
| R base image | [rocker/r-ver](https://github.com/rocker-org/rocker-versioned2) | `rocker/r-ver:4.4.1` |

The AoU analysis repository and the Stephens lab R package repositories are distinct
sources. The original environment specified `r-mvsusier=0.3.0` and `r-susier>=0.15`;
it did not pin an exact `susieR` version. This image's GitHub defaults do not guarantee
the same resolved package versions as that environment.

[get_scripts.sh](get_scripts.sh) copies the upstream `scripts/` directory, overwriting
matching local files, and prints the downloaded commit SHA. The Docker build itself
copies the local scripts rather than downloading them. Their exact upstream revision
is not currently recorded in this directory; retain the printed SHA when refreshing
scripts and document any local changes. Missing required scripts cause the build to fail.

The image includes OCI source labels for this packaging repository and custom
`org.broadinstitute.warp-tools.*.source` labels for the script and R package sources.
These URLs identify origins, not pinned revisions or a license grant; consult each
upstream repository for its license and attribution requirements.

## Build and push (for Terra on Google Cloud)

From this directory, with Docker running and GCR authentication configured:

```bash
docker build --platform linux/amd64 \
  -t us.gcr.io/broad-gotc-prod/aou-mvsusie:1.0.0-20261029 .
docker push us.gcr.io/broad-gotc-prod/aou-mvsusie:1.0.0-20261029
```

Building on Apple Silicon uses amd64 emulation and can be slow. Use a new versioned
tag for each release, and pin the published `@sha256` digest in the consuming WARP
pipeline. See [BUILDING.md](../../BUILDING.md) for registry access and versioning.

The installed package versions are saved inside the image at `/opt/mvsusie/VERSIONS.txt`.

## Versions to check

Pin compatible upstream commits before production use; `HEAD` can change between builds:

    docker buildx build --platform linux/amd64 \
  --build-arg MVSUSIER_REF=<commit-sha> --build-arg SUSIER_REF=<commit-sha> \
  -t us.gcr.io/broad-gotc-prod/aou-mvsusie:<versioned-tag> --push .

`SUSIER_VERSION` overrides `SUSIER_REF` and must identify an available CRAN release.
The build checks `susieR >= 0.15.54`; the scripts must also be tested against the
selected package APIs.

The collaborator's `session_info.txt` output (from `FitMvSusie`) lists the versions they used.
