# Imputation GLIMPSE

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/imputation-glimpse:1.0.0-2a1a895-1777494633`

- __What is this image:__ This image is built using a local two-stage build script ([docker_build.sh](docker_build.sh)) for running GLIMPSE2 in imputation pipelines.
- __What is GLIMPSE:__ GLIMPSE is a set of tools for phasing and imputation of low-coverage sequencing datasets. See [here](https://odelaneau.github.io/GLIMPSE/) for more information.

## What's included in this image

This image contains:
- **GLIMPSE2 tools**: `GLIMPSE2_chunk`, `GLIMPSE2_split_reference`, `GLIMPSE2_phase`, `GLIMPSE2_ligate`, `GLIMPSE2_concordance` ([based on GLIMPSE repository](https://github.com/odelaneau/GLIMPSE#features))
- **bcftools**: For VCF/BCF file manipulation
- **Picard**: For sequence dictionary updates
- **Google Cloud SDK**: For cloud storage operations
- **htslib**: For high-throughput sequencing data manipulation

## Versioning

The Imputation GLIMPSE image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/imputation-glimpse:<image-version>-<short-git-commit-hash>-<unix-timestamp>`

We keep track of all past versions in [docker_versions.tsv](docker_versions.tsv) with the last image listed being the currently used version in WARP.

_Note: The commit hash comes from the [GLIMPSE GitHub repository](https://github.com/odelaneau/GLIMPSE)._

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/imputation-glimpse:<tag-placeholder>
$ docker inspect us.gcr.io/broad-gotc-prod/imputation-glimpse:<tag-placeholder>
```

## Building the Image

This image uses a **two-stage build process** implemented in the local [`docker_build.sh`](docker_build.sh) script:

1. **Base image stage**: Clones the [GLIMPSE repository](https://github.com/odelaneau/GLIMPSE) and builds the base image using their Dockerfile, which includes GLIMPSE2 tools and htslib
2. **Extension image stage**: Builds on top of the base image using the local [Dockerfile](Dockerfile) in this directory, adding:
   - bcftools (for VCF/BCF manipulation)
   - Picard (for sequence dictionary updates)
   - Google Cloud SDK (for cloud storage operations)

The build is automated via GitHub Actions. See [`.github/workflows/build-glimpse.yml`](../../.github/workflows/build-glimpse.yml) for the automation workflow.

### Building Locally

To manually build the image locally, use the provided `docker_build.sh` script:

```bash
# Navigate to the glimpse directory
cd 3rd-party-tools/glimpse

# Build with auto-generated tag (no push to GCR)
./docker_build.sh --no-push

# Build with custom tag
./docker_build.sh -t "my-custom-tag" --no-push

# Build from specific GLIMPSE commit
./docker_build.sh -c 2a1a895 --no-push

# Build from different GLIMPSE branch
./docker_build.sh -b develop --no-push

# Build from different GLIMPSE repository
./docker_build.sh -r "https://github.com/user/GLIMPSE.git" -b feature --no-push

# Skip confirmation prompts (useful for scripts)
./docker_build.sh -y --no-push

# Combine options
./docker_build.sh -c 2a1a895 -t "1.0.0-2a1a895-custom" -y --no-push

# Build and push to GCR (requires Docker to be authenticated with GCR)
./docker_build.sh -t "1.0.0-test"

# For help and all options
./docker_build.sh --help
```

**Note:** The script requires `docker`, `gcloud`, and `git` to be installed. It will:
1. Clone the GLIMPSE repository to a temporary directory
2. Build the GLIMPSE base image from their Dockerfile
3. Build the extension image using the local Dockerfile
4. Clean up temporary files
5. Optionally push to GCR (requires appropriate permissions)
