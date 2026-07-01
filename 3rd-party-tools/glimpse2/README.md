# Imputation GLIMPSE2

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/imputation-glimpse2:1.1.0-c276764-1782839248`

- __What is this image:__ This image is built using a local two-stage build script ([docker_build.sh](docker_build.sh)) for running GLIMPSE2 in imputation pipelines.
- __What is GLIMPSE2:__ GLIMPSE2 is a set of tools for phasing and imputation of low-coverage sequencing datasets. See [here](https://odelaneau.github.io/GLIMPSE/) for more information.

## What's included in this image

This image contains:
- **GLIMPSE2 tools**: `GLIMPSE2_chunk`, `GLIMPSE2_split_reference`, `GLIMPSE2_phase`, `GLIMPSE2_ligate`, `GLIMPSE2_concordance` ([based on GLIMPSE repository](https://github.com/odelaneau/GLIMPSE#features))
- **bcftools**: For VCF/BCF file manipulation
- **Picard**: For sequence dictionary updates
- **Google Cloud SDK**: For cloud storage operations
- **htslib**: For high-throughput sequencing data manipulation

## Versioning

The Imputation GLIMPSE2 image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/imputation-glimpse2:<image-version>-<short-git-commit-hash>-<unix-timestamp>`

We keep track of all past versions in [docker_versions.tsv](docker_versions.tsv) with the last image listed being the currently used version in WARP.

_Note: The commit hash comes from the [GLIMPSE GitHub repository](https://github.com/odelaneau/GLIMPSE)._

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/imputation-glimpse2:<tag-placeholder>
$ docker inspect us.gcr.io/broad-gotc-prod/imputation-glimpse2:<tag-placeholder>
```

## Building the Image

This image uses a **two-stage build process** implemented in the local [`docker_build.sh`](docker_build.sh) script:

1. **Base image stage**: Clones the [GLIMPSE repository](https://github.com/odelaneau/GLIMPSE) and builds the base image using their Dockerfile, which includes GLIMPSE2 tools and htslib
2. **Extension image stage**: Builds on top of the base image using the local [Dockerfile](Dockerfile) in this directory, adding:
   - bcftools (for VCF/BCF manipulation)
   - Picard (for sequence dictionary updates)
   - Google Cloud SDK (for cloud storage operations)

The build is automated via GitHub Actions. See [`.github/workflows/build-glimpse2.yml`](../../.github/workflows/build-glimpse2.yml) for the automation workflow.

### Building Locally

To manually build the image locally, use the provided `docker_build.sh` script:

```bash
# Navigate to the glimpse directory
cd 3rd-party-tools/glimpse2

# Basic usage - build locally without pushing (commit required)
./docker_build.sh -c 2a1a895 --no-push

# Build, record tag and push to GCR (requires Docker authentication)
./docker_build.sh -c 2a1a895 -t "1.0.0-test" --record-tag -y

# For help and all options
./docker_build.sh --help
```

**Common options:**
- `-c, --commit <hash>` - GLIMPSE commit hash (required)
- `-r, --repo <url>` - GLIMPSE repository URL (default: official GLIMPSE repo)
- `-t, --tag <tag>` - Custom image tag (default: auto-generated with format `<image-version>-<short-git-commit-hash>-<unix-timestamp>`)
- `--no-push` - Build locally without pushing to GCR
- `--record-tag` - Save the tag to docker_versions.tsv (useful for CI/CD)
- `-y, --yes` - Skip confirmation prompts

**Note:** The script requires `docker` and `git` to be installed. It will:
1. Clone the GLIMPSE repository to a temporary directory
2. Build the GLIMPSE2 base image from their Dockerfile
3. Build the extension image using the local Dockerfile
4. Clean up temporary files
5. Optionally push to GCR (if `--no-push` is not specified)
6. Optionally record the image tag to `docker_versions.tsv` (if `--record-tag` is specified)

**The `--record-tag` flag** appends the final image tag to `docker_versions.tsv`. By default, local builds do NOT 
record to `docker_versions.tsv` unless you explicitly use `--record-tag`. GHA does not use `--record-tag` since we 
want to control which images are recorded in `docker_versions.tsv` and avoid recording test builds or duplicate tags.
