# SV Imputation Rust Tools

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/sv-imputation-rust-tools:1.0.0-5dc0f19-1784328222`

- __What is this image:__ This image is built using a local build script ([docker_build.sh](docker_build.sh)) for running structural variant (SV) imputation tools in genomics pipelines.
- __What are these tools:__ This image contains three Rust-based tools for SV imputation workflows: `pop-glimpse2`, `paste-vcfs`, and `extract-bubble-PLs`, along with BCFtools for VCF/BCF file manipulation.

## What's included in this image

This image contains:
- **Three Rust binaries** from the [lrma-sv-imputation-utils repository](https://github.com/broadinstitute/lrma-sv-imputation-utils):
  - `pop-glimpse2`: Population-based GLIMPSE2 processing tool
  - `paste-vcfs`: VCF file merging and manipulation tool
  - `extract-bubble-PLs`: Extraction tool for variant bubbles
- **BCFtools**: For VCF/BCF file manipulation and variant calling utilities
- **Google Cloud SDK**: For cloud storage operations
- **Runtime libraries**: Including htslib, bzip2, lzma, and other required dependencies

## Versioning

The SV Imputation Rust Tools image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/sv-imputation-rust-tools:<image-version>-<short-git-commit-hash>-<unix-timestamp>`

We keep track of all past versions in [docker_versions.tsv](docker_versions.tsv) with the last image listed being the currently used version in WARP.

_Note: The commit hash comes from the [lrma-sv-imputation-utils GitHub repository](https://github.com/broadinstitute/lrma-sv-imputation-utils)._

You can see more information about the image, including the tool versions and commit hash, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/sv-imputation-rust-tools:1.0.0-5dc0f19-1784328222
$ docker inspect us.gcr.io/broad-gotc-prod/sv-imputation-rust-tools:1.0.0-5dc0f19-1784328222
```

## Building the Image

This image uses a **two-stage build process** defined in the [Dockerfile](Dockerfile):

1. **Builder stage**: 
   - Compiles BCFtools 1.24 from source with optimized configuration
   - Clones the [lrma-sv-imputation-utils repository](https://github.com/broadinstitute/lrma-sv-imputation-utils) at a specific commit
   - Compiles the three Rust binaries (`pop-glimpse2`, `paste-vcfs`, `extract-bubble-PLs`) in release mode

2. **Runtime stage**: 
   - Creates a minimal Ubuntu 24.04 image with only runtime dependencies
   - Copies the compiled BCFtools binary and plugins
   - Copies the three compiled Rust binaries to the system PATH
   - Installs Google Cloud SDK for cloud operations

The build is automated via GitHub Actions. See [`.github/workflows/build-sv-imputation-rust-tools.yml`](../../.github/workflows/build-sv-imputation-rust-tools.yml) for the automation workflow.

### Building Locally

To manually build the image locally, use the provided `docker_build.sh` script:

```bash
# Navigate to the sv-imputation-rust-tools directory
cd 3rd-party-tools/sv-imputation-rust-tools

# Basic usage - build with latest commit from main (auto-detected)
./docker_build.sh --no-push

# Build with a specific commit hash (must be full 40-character hash)
./docker_build.sh -c 1234567890abcdef1234567890abcdef12345678 --no-push

# Build, record tag and push to GCR (requires Docker authentication)
./docker_build.sh -c 1234567890abcdef1234567890abcdef12345678 -t "1.0.0-test" --record-tag -y

# For help and all options
./docker_build.sh --help
```

**Common options:**
- `-c, --commit <hash>` - Full 40-character commit hash from lrma-sv-imputation-utils (optional, fetches latest from main if not provided)
- `-t, --tag <tag>` - Custom image tag (default: auto-generated with format `<image-version>-<short-git-commit-hash>-<unix-timestamp>`)
- `--no-push` - Build locally without pushing to GCR
- `--record-tag` - Save the tag to docker_versions.tsv (useful for CI/CD)
- `-y, --yes` - Skip confirmation prompts

**Note:** The script requires `docker` and `git` to be installed. It will:
1. Validate the commit hash format (if provided)
2. Fetch the latest commit from the main branch (if no commit is provided)
3. Build the Docker image with the specified commit hash
4. Display the commit hash being used during the build
5. Optionally push to GCR (if `--no-push` is not specified)
6. Optionally record the image tag to `docker_versions.tsv` (if `--record-tag` is specified)

**The `--record-tag` flag** appends the final image tag to `docker_versions.tsv`. By default, local builds do NOT 
record to `docker_versions.tsv` unless you explicitly use `--record-tag`. GHA does not use `--record-tag` since we 
want to control which images are recorded in `docker_versions.tsv` and avoid recording test builds or duplicate tags.
