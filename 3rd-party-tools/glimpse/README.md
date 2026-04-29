# Imputation GLIMPSE

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/imputation-glimpse:1.0.0-2a1a895-1777494633`

- __What is this image:__ This image is built using the [palantir-workflows build script](https://github.com/broadinstitute/palantir-workflows/blob/main/GlimpseImputationPipeline/glimpse_docker/build_base_and_extension_docker.sh) for running GLIMPSE2 in imputation pipelines.
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
$ docker pull us.gcr.io/broad-gotc-prod/imputation-glimpse:1.0.0-2a1a895-1777494633
$ docker inspect us.gcr.io/broad-gotc-prod/imputation-glimpse:1.0.0-2a1a895-1777494633
```

## Building the Image

This image is built using the [build script in palantir-workflows repository](https://github.com/broadinstitute/palantir-workflows/blob/main/GlimpseImputationPipeline/glimpse_docker/build_base_and_extension_docker.sh), 
which contains a two-stage build:

1. **Base image**: Built from [this Dockerfile](https://github.com/odelaneau/GLIMPSE/blob/master/Dockerfile) in GLIMPSE repository which includes GLIMPSE2 tools and htslib
2. **Extension image**: Built from [this Dockerfile](https://github.com/broadinstitute/palantir-workflows/blob/main/GlimpseImputationPipeline/glimpse_docker/Dockerfile) in palantir-workflows repository and adds bcftools, Picard, and Google Cloud SDK on top of the above base image

The build is automated via GitHub Actions. See [`.github/workflows/build-glimpse.yml`](../../.github/workflows/build-glimpse.yml) for the automation workflow.

To manually build the image, you can use the palantir build script:

```bash
# Clone the palantir-workflows repository
git clone https://github.com/broadinstitute/palantir-workflows.git
cd palantir-workflows/GlimpseImputationPipeline/glimpse_docker

# Run the build script
bash build_base_and_extension_docker.sh \
  -r "https://github.com/odelaneau/GLIMPSE.git" \
  -b master \
  -t "us.gcr.io/broad-gotc-prod/imputation-glimpse:<your-tag>"
```
