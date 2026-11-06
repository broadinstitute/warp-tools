# GATK BCFtools gcloud

## Quick reference

Copy and paste to pull this image (see [docker_versions.tsv](docker_versions.tsv) for the current tag once published)

```bash
docker pull us.gcr.io/broad-gotc-prod/gatk-bcftools-gcloud:<image-version>-<gatk4-version>-<bcftools-version>-<unix-timestamp>
```

- **What is this image:** A lightweight Alpine-based image bundling GATK4, BCFtools, and the Google Cloud CLI (`gcloud`/`gsutil`) in a single execution environment for WDL tasks that need all three.
- **What is `GATK`:** GATK (Genome Analysis Toolkit) is the industry standard for identifying SNPs and indels in germline DNA and RNAseq data. See [here](https://gatk.broadinstitute.org/hc/en-us) for more information.
- **What is `BCFtools`:** A suite of tools for variant calling and manipulating BCFs and VCFs. See [here](https://github.com/samtools/bcftools) for more information.
- **What is `gcloud`:** The Google Cloud CLI, providing `gcloud` and `gsutil` for interacting with Google Cloud Storage and other GCP services from within a task. See [here](https://cloud.google.com/sdk/gcloud) for more information.

## Synopsis

This image combines [gatk](../gatk) and the BCFtools portion of [bcftools-vcftools](../bcftools-vcftools) with the Google Cloud CLI so that a single WDL task can run GATK4, BCFtools, and `gcloud`/`gsutil` commands without needing to localize/delocalize files through a separate task. It's intended for pipeline steps that mix variant-calling/manipulation with direct GCS access (e.g. streaming a file from a bucket mid-task).

## Image contents

- Base image: `eclipse-temurin:8-jdk-alpine`
- `gatk4` `4.2.6.1`
- `bcftools` `1.24`
- `htslib` `1.24`
- Google Cloud CLI (`gcloud`, `gsutil`)

## Versioning

This image follows the WARP-Tools tag convention
`us.gcr.io/broad-gotc-prod/gatk-bcftools-gcloud:<image-version>-<gatk4-version>-<bcftools-version>-<unix-timestamp>`
(see [BUILDING.md → Versioning strategy](../../BUILDING.md#versioning-strategy)).
Past versions are tracked in [docker_versions.tsv](docker_versions.tsv); the last
line is the version currently used by WARP. Inspect an image with:

```bash
docker inspect us.gcr.io/broad-gotc-prod/gatk-bcftools-gcloud:<tag>
```

## Usage

Run an interactive shell in the image:

```bash
docker run --rm -it us.gcr.io/broad-gotc-prod/gatk-bcftools-gcloud:<tag> bash
```

### GATK4

```bash
$ docker run --rm -it \
    -v /gatk-files:/gatk-files \
    us.gcr.io/broad-gotc-prod/gatk-bcftools-gcloud:<tag> gatk --java-options "-Xms2000m -Xmx2500m" \
    PrintReads \
    -I /gatk-files/input_bam \
    --interval-padding 500 \
    -L /gatk-files/interval_list \
    -O /gatk-files/local.sharded.bam
```

### BCFtools

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/gatk-bcftools-gcloud:<tag> bcftools
```

### gcloud / gsutil

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/gatk-bcftools-gcloud:<tag> gsutil ls gs://<bucket>/
```

## Task descriptions

- **Pipeline:** TBD — fill in once a WARP pipeline task is wired to this image.

## Troubleshooting

To run the image independently of a WDL for testing, start a shell explicitly
(`docker run -it --rm <image-url> bash`). See
[BUILDING.md → Troubleshooting and running standalone](../../BUILDING.md#troubleshooting-and-running-standalone).
