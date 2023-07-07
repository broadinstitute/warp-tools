# HISAT2-3N

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/hisat3n:1.0.0-2.2.1`

- __What is this image:__ This image is a lightweight alpine-based image for running HISAT2.
- __What is HISAT2-3N:__ HISAT2 is a fast and sensitive alignment program for mapping next-generation sequencing reads to the human genome. See [here](https://github.com/DaehwanKimLab/hisat2) more information.
- __How to see HISAT2-3N version used in image:__ Please see below.

## Versioning

The HISAT2-3N image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/hisat3n:<image-version>-<hisat3n-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/hisat3n:1.0.0-2.2.1
$ docker inspect us.gcr.io/broad-gotc-prod/hisat3n:1.0.0-2.2.1
```

## Usage

### Display default menu

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/hisat3n:1.0.0-2.2.1 hisat3n(?) 
```
