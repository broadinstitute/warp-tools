# Dropseq Picard

## Quick reference

Copy and paste to pull this image

#### `docker pull <IMAGE_tag>`


- __What is this image:__ This image is a custom image for running Dropseq tools with SingleCellRnaSeqMetricsCollector.java, a tool for per barcode metrics for single-cell RNA datasets.
- __What Dropseq:__ Dropseq is a software package for , [more info](). 
- __How to see tool version used in image:__ Please see below.

## Versioning

Dropseq-Picard image uses the following convention for verisoning:

#### `us.gcr.io/broad-gotc-prod/Dropseq-picard:<image-version>-<Picard-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull <docker image>
$ docker inspect <docer image>
```

## Usage

### SnapATAC2 

```bash
$ docker run --rm -it \
    -v /files:/files \
    <docker image bash
```
Dropseq-picard is a java module. To learn more about it, see the [Dropseq code](https://github.com/broadinstitute/Drop-seq/blob/master/src/java/org/broadinstitute/dropseqrna/barnyard/SingleCellRnaSeqMetricsCollector.java). 