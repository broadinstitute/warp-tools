# UPSTOOLS

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/upstools:1.0.0-2023.03.03-1701813153`


- __What is this image:__ This image is a custom image for running SnapATAC2.
- __What upstools:__ UPSTOOLS is a software package for paired-tag demultipexing by 3 bp barcode. 
- __How to see tool version used in image:__ Please see below.

## Versioning

UPSTOOLS image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/snapatac2:<image-version>-<upstools-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/upstools:1.0.0-2023.03.03-1701813153
$ docker inspect us.gcr.io/broad-gotc-prod/upstools:1.0.0-2023.03.03-1701813153
```

## Usage

### Upstools 

```bash
$ docker run --rm -it \
    -v /files:/files \
    us.gcr.io/broad-gotc-prod/upstools:1.0.0-2023.03.03-1701813153 bash
```
Upstools is a c++ module. To handle the 1st round barcode (pre-index) for cellranger-atac pipeline, we extracted the barcode from the beginning of R2 fastq sequence and append it to read name, which can be kept during the process of alignment.