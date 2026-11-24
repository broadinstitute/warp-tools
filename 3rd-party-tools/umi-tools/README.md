# RNA Seq UMI-tools

## Quick reference

Copy and paste to pull this image

#### `us.gcr.io/broad-gotc-prod/umi_tools:1.0.0-1.1.1-1638821470`

- __What is this image:__ This image is a lightweight debian based image for running UMI-tools within our RNA sequencing pipeline.
- __What is fgbio:__ UMI-tools is set of tools for dealing with Unique Molecular Identifiers (UMIS)/Random Molecular Tags (RMTS) and single-cell RNA-seq barcodes. See [here](https://github.com/CGATOxford/UMI-tools) for more information.
- __How to see tool version used in image:__ Please see below.

## Versioning

fgbio uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/umi_tools:<image-version>-<umi_tools-version>-<unix-timestamp>` 

WARP pins the version it runs in its WDLs; see [Pinning an image in WARP](../../BUILDING.md#pinning-an-image-in-warp).

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/umi_tools:1.0.0-1.1.1-1638821470
$ docker inspect us.gcr.io/broad-gotc-prod/umi_tools:1.0.0-1.1.1-1638821470
```

## Usage


```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/umi_tools:1.0.0-1.1.1-1638821470 umi_tools
```