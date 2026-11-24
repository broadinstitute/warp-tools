# SNAPATAC2

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/snapatac2:1.0.2-2.2.0-1679678908`


- __What is this image:__ This image is a custom image for running SnapATAC2.
- __What SanpATAC2:__ SnatATAC2 is a software package for , [more info](). 
- __How to see tool version used in image:__ Please see below.

## Versioning

SnapATAC2 image uses the following convention for verisoning:

#### `us.gcr.io/broad-gotc-prod/snapatac2:<image-version>-<snapatac2-version>-<unix-timestamp>` 

WARP pins the version it runs in its WDLs; see [Pinning an image in WARP](../../BUILDING.md#pinning-an-image-in-warp).

You can see more information about the image, including the tool versions, by running the following commands:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/snapatac2:1.0.2-2.2.0-1679678908
$ docker inspect us.gcr.io/broad-gotc-prod/snapatac2:1.0.2-2.2.0-1679678908
```

## Usage

### SnapATAC2 

```bash
$ docker run --rm -it \
    -v /files:/files \
    us.gcr.io/broad-gotc-prod/snapatac2:1.0.2-2.2.0-1679678908 bash
```
Snapatac2 is a python module. To learn more about it, see the [Snapatac2 documentation](https://kzhang.org/SnapATAC2/). 