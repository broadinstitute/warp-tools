# Imputation Beagle

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/imputation-beagle:0.0.2-17Dec24.224-1739980713`

- __What is this image:__ This image is a lightweight alpine-based image for running Beagle in the [ImputationBeagle pipeline](../../../../pipelines/broad/arrays/imputation_beagle/ImputationBeagle.wdl).
- __What is Beagle:__ Beagle is a software package for phasing genotypes and imputing ungenotyped markers. Beagle version 5.4 has improved memory and computational efficiency when analyzing large sequence data sets. See [here](https://faculty.washington.edu/browning/beagle/beagle.html) for more information.
- __How to see Beagle version used in image:__ Please see below.

## Versioning

The Imputation Beagle image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/imputation-beagle:<image-version>-<beagle-version>-<unix-timestamp>`

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/imputation-beagle:0.0.2-17Dec24.224-1739980713
$ docker inspect us.gcr.io/broad-gotc-prod/imputation-beagle:0.0.2-17Dec24.224-1739980713
```

## Usage

### Display default menu

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/imputation-beagle:0.0.2-17Dec24.224-1739980713 java -jar /usr/gitc/beagle.17Dec24.224.jar
```