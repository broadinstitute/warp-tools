# Imputation Beagle

## Quick reference

Copy and paste to pull this image

#### `us-central1-docker.pkg.dev/morgan-fieldeng-gcp/imputation-beagle-development/imputation-beagle:0.0.1-01Mar24.d36-wip-temp-20240301`

- __What is this image:__ This image is a lightweight alpine-based image for running Beagle in the [ImputationBeagle pipeline](../../../../pipelines/broad/arrays/imputation_beagle/ImputationBeagle.wdl).
- __What is Beagle:__ Beagle is a software package for phasing genotypes and imputing ungenotyped markers. Beagle version 5.4 has improved memory and computational efficiency when analyzing large sequence data sets. See [here](https://faculty.washington.edu/browning/beagle/beagle.html) for more information.
- __How to see Beagle version used in image:__ Please see below.

## Versioning

The Imputation Beagle image uses the following convention for versioning:

#### `us-central1-docker.pkg.dev/morgan-fieldeng-gcp/imputation-beagle-development/imputation-beagle:<image-version>-<beagle-version>-<manual-timestamp>`

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us-central1-docker.pkg.dev/morgan-fieldeng-gcp/imputation-beagle-development/imputation-beagle:0.0.1-01Mar24.d36-wip-temp-20240301
$ docker inspect us-central1-docker.pkg.dev/morgan-fieldeng-gcp/imputation-beagle-development/imputation-beagle:0.0.1-01Mar24.d36-wip-temp-20240301
```

## Usage

### Display default menu

```bash
$ docker run --rm -it \
    us-central1-docker.pkg.dev/morgan-fieldeng-gcp/imputation-beagle-development/imputation-beagle:0.0.1-01Mar24.d36-wip-temp-20240301 /usr/gitc/beagle
```