# Imputation Beagle

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/imputation-beagle:2.0.0-d820c4e-1763720636`

- __What is this image:__ This image is a lightweight alpine-based image for running Beagle in the [ImputationBeagle pipeline](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/ImputationBeagleTasks.wdl).
- __What is Beagle:__ Beagle is a software package for phasing genotypes and imputing ungenotyped markers. Beagle version 5.4 has improved memory and computational efficiency when analyzing large sequence data sets. See [here](https://faculty.washington.edu/browning/beagle/beagle.html) for more information.
- __How to see Beagle version used in image:__ Please see below.

## Versioning

The Imputation Beagle image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/imputation-beagle:<image-version>-<short-git-commit-hash>-<unix-timestamp>`

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

_Note: The commit hash comes from GitHub repo [tmp-sharing/imp-server](https://github.com/tmp-sharing/imp-server/tree/master)._

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/imputation-beagle:2.0.0-d820c4e-1763720636
$ docker inspect us.gcr.io/broad-gotc-prod/imputation-beagle:2.0.0-d820c4e-1763720636
```

## Usage

### Display default menu

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/imputation-beagle:2.0.0-d820c4e-1763720636 java -jar /usr/gitc/beagle.d820c4e.jar
```