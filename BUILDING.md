# Building and Maintaining WARP-Tools Images

How to build, version, publish, and add Docker images in this repository, plus the Dockerfile style guide. For a high-level overview of the repo see [README.md](README.md); for agent/automation-specific guidance see [AGENTS.md](AGENTS.md).

## Table of Contents

- [How images are built and published](#how-images-are-built-and-published)
- [Versioning strategy](#versioning-strategy)
- [Pinning an image in WARP](#pinning-an-image-in-warp)
- [Adding or updating an image](#adding-or-updating-an-image)
- [Docker style guide](#docker-style-guide)
- [Troubleshooting and running standalone](#troubleshooting-and-running-standalone)

## How images are built and published

All images publish to **Google Container Registry (GCR)** under `us.gcr.io/broad-gotc-prod/<image>`. There is **no Quay (or other) mirror** — GCR is the single source of truth.

There are two build paths:

1. **CI (canonical).** Each image has a `.github/workflows/build-<tool>.yml` workflow that triggers on a pull request touching that image's directory (e.g. `paths: ['3rd-party-tools/<tool>/**']`), and can also be run manually via **workflow_dispatch** (optional `image_tag` input). Changes to the image's `README`, `docker_build.sh` or `docker_versions.tsv` alone start no build, since CI doesn't use them (except glimpse2 and sv-imputation-rust-tools, whose CI runs `docker_build.sh`). It builds the Dockerfile and pushes to `us.gcr.io/broad-gotc-prod/<tool>`, **tagged by the branch name** (`github.head_ref`). Authentication uses the `GCR_CI_KEY` repository secret — no local credentials are needed. **Opening a PR is the normal way to build an image.**
2. **Local (`docker_build.sh`).** Builds the image on your machine for testing. Its `DOCKER_IMAGE_VERSION`, tool versions and `IMAGE_TAG` line also set the version tag WARP gives the CI build (next sections). Most scripts also still push to GCR; don't use that part, and remove it when you next change the image ([Tidying an image](#tidying-an-image)).

> Image vulnerability scanning runs on every PR via [trivy](https://github.com/aquasecurity/trivy); when you add a new image, register it in [.github/workflows/trivy.yml](.github/workflows/trivy.yml).

## Versioning strategy

Images use **semantic tags**, never rolling tags like `latest` or `master` (rolling tags make it impossible to trace which image content a pipeline actually ran):

```
us.gcr.io/broad-gotc-prod/<image>:<image-version>[-<tool-version>]-<unix-timestamp>
```

- **`image-version`** — `major.minor.patch` of the image itself; bump it when the image changes for reasons unrelated to the wrapped tool (base OS, system packages, baked-in scripts).
- **`tool-version`** — the version of the wrapped tool (e.g. `samtools` `1.19`), so consumers can identify it without inspecting the image. Left out for images with no single wrapped tool, such as `warp-tools`.
- **`unix-timestamp`** — guarantees a unique tag and avoids Cromwell image-cache collisions.

Version tags are added to a PR's CI build by WARP's `image_pins.py` (next section). While the PR is open, each run moves the tag to the build it pins; once the PR merges, the tag is fixed, and a fix is a new image version with a new tag. `DOCKER_IMAGE_VERSION` and the tool versions in each image's `docker_build.sh` are the `<image-version>[-<tool-version>]` part of the next version tag. `image_pins.py` refuses a `DOCKER_IMAGE_VERSION` below a version already in the registry, and a tool version that differs from the Dockerfile's `ARG` default, because CI builds the Dockerfile without build args.

## Pinning an image in WARP

WARP pins each image as `<image>:<version-tag>@sha256:<digest>`: the digest fixes what runs and the tag says which image version it is. Never pin a branch tag, `latest`, a bare tag or a bare digest. A change takes four steps; the WARP side is in WARP's [Updating a warp-tools image](https://broadinstitute.github.io/warp/docs/contribution/contribute_to_warp/update_warp_tools_image) page.

1. **Change the image here and open a PR to `develop`.** Bump `DOCKER_IMAGE_VERSION` (and the tool version) in `docker_build.sh`, and [tidy the image](#tidying-an-image) in the same PR. CI builds the image and pushes it under the branch name.
2. **From a WARP branch, tag and pin the build:** `python3 scripts/image_pins.py pin <PR number>`. It finds the image, its version and the PR's latest build, gives that build the version tag `<version>-<unix-timestamp>`, and rewrites WARP's pins. Rerun after another build and the same tag moves to it, so a PR leaves one version tag per image. It needs write access to the registry (`gcloud auth`).
3. **Bump, test and review in WARP.** If the image has to change again, push here and repeat step 2 once it has built.
4. **Merge this PR, then the WARP PR.** The merge fixes the version tag where it is.

### Tidying an image

Each image drops what's left of the old publishing process the next time it changes, in the PR that changes it. After pinning, `image_pins.py` lists what's left for the image:

- Make `docker_build.sh` build locally only: no `docker push` and no `--build-arg` (so a local build matches CI's). Keep `DOCKER_IMAGE_VERSION`, the tool versions and the `IMAGE_TAG` line:

  ```bash
  #!/bin/bash
  # Builds <tool> locally. CI builds and pushes the image; WARP's image_pins.py tags that build.
  set -e
  DOCKER_IMAGE_VERSION=1.0.0
  TOOL_VERSION="1.2.3"  # the Dockerfile's ARG TOOL_VERSION default
  TIMESTAMP=$(date +"%s")
  IMAGE_TAG="$DOCKER_IMAGE_VERSION-$TOOL_VERSION-$TIMESTAMP"
  docker build -t "<tool>:$IMAGE_TAG" "$(dirname "$0")"
  ```
- In its `README.md`, replace a copied version tag with a link to the WARP WDL that pins the image.

## Adding or updating an image

1. Create `3rd-party-tools/<tool>/` (third-party) or extend `tools/` (in-house) with a `Dockerfile`, a `docker_build.sh` in the [local-only form](#tidying-an-image), and a `README.md`.
2. Write the `README.md` from the standard [tool README template](TOOL_README_TEMPLATE.md) (synopsis, inputs, outputs, usage, task descriptions).
3. Keep any baked-in scripts **importable and stable**: WARP pipelines often import specific functions from a script in the image (e.g. `from multiome_label_transfer import train_scvi`) rather than calling its `main()`. The imported function API — not `main()` — is the production contract, so preserve signatures/behavior and document hidden preconditions in docstrings. When adding a variant, factor shared internals into a helper to avoid drift.
4. Add a `.github/workflows/build-<tool>.yml` (copy an existing one; set `GCR_PATH: broad-gotc-prod/<tool>` and the `paths` trigger) and register the image in [.github/workflows/trivy.yml](.github/workflows/trivy.yml).
5. Open a PR to build the image, then tag and pin it in WARP as above.

## Docker style guide

Guidelines and best practices for writing Dockerfiles in WARP.

### Small images

Smaller images mean faster upload/download, lower storage cost, and a reduced attack surface. The two easiest wins are a small base image and fewer layers.

**Alpine base.** Prefer an [Alpine](https://alpinelinux.org/) base — it's built for size and security, deletes its package index automatically, and provides [tini](https://github.com/krallin/tini) via APK. Use a Debian base only as a last resort, when a dependency isn't available in APK.

```dockerfile
# OKAY, NOT GREAT - uses debian (must clean apt cache manually)
FROM python:debian
RUN set -eux; \
        apt-get update; \
        apt-get install -y curl bash; \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# GOOD - uses alpine
FROM alpine:3.9
RUN set -eux; \
        apk add --no-cache curl bash
```

**Specify the platform.** Images built on ARM machines (e.g. Apple Silicon) can fail the automated PR test suite. Pin `linux/amd64`:

```dockerfile
FROM --platform="linux/amd64" alpine
```

**Minimal RUN steps.** Each instruction creates a layer; collapse installs into a single `RUN`. (Multi-stage builds help for statically linked binaries, but most WARP images need system-level dependencies, so multi-stage is uncommon here.)

```dockerfile
# GOOD - single RUN
RUN set -eux; \
        apk add --no-cache curl bash wget; \
    wget https://example.com/zip; \
    unzip zip
```

### Publicly accessible

WARP pipelines are designed for public use, so the images should be too:

- **Anybody can pull them** — they are hosted in the public GCR repo `us.gcr.io/broad-gotc-prod`.
- **Anybody can build them** — all functionality is encapsulated in the Dockerfile. Download custom software/dependencies from public links; never copy files from internal Broad infrastructure into an image.

### Semantic tagging

Use the semantic tag described in [Versioning strategy](#versioning-strategy); avoid rolling tags.

### Proper process reaping

In a container, PID 1 is responsible for reaping orphaned/zombie processes, but the default `/bin/sh` won't do it — risking resource leaks. Install [tini](https://github.com/krallin/tini/issues/8) (available via APK) and set it as the entrypoint:

```dockerfile
FROM alpine:3.9
RUN set -eux; \
        apk add --no-cache tini
ENTRYPOINT ["/sbin/tini", "--"]
```

### Build scripts and README

Each image keeps a `docker_build.sh` next to its `Dockerfile` that builds it locally and holds its version ([Tidying an image](#tidying-an-image) shows the form). Every image also has a `README.md` (use the [tool README template](TOOL_README_TEMPLATE.md)).

### Formatting

- `ARG`, `ENV`, `LABEL` in that order
- Always include tool versions in `LABEL`
- Single `RUN` step
- Alphabetize package installs
- Clean up the package-index cache
- Use `;` (not `&&`) for line continuation
- Logically separate steps within a `RUN`
- Four-space indentation
- Short comments describing each step
- `tini` is always the default entrypoint

```dockerfile
# Debian base (packages here are not available in Alpine)
FROM us.gcr.io/broad-dsp-gcr-public/base/python:debian

ARG GIT_HASH=c1cba76e979904eb69c31520a0d7f5be63c72253
ENV TERM=xterm-256color \
    BAMID_URL=https://github.com/Griffan/VerifyBamID/archive \
    TINI_VERSION=v0.19.0
LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        GIT_HASH=${GIT_HASH}
WORKDIR /usr/gitc

RUN set -eux; \
        apt-get update; \
        apt-get install -y autoconf cmake g++ gcc git unzip wget zlib1g-dev; \
# Install tini
    wget https://github.com/krallin/tini/releases/download/$TINI_VERSION/tini -O /sbin/tini; \
    chmod +x /sbin/tini; \
# Clean up cached files
    apt-get clean && rm -rf /var/lib/apt/lists/*

ENTRYPOINT [ "/sbin/tini", "--" ]
```

## Troubleshooting and running standalone

WARP images are designed to run from their WDL pipelines. To run one independently for testing, start a shell explicitly:

```bash
docker run -it --rm <image-url> bash
```

Questions? File a [GitHub issue in WARP](https://github.com/broadinstitute/warp/issues/new).
