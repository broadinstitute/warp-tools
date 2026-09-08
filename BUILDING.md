# Building and Maintaining WARP-Tools Images

How to build, version, publish, and add Docker images in this repository, plus the Dockerfile style guide. For a high-level overview of the repo see [README.md](README.md); for agent/automation-specific guidance see [AGENTS.md](AGENTS.md).

## Table of Contents

- [How images are built and published](#how-images-are-built-and-published)
- [Versioning strategy](#versioning-strategy)
- [How images are consumed by WARP (digest pinning)](#how-images-are-consumed-by-warp-digest-pinning)
- [Adding or updating an image](#adding-or-updating-an-image)
- [Docker style guide](#docker-style-guide)
- [Troubleshooting and running standalone](#troubleshooting-and-running-standalone)

## How images are built and published

All images publish to **Google Container Registry (GCR)** under `us.gcr.io/broad-gotc-prod/<image>`. There is **no Quay (or other) mirror** — GCR is the single source of truth.

There are two build paths:

1. **CI (canonical).** Each image has a `.github/workflows/build-<tool>.yml` workflow that triggers on a pull request touching that image's directory (e.g. `paths: ['3rd-party-tools/<tool>/**']`), and can also be run manually via **workflow_dispatch** (optional `image_tag` input). It builds the Dockerfile and pushes to `us.gcr.io/broad-gotc-prod/<tool>`, **tagged by the branch name** (`github.head_ref`). Authentication uses the `GCR_CI_KEY` repository secret — no local credentials are needed. **Opening a PR is the normal way to build an image.**
2. **Local/manual (`docker_build.sh`).** A convenience script next to each Dockerfile for building locally. It tags the image per the [versioning strategy](#versioning-strategy), pushes to GCR, and appends the resulting reference to `docker_versions.tsv`. Running it requires a local container engine, `gcloud`, and GCR push access (`gcloud auth configure-docker`).

> Image vulnerability scanning runs on every PR via [trivy](https://github.com/aquasecurity/trivy); when you add a new image, register it in [.github/workflows/trivy.yml](.github/workflows/trivy.yml).

## Versioning strategy

Images use **semantic tags**, never rolling tags like `latest` or `master` (rolling tags make it impossible to trace which image content a pipeline actually ran):

```
us.gcr.io/broad-gotc-prod/<image>:<image-version>-<tool-version>-<unix-timestamp>
```

- **`image-version`** — `major.minor.patch` of the image itself; bump it when the image changes for reasons unrelated to the wrapped tool (base OS, system packages, baked-in scripts).
- **`tool-version`** — the version of the wrapped tool (e.g. `samtools` `1.19`), so consumers can identify it without inspecting the image.
- **`unix-timestamp`** — guarantees a unique tag and avoids Cromwell image-cache collisions.

Each image directory keeps a `docker_versions.tsv` recording every published tag; **the last line is the version currently used by WARP**. The `docker_build.sh` script appends to it automatically. When bumping a manual build, also bump `DOCKER_IMAGE_VERSION` in that image's `docker_build.sh`.

## How images are consumed by WARP (digest pinning)

WARP pipelines pin images by **immutable `@sha256` digest**, not by tag, so a pipeline always runs an exact image. After changing an image you must update the digest in the consuming pipeline or it keeps running the old image:

1. Make the change here and open a PR → CI builds & pushes the branch-tagged image.
2. Resolve the digest for that tag:
   ```bash
   gcloud container images list-tags us.gcr.io/broad-gotc-prod/<image> \
     --filter="tags=<branch-or-tag>" --format="get(digest)"
   ```
   (or read it from the CI "Push image" log).
3. In the [warp](https://github.com/broadinstitute/warp) repo, update the pipeline's `String docker = "us.gcr.io/broad-gotc-prod/<image>@sha256:<digest>"` (and any docs that cite the digest), and validate the WDL.

## Adding or updating an image

1. Create `3rd-party-tools/<tool>/` (third-party) or extend `tools/` (in-house) with a `Dockerfile`, a `docker_build.sh`, a `docker_versions.tsv`, and a `README.md`.
2. Write the `README.md` from the standard [tool README template](TOOL_README_TEMPLATE.md) (synopsis, inputs, outputs, usage, task descriptions).
3. Keep any baked-in scripts **importable and stable**: WARP pipelines often import specific functions from a script in the image (e.g. `from multiome_label_transfer import run_multi_model`) rather than calling its `main()`. The imported function API — not `main()` — is the production contract, so preserve signatures/behavior and document hidden preconditions in docstrings. When adding a variant, factor shared internals into a helper to avoid drift.
4. Add a `.github/workflows/build-<tool>.yml` (copy an existing one; set `GCR_PATH: broad-gotc-prod/<tool>` and the `paths` trigger) and register the image in [.github/workflows/trivy.yml](.github/workflows/trivy.yml).
5. Open a PR to build the image, then pin its digest in the consuming pipeline as above.

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

Each image keeps an easy-to-use `docker_build.sh` next to its `Dockerfile`, with configurable inputs for tool versions, and writes published tags to an accompanying `docker_versions.tsv` automatically. Every image also has a `README.md` (use the [tool README template](TOOL_README_TEMPLATE.md)). See [3rd-party-tools/samtools](3rd-party-tools/samtools) for an example (`docker_build.sh`, `docker_versions.tsv`, `README.md`).

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
