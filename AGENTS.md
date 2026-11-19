# Agent Instructions (warp-tools)

Single entry point for agent-driven and AI-coding work in **warp-tools** — the repository that builds and publishes the Docker images consumed by the [WARP](https://github.com/broadinstitute/warp) pipelines.

This file holds **container-operational** guidance only. It does **not** repeat the human-facing build guide or the WARP pipeline/WDL rules — follow the links instead of restating them.

## Authoritative References

| Topic | Document |
| --- | --- |
| Repo overview, directory structure, where images live, versioning strategy | [README.md](README.md) |
| How to build/publish/add images, digest pinning, Dockerfile style guide | [BUILDING.md](BUILDING.md) (this repo; canonical for build mechanics) |
| Standard structure for a tool's README | [TOOL_README_TEMPLATE.md](TOOL_README_TEMPLATE.md) |
| Pipeline/WDL rules, changelog & **pipeline** versioning, when to package logic into a warp-tools image vs inline WDL | [warp `AGENTS.md`](https://github.com/broadinstitute/warp/blob/develop/AGENTS.md) |

## Build & publish, in one line

Open a PR: CI builds the image into GCR (`us.gcr.io/broad-gotc-prod/<image>`), and WARP's `scripts/image_pins.py pin <PR number>` tags that build and pins it. The flow is in [BUILDING.md](BUILDING.md#pinning-an-image-in-warp).

> **Never push from `docker_build.sh`.** Its `docker push` lines, and version tags copied into READMEs, are leftovers of the old process; remove an image's in the PR that next changes it ([Tidying an image](BUILDING.md#tidying-an-image)). There is no Quay mirror and no `docker_versions.tsv` any more.

## Keep scripts importable

WARP pipelines frequently import specific functions from a baked-in script rather than calling its `main()` (e.g. scANVI's `label_transfer_from_preprocessed.py` runs `from multiome_label_transfer import train_scvi, train_scanvi`). The imported function API — not `main()` — is the production contract. Preserve signatures/behavior, factor shared internals into a helper when adding a variant (scvi-scanvi: `train_scvi` serves both multiome and GEX-only mode), and document hidden preconditions in docstrings. See [BUILDING.md → Adding or updating an image](BUILDING.md#adding-or-updating-an-image).

## Local build & GPU testing

You can build and run these images locally — including GPU images — without CI.

- **Sandbox note.** In a Flatpak/dev-container sandbox the host tools live under `/run/host` and are not on `PATH`. Run the host engine via `flatpak-spawn --host podman …` (direct exec of `/run/host/usr/bin/podman` fails on host libs). **Check for an existing host engine before installing one.**
- **Build:** `flatpak-spawn --host podman build -t <tool>:local 3rd-party-tools/<tool>`.
- **GPU passthrough:** the host NVIDIA Container Toolkit (`nvidia-ctk`) provides a CDI spec at `/etc/cdi/nvidia.yaml`; attach the GPU with `--device nvidia.com/gpu=all`. Smoke-test with `podman run --rm --device nvidia.com/gpu=all <img> nvidia-smi` and an in-container `python -c "import torch; print(torch.cuda.is_available())"`.
- **Test edited code without rebuilding** by volume-mounting the changed script over the image copy: `-v "$PWD/3rd-party-tools/<tool>/script.py:/usr/local/script.py:Z"`.
- **Verify a published image's contents** by importing the module and checking expected symbols, e.g. `python -c "import multiome_label_transfer as m; print(hasattr(m,'train_scanvi'))"`.
- **Pull auth (GCR):** `gcloud auth print-access-token | podman login -u oauth2accesstoken --password-stdin us.gcr.io`.

## Image-specific notes: scvi-scanvi (GPU)

- `scvi-tools` is pinned via a build `ARG` and `snapatac2` is hardcoded — **keep them pinned**; changing `scvi-tools` (currently 1.5.1) or `snapatac2` (currently 2.10.0) changes model behavior for the consuming pipeline.
- No genome annotation is baked in: snapatac2 downloads `snap.genome.hg38` (or mm10/mm39) at runtime, so the container needs network during the ATAC gene-activity step.
- SCVI/SCANVI training is **stochastic** — outputs are not bit-reproducible, so the consuming pipeline verifies **tolerantly** (warp `VerifyScANVI` / `CompareScanviH5ad`), not by exact match. See the image's [README](3rd-party-tools/scvi-scanvi/README.md).
- **Keep GPU images device-agnostic — don't hardcode a device.** scvi-tools defaults to `accelerator="auto"`, using a GPU when `torch.cuda.is_available()` and CPU otherwise. This lets ONE image back both a GPU and a CPU-only WDL task; the WDL `runtime` (not the container) decides whether a GPU is attached. Adding explicit `accelerator=`/`devices=` risks changing device-*count* behavior on the multi-GPU path — leave it default unless you mean to. (The warp side needs two tasks — GPU and CPU — for this; see the warp `AGENTS.md` GPU note.)
- **Prefer a container entry point over a big inline-Python WDL heredoc.** scANVI's post-preprocessing step lives in `label_transfer_from_preprocessed.py` (a CLI over the same imported functions), so the WDL task is a one-line `python3 …` call. This keeps the command and the functions it calls in **one image** — no WDL↔image kwarg skew (an image predating a new arg like `max_epochs`/`batch_size` otherwise raises `unexpected keyword argument`) — and lets the WDL split cheaply into GPU/CPU task variants. Still honor the importable-function contract above.

## Testing an image in its pipeline

Pipeline tests live in WARP; see its [agents/ci-tests.md](https://github.com/broadinstitute/warp/blob/develop/agents/ci-tests.md). For a stochastic image such as scvi-scanvi, the pipeline's `Verify*` task compares outputs tolerantly (scANVI's `CompareScanviH5ad`), never by exact match.

## Agent-Specific Notes

Record recurring container-build learnings/mistakes here. Keep entries short and **link** to [BUILDING.md](BUILDING.md) or the warp AGENTS.md rather than duplicating them. Before adding a note, check it isn't already covered above or in a linked doc.

- **Dockerfiles `COPY` each script explicitly (no wildcard).** Adding a new `.py` to a tool dir does nothing until you add a matching `COPY <file> .` line to that tool's `Dockerfile` — otherwise it's silently absent from the built image and a WDL calling it fails at **runtime**, not build (scvi-scanvi hit exactly this when `label_transfer_from_preprocessed.py` was added). Check the `COPY` list whenever you add files.
- **A PR's builds share its branch tag**, so each new commit's build leaves the previous one untagged in GCR. That's expected, not a double build. Never pin the branch tag; `image_pins.py` pins the build by digest.
- **CI builds the Dockerfile without build args**, so a tool version in `docker_build.sh` must equal the Dockerfile's `ARG` default; `image_pins.py` refuses a mismatch (samtools-star-python and samtools-bwa-mem-2-lisa have one today).
