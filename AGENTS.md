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

Images publish to **GCR only** (`us.gcr.io/broad-gotc-prod/<image>`); CI builds on a PR touching an image's directory; then pin the resulting `@sha256` digest in the consuming WARP pipeline. Commands and the full flow are in [BUILDING.md](BUILDING.md#how-images-are-consumed-by-warp-digest-pinning).

> There is **no Quay** publishing. CI pushes only to GCR. Some `docker_build.sh` scripts still contain `quay.io/...` tag/push lines — these are **stale and nonfunctional** (the newer images were never pushed to Quay); treat them as dead code to be removed, not a live target.

## Keep scripts importable

WARP pipelines frequently import specific functions from a baked-in script rather than calling its `main()` (e.g. scANVI runs `from multiome_label_transfer import run_multi_model, run_gex_only_model, …`). The imported function API — not `main()` — is the production contract. Preserve signatures/behavior, factor shared internals into a helper when adding a variant (scvi-scanvi: `_concat_and_train_models` is shared by `run_multi_model` and `run_gex_only_model`), and document hidden preconditions in docstrings. See [BUILDING.md → Adding or updating an image](BUILDING.md#adding-or-updating-an-image).

## Local build & GPU testing

You can build and run these images locally — including GPU images — without CI.

- **Sandbox note.** In a Flatpak/dev-container sandbox the host tools live under `/run/host` and are not on `PATH`. Run the host engine via `flatpak-spawn --host podman …` (direct exec of `/run/host/usr/bin/podman` fails on host libs). **Check for an existing host engine before installing one.**
- **Build:** `flatpak-spawn --host podman build -t <tool>:local 3rd-party-tools/<tool>`.
- **GPU passthrough:** the host NVIDIA Container Toolkit (`nvidia-ctk`) provides a CDI spec at `/etc/cdi/nvidia.yaml`; attach the GPU with `--device nvidia.com/gpu=all`. Smoke-test with `podman run --rm --device nvidia.com/gpu=all <img> nvidia-smi` and an in-container `python -c "import torch; print(torch.cuda.is_available())"`.
- **Test edited code without rebuilding** by volume-mounting the changed script over the image copy: `-v "$PWD/3rd-party-tools/<tool>/script.py:/usr/local/script.py:Z"`.
- **Verify a published image's contents** by importing the module and checking expected symbols, e.g. `python -c "import multiome_label_transfer as m; print(hasattr(m,'run_gex_only_model'))"`.
- **Pull auth (GCR):** `gcloud auth print-access-token | podman login -u oauth2accesstoken --password-stdin us.gcr.io`.

## Image-specific notes: scvi-scanvi (GPU)

- `scvi-tools` and `snapatac2` are pinned via build `ARG`s — **keep them pinned**; changing `scvi-tools` changes model behavior for the consuming pipeline.
- The gencode annotation is baked in; `snap.genome.hg38` is downloaded at runtime (the container needs network during the ATAC gene-activity step).
- SCVI/SCANVI training is **stochastic** — outputs are not bit-reproducible, so the consuming pipeline verifies **tolerantly** (warp `VerifyScANVI` / `CompareScanviH5ad`), not by exact match. See the image's [README](3rd-party-tools/scvi-scanvi/README.md).

## Agent-Specific Notes

Record recurring container-build learnings/mistakes here. Keep entries short and **link** to [BUILDING.md](BUILDING.md) or the warp AGENTS.md rather than duplicating them. Before adding a note, check it isn't already covered above or in a linked doc.
