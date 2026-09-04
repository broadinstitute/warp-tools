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
- **Keep GPU images device-agnostic — don't hardcode a device.** scvi-tools defaults to `accelerator="auto"`, using a GPU when `torch.cuda.is_available()` and CPU otherwise. This lets ONE image back both a GPU and a CPU-only WDL task; the WDL `runtime` (not the container) decides whether a GPU is attached. Adding explicit `accelerator=`/`devices=` risks changing device-*count* behavior on the multi-GPU path — leave it default unless you mean to. (The warp side needs two tasks — GPU and CPU — for this; see the warp `AGENTS.md` GPU note.)
- **Prefer a container entry point over a big inline-Python WDL heredoc.** scANVI's post-preprocessing step lives in `label_transfer_from_preprocessed.py` (a CLI over the same imported functions), so the WDL task is a one-line `python3 …` call. This keeps the command and the functions it calls in **one image** — no WDL↔image kwarg skew (an image predating a new arg like `max_epochs`/`batch_size` otherwise raises `unexpected keyword argument`) — and lets the WDL split cheaply into GPU/CPU task variants. Still honor the importable-function contract above.

## Pipeline verification & CI (in the warp consuming repo)

How this image's outputs get tested lives in the [warp](https://github.com/broadinstitute/warp) repo, but these cross-repo lessons are easy to trip over:

- **Tolerant verification for stochastic models.** scvi-scanvi training is non-deterministic, so the scANVI pipeline does **not** compare outputs to truth by exact equality. Its `CompareScanviH5ad` task checks structure + distribution: cell counts match, the predicted-label vocabulary is a subset of truth's, and per-cell-type proportions correlate with truth above a threshold. Use this pattern for any stochastic/ML image.
- **Keep a pipeline's verification task out of the shared `verification/VerifyTasks.wdl`.** That shared file is in **~27** pipelines' GitHub Actions `paths:` filters, so editing it triggers **every** pipeline's (Terra-based) test suite — slow, costly, noisy. Put the compare task in the pipeline's own `verification/Verify<Pipeline>.wdl` instead (scANVI keeps `CompareScanviH5ad` in `VerifyScANVI.wdl`). The same "don't touch in a feature PR" caution applies to other shared trigger files: `tasks/wdl/Utilities.wdl`, `tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl`, `.github/workflows/warp_test_workflow.yml`, and `scripts/firecloud_api/firecloud_api.py` (each in ~24–27 path filters).
- **New `test_<pipeline>.yml` needs an explicit `permissions:` block.** The reusable `warp_test_workflow.yml` job requests `contents: read`, `id-token: write` (Terra/GCP auth), and `actions: write`; the caller must grant at least those (CodeQL also flags a missing block).

### Future work (stashed, separate PR)

- **`firecloud_api.py` robustness.** Terra/Dockstore intermittently return non-JSON 5xx (e.g. an HTML `502`) during submission polling and method-config cleanup; `firecloud_api.py` parses the body unconditionally and crashes with `requests.exceptions.JSONDecodeError`. This surfaces as widespread, transient red CI across unrelated pipelines — infrastructure flakiness, not real breakage. A future PR should detect non-JSON / 5xx responses and retry-with-backoff instead of crashing. Note: `firecloud_api.py` is itself in ~24 pipeline `paths:` filters, so that fix is a CI-infra PR that triggers every pipeline — keep it out of feature PRs.

## Agent-Specific Notes

Record recurring container-build learnings/mistakes here. Keep entries short and **link** to [BUILDING.md](BUILDING.md) or the warp AGENTS.md rather than duplicating them. Before adding a note, check it isn't already covered above or in a linked doc.

- **Dockerfiles `COPY` each script explicitly (no wildcard).** Adding a new `.py` to a tool dir does nothing until you add a matching `COPY <file> .` line to that tool's `Dockerfile` — otherwise it's silently absent from the built image and a WDL calling it fails at **runtime**, not build (scvi-scanvi hit exactly this when `label_transfer_from_preprocessed.py` was added). Check the `COPY` list whenever you add files.
- **Per-branch image tags are reused, so untagged digests pile up.** `build-<tool>.yml` triggers on `pull_request` (to develop/master, path-filtered) and `workflow_dispatch` — **no `push` trigger** — building once per commit and tagging with the branch name (`head_ref`), not per-commit. Each new commit's build **moves** that tag to a new digest, leaving prior builds untagged in GCR. So a PR shows several untagged digests with only the latest tagged — expected, not a double-build. Always pin the consuming pipeline to the `@sha256` digest, never the branch tag.
