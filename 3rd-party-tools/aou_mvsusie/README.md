# cis mvSuSiE image

One Docker image for every task in `CisMvSusie.wdl`.

## What is in this zip

- `Dockerfile`: R 4.4.1 (rocker/r-ver), tabix/bgzip, the R packages, mashr, mvsusieR, and the scripts.
- `scripts/`: five scripts shared in chat (`prepare_trans_window.R`, `prepare_window.R`,
  `fit_window.R`, `trans_window_io.R`, `trans_window_preprocess.R`).
- `get_scripts.sh`: copies the whole `scripts/` folder from
  https://github.com/AoU-Multiomics-Analysis/mvSusieR into `scripts/`.

## Six scripts are not in this zip

They were not shared, so the Dockerfile stops with a clear message if they are missing:
`trans_window_cli.R`, `trans_window_logging.R`, `trans_window_model.R`, `trans_window_prior.R`,
`summarize_window.R`, `plot_window_mvsusie.R`.
Run `./get_scripts.sh` (needs git and internet) or copy them into `scripts/` by hand.

## Build and push (for Terra on Google Cloud)

Run `./get_scripts.sh` first, then build natively on Google Cloud:

    gcloud builds submit --timeout=60m --machine-type=e2-highcpu-8 \
      --tag <region>-docker.pkg.dev/<project>/<repo>/cis-mvsusie:0.1 .

Terra must be able to pull the image. Either make the Artifact Registry repository public, or
grant `roles/artifactregistry.reader` to the Terra proxy group or your pet service account.
Use a version tag, not `latest`. Then give the image name to the WDL input `docker_image`.

Building locally also works but is slow on Apple Silicon:

    docker buildx build --platform linux/amd64 -t <registry>/<name>:0.1 --push .

The installed package versions are saved inside the image at `/opt/mvsusie/VERSIONS.txt`.

## Versions to check

`trans_window_model.R` calls mvsusieR, whose API changed in its 2026 refactor. If the fit step
fails with unused or unknown arguments, rebuild with a pinned version:

    docker buildx build --platform linux/amd64 \
      --build-arg MVSUSIER_REF=<commit-or-tag> --build-arg SUSIER_VERSION=<x.y.z> \
      -t ghcr.io/<you>/cis-mvsusie:0.1 --push .

The collaborator's `session_info.txt` output (from `FitMvSusie`) lists the versions they used.
