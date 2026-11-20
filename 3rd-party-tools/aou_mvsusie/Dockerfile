# One image for every task in CisMvSusie.wdl.
#
# Terra/Cromwell runs on Google Compute Engine VMs (Intel/AMD, linux/amd64) as root,
# starts the task with its own command (bash), and mounts the task directory at
# /cromwell_root. Nothing here needs an entrypoint or a working directory.
#
# Previous images: us.gcr.io/broad-gotc-prod/aou-mvsusie:1.1.0 and :1.0.0-20261029
#
# Build from this directory and this checkout's scripts/.
# Do not run ./get_scripts.sh first; it overwrites the local script changes.
# On Google Cloud (native amd64, no emulation):
#   gcloud builds submit --timeout=60m --machine-type=e2-highcpu-8 --no-cache \
#     --tag us.gcr.io/broad-gotc-prod/aou-mvsusie:1.1.1 .
# On a Mac (amd64 emulation can be slow on Apple Silicon):
#   docker buildx build --platform linux/amd64 --pull --no-cache \
#     -t us.gcr.io/broad-gotc-prod/aou-mvsusie:1.1.1 --load .
# Rescan the rebuilt image before publishing or using it on Terra.
#   docker push us.gcr.io/broad-gotc-prod/aou-mvsusie:1.1.1
#
# Pinning: the R packages come from the dated Posit package snapshot that rocker/r-ver
# uses for this R version, so CRAN versions are reproducible. mvsusieR and susieR
# are installed from GitHub by default. Set MVSUSIER_REF and SUSIER_REF (or a released
# CRAN SUSIER_VERSION) to match what trans_window_model.R was
# written against. The "R session information" file written by FitMvSusie in the
# collaborator's runs (session_info.txt) lists the exact versions.

FROM --platform=linux/amd64 rocker/r-ver:4.4.1

ARG DEBIAN_FRONTEND=noninteractive
ARG MVSUSIER_REF=HEAD
ARG SUSIER_REF=HEAD
ARG SUSIER_VERSION=""

LABEL org.opencontainers.image.title="WARP AoU mvSuSiE" \
  org.opencontainers.image.description="R environment and AoU analysis scripts for CisMvSusie.wdl" \
  org.opencontainers.image.source="https://github.com/broadinstitute/warp-tools" \
  org.opencontainers.image.url="https://github.com/broadinstitute/warp-tools/tree/develop/3rd-party-tools/aou_mvsusie" \
  org.broadinstitute.warp-tools.scripts.source="https://github.com/AoU-Multiomics-Analysis/mvSusieR" \
  org.broadinstitute.warp-tools.mvsusier.source="https://github.com/stephenslab/mvsusieR" \
  org.broadinstitute.warp-tools.susier.source="https://github.com/stephenslab/susieR"

# Upgrade packages inherited from the base image as well as installing workflow tools.
RUN apt-get update -o APT::Update::Error-Mode=any \
    && apt-get upgrade -y --with-new-pkgs \
    && apt-get install -y --no-install-recommends \
      tabix gzip gawk coreutils procps ca-certificates curl git openssl \
      build-essential gfortran cmake \
      libgsl-dev libblas-dev liblapack-dev \
      libcurl4-openssl-dev libssl-dev libxml2-dev \
      libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
      libpng-dev libjpeg-dev libtiff-dev libcairo2-dev \
    && for package in openssl libssl3 libssl-dev; do \
      installed="$(dpkg-query -W -f='${Version}' "$package")" \
        && dpkg --compare-versions "$installed" ge 3.0.2-0ubuntu1.30 \
        || { echo "Expected $package >= 3.0.2-0ubuntu1.30 (Jammy); found ${installed:-not installed}." >&2; exit 1; }; \
    done \
    && rm -rf /var/lib/apt/lists/*

# R packages used by the scripts (and likely by the summary and plot scripts)
RUN R -q -e 'options(Ncpus = 4L); \
      install.packages(c("remotes", "optparse", "data.table", "R.utils", \
        "dplyr", "purrr", "readr", "tibble", "tidyr", "stringr", \
        "ggplot2", "patchwork", "cowplot", "jsonlite", "ragg", \
        "Matrix", "matrixStats", "ashr", "mashr", "susieR"))' \
 && R -q -e 'pkgs <- c("remotes", "optparse", "data.table", "R.utils", "dplyr", "purrr", \
        "readr", "tibble", "tidyr", "stringr", "ggplot2", "patchwork", "cowplot", \
        "jsonlite", "ragg", "Matrix", "matrixStats", "ashr", "mashr", "susieR"); \
      missing <- setdiff(pkgs, rownames(installed.packages())); \
      if (length(missing)) stop("Failed to install: ", paste(missing, collapse = ", "))'

RUN if [ -n "${SUSIER_VERSION}" ]; then \
      R -q -e "remotes::install_version('susieR', version = '${SUSIER_VERSION}', repos = 'https://cloud.r-project.org', upgrade = 'never'); stopifnot(as.character(packageVersion('susieR')) == '${SUSIER_VERSION}')"; \
    else \
      R -q -e 'remotes::install_github("stephenslab/susieR", ref = Sys.getenv("SUSIER_REF", "HEAD"), upgrade = "never")'; \
    fi \
 && R -q -e 'library(susieR); stopifnot(packageVersion("susieR") >= "0.15.54")'

RUN R -q -e 'remotes::install_github("stephenslab/mvsusieR", ref = Sys.getenv("MVSUSIER_REF", "HEAD"), upgrade = "never")' \
 && R -q -e 'library(mvsusieR); library(mashr); library(susieR); \
      cat("mvsusieR", as.character(packageVersion("mvsusieR")), \
          "susieR", as.character(packageVersion("susieR")), \
          "mashr", as.character(packageVersion("mashr")), "\n")'

# Record the exact package versions inside the image for traceability
RUN mkdir -p /opt/mvsusie && R -q -e 'writeLines(capture.output(sessionInfo()), "/opt/mvsusie/VERSIONS.txt")'

# Scripts: the WDL calls them as /opt/mvsusie/scripts/<name>.R
COPY scripts/ /opt/mvsusie/scripts/
RUN chmod +x /opt/mvsusie/scripts/*.R \
 && for f in prepare_trans_window.R prepare_window.R fit_window.R summarize_window.R \
             plot_window_mvsusie.R trans_window_cli.R trans_window_io.R trans_window_logging.R \
             trans_window_model.R trans_window_preprocess.R trans_window_prior.R \
             trans_window_prior_modes.R; do \
      test -s "/opt/mvsusie/scripts/$f" || { echo "Missing /opt/mvsusie/scripts/$f in the build context." >&2; exit 1; }; \
    done \
 && R -q -e 'for (f in list.files("/opt/mvsusie/scripts", "[.]R$", full.names = TRUE)) invisible(parse(f))'

CMD ["/bin/bash"]
