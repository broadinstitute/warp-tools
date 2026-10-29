# One image for every task in CisMvSusie.wdl.
#
# Terra/Cromwell runs on Google Compute Engine VMs (Intel/AMD, linux/amd64) as root,
# starts the task with its own command (bash), and mounts the task directory at
# /cromwell_root. Nothing here needs an entrypoint or a working directory.
#
# Build on Google Cloud (native amd64, no emulation), after running ./get_scripts.sh:
#   gcloud builds submit --timeout=60m --machine-type=e2-highcpu-8 \
#     --tag <region>-docker.pkg.dev/<project>/<repo>/cis-mvsusie:0.1 .
# Or locally (slow on Apple Silicon):
#   docker buildx build --platform linux/amd64 -t <registry>/<name>:<tag> --push .
#
# Pinning: the R packages come from the dated Posit package snapshot that rocker/r-ver
# uses for this R version, so CRAN versions are reproducible. mvsusieR is installed from
# GitHub. Its API changed in its 2026 refactor (it now builds on susieR 2.0), so set
# MVSUSIER_REF (and SUSIER_VERSION if needed) to match what trans_window_model.R was
# written against. The "R session information" file written by FitMvSusie in the
# collaborator's runs (session_info.txt) lists the exact versions.

FROM --platform=linux/amd64 rocker/r-ver:4.4.1

ARG DEBIAN_FRONTEND=noninteractive
ARG MVSUSIER_REF=HEAD
ARG SUSIER_VERSION=""

# Command-line tools used by the WDL tasks: tabix and bgzip (htslib), gzip, awk, sort, comm
RUN apt-get update && apt-get install -y --no-install-recommends \
      tabix gzip gawk coreutils procps ca-certificates curl git \
      build-essential gfortran cmake \
      libgsl-dev libblas-dev liblapack-dev \
      libcurl4-openssl-dev libssl-dev libxml2-dev \
      libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
      libpng-dev libjpeg-dev libtiff-dev libcairo2-dev \
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
      R -q -e "remotes::install_version('susieR', version = '${SUSIER_VERSION}', upgrade = 'never')"; \
    fi

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
             trans_window_model.R trans_window_preprocess.R trans_window_prior.R; do \
      test -s "/opt/mvsusie/scripts/$f" || { echo "Missing /opt/mvsusie/scripts/$f. Run ./get_scripts.sh first." >&2; exit 1; }; \
    done \
 && R -q -e 'for (f in list.files("/opt/mvsusie/scripts", "[.]R$", full.names = TRUE)) invisible(parse(f))'

CMD ["/bin/bash"]
