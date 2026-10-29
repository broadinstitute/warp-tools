#!/usr/bin/env bash
# Copies the complete scripts/ folder from the collaborator's repository into ./scripts.
# The repository copy replaces the five files in this zip, which are the ones shared in chat;
# the repository also has the six this zip does not include (trans_window_cli.R,
# trans_window_logging.R, trans_window_model.R, trans_window_prior.R, summarize_window.R,
# plot_window_mvsusie.R).
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
git clone --depth 1 https://github.com/AoU-Multiomics-Analysis/mvSusieR.git "$tmp/repo"
mkdir -p "$here/scripts"
cp -R "$tmp/repo/scripts/." "$here/scripts/"
echo "Copied: $(ls "$here/scripts" | tr '\n' ' ')"
echo "Repository commit: $(git -C "$tmp/repo" rev-parse HEAD)"
