#!/bin/bash
## Unit test for emptyDropsWrapper.R: the same --seed must give identical output,
## and running without --seed must still work.
set -euo pipefail
cd "$(dirname "$0")"
wrapper=../emptyDropsWrapper.R
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

## Synthetic gene x droplet matrix: 200 cell-like droplets plus ambient ones
Rscript -e 'suppressMessages(library(Matrix)); set.seed(1); g <- 500
  tot <- c(rpois(200, 2000), rpois(2800, 30))
  m <- as(sapply(tot, function(t) tabulate(sample(g, t, TRUE), g)), "CsparseMatrix")
  dimnames(m) <- list(paste0("g", 1:g), paste0("c", seq_along(tot)))
  saveRDS(m, commandArgs(TRUE)[1])' "$tmp/m.rds"

run() { $wrapper -i "$tmp/m.rds" -o "$tmp/$1.csv" --emptydrops-niters 1000 --min-molecules 100 "${@:2}"; }
run seeded1 --seed 42
run seeded2 --seed 42
run unseeded

if cmp -s "$tmp/seeded1.csv" "$tmp/seeded2.csv"; then
  echo "PASSED: --seed 42 output is identical across runs"
else
  echo "FAILED: --seed 42 output differs across runs"; exit 1
fi
