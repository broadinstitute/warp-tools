"""Check every query input before any sample is processed, and report every problem at once.

  python3 validate_inputs.py --input-ids a b --gex a.h5ad b.h5ad [--atac a_atac.h5ad b_atac.h5ad] \
      --outputs a_annotated.h5ad b_annotated.h5ad ...
"""
import argparse
import os
import sys

import h5py


def file_problems(path, kind):
    if not os.path.isfile(path):
        return [f"{kind} file not found: {path}"]
    if os.path.getsize(path) == 0:
        return [f"{kind} file is empty: {path}"]
    try:
        with h5py.File(path, "r") as f:
            missing = [group for group in ("X", "obs", "var") if group not in f]
    except OSError as e:
        return [f"{kind} file is not an h5ad (HDF5) file: {path} ({e})"]
    return [f"{kind} file {path} has no '{group}'" for group in missing]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-ids", nargs="+", required=True)
    p.add_argument("--gex", nargs="+", required=True)
    p.add_argument("--atac", nargs="*", default=[])
    p.add_argument("--outputs", nargs="*", default=[], help="Every output file name of the run.")
    a = p.parse_args()

    errors = []
    if len(a.gex) != len(a.input_ids) or (a.atac and len(a.atac) != len(a.gex)):
        errors.append(f"{len(a.input_ids)} input ids, {len(a.gex)} GEX files and {len(a.atac)} ATAC "
                      "files do not line up (ATAC must be empty or one per GEX file)")
    duplicates = sorted({i for i in a.input_ids if a.input_ids.count(i) > 1})
    if duplicates:
        errors.append(f"input ids must be unique: {duplicates}")
    if not all(a.input_ids):
        errors.append("input ids cannot be empty")
    # Outputs are named after the input files, so two inputs with the same file name collide.
    clashes = sorted({o for o in a.outputs if a.outputs.count(o) > 1})
    if clashes:
        errors.append(f"input files must have distinct names, because outputs are named after them: "
                      f"{clashes} would be written twice")
    for i, input_id in enumerate(a.input_ids):
        files = [(a.gex[i] if i < len(a.gex) else None, "GEX"),
                 (a.atac[i] if i < len(a.atac) else None, "ATAC")]
        for path, kind in files:
            if path:
                errors += [f"{input_id}: {problem}" for problem in file_problems(path, kind)]

    if errors:
        print("ERROR: input check failed; no sample was processed:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        sys.exit(1)
    print(f"Input check passed for {len(a.input_ids)} sample(s).", flush=True)


if __name__ == "__main__":
    main()
