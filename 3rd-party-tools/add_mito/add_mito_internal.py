#!/usr/bin/env python3
"""
Author: Matthew Schmitz, Allen Institute, 2024
Modified to work inside Docker container with MitoFinder installed

Append mitochondrial annotations (via MitoFinder) to a genome assembly.
This version runs MitoFinder directly inside the container.

usage:
    add_mito -a NC_012920.1 \
        -r /data/mito_ref.gbk \
        -g /data/genome.fa \
        -t /data/annotations.gtf \
        -n rhemac10 \
        -e your.email@example.com
"""

import sys
import os
import argparse
import requests
import subprocess
import re

def download_fasta_http(accession: str,
                        out_path: str = None,
                        email: str = "your.email@example.com",
                        api_key: str = None) -> str:
    """
    Fetches the FASTA via HTTP GET to NCBI E-utilities.
    """

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nucleotide",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        "email": email,
    }
    if api_key:
        params["api_key"] = api_key
    r = requests.get(base, params=params)
    r.raise_for_status()
    if out_path is None:
        out_path = accession + ".fasta"
    with open(out_path, "w") as f:
        f.write(r.text)
    return out_path

def gff_to_gtf_with_transcripts(gff_path: str, gtf_path: str):
    """
    Convert a MitoFinder GFF → GTF, adding:
      • gene features with gene_id/gene/gene_name = Name
      • transcript features (type 'transcript'), same span as each gene,
        with transcript_id/gene_id/gene_name/transcript_name = Name(.t1)
    All other lines are converted 1:1 (just with gene_id/gene/gene_name attrs).
    """
    with open(gff_path) as gf, open(gtf_path, "w") as out:
        for line in gf:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            feature_type = cols[2]
            raw_attrs    = cols[8]

            # pull out the Name= tag
            name = None
            for tok in re.split(r"[;\s]+", raw_attrs):
                if tok.startswith("Name="):
                    _, name = tok.split("=", 1)
                    break
            if not name:
                name = "."

            # build the common gene‐level attribute string
            gene_attrs = f'gene_id "{name}"; gene "{name}"; gene_name "{name}";'

            # write the gene line for any 'gene' feature
            if feature_type == "gene":
                # 1) gene
                cols[8] = gene_attrs
                out.write("\t".join(cols) + "\n")

                # 2) transcript
                tx_cols = cols.copy()
                tx_cols[2] = "transcript"
                tx_cols[8] = (
                    f'gene_id "{name}"; '
                    f'transcript_id "{name}.t1"; '
                    f'gene_name "{name}"; '
                    f'transcript_name "{name}.t1";'
                )
                out.write("\t".join(tx_cols) + "\n")

            else:
                # everything else stays a single line, tagged as gene
                cols[8] = gene_attrs
                out.write("\t".join(cols) + "\n")

def run_mitofinder_and_append(
        mito_accession: str,
        ref_gbk: str,
        genomic_fasta: str,
        transcript_gtf: str,
        spec_name: str,
        mf_opts: list,
        email: str = None,
        api_key: str = None,
) -> (str, str):
    """
    1) Download mitochondrial FASTA by accession.
    2) Run MitoFinder (installed in container).
    3) Convert its GFF → GTF, append to transcript_gtf → *_mito.gtf.
    4) Append mito.fa to genomic_fasta → *_mito.fasta.
    Returns (new_genome_fasta, new_transcript_gtf).
    """
    # --- 1) download mito fasta
    mito_fasta = download_fasta_http(
        mito_accession,
        out_path=f"{mito_accession}.fasta",
        email=email,
        api_key=api_key
    )

    # --- 2) run MitoFinder directly (assuming it's installed in the container)
    cmd = [
              "mitofinder",
              "-j", spec_name,
              "-a", mito_fasta,
              "-r", ref_gbk
          ] + mf_opts

    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # --- 3) locate and convert its GFF
    results_dir = os.path.join(
        spec_name, f"{spec_name}_MitoFinder_mitfi_Final_Results"
    )
    gff_file = os.path.join(results_dir, f"{spec_name}_mtDNA_contig.gff")
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"Expected {gff_file}")

    mito_gtf = os.path.join(results_dir, f"{spec_name}_mtDNA_contig.gtf")
    gff_to_gtf_with_transcripts(gff_file, mito_gtf)

    # --- 4) write new transcriptome GTF
    base_gtf = os.path.splitext(transcript_gtf)[0]
    out_gtf = f"{base_gtf}_mito.gtf"
    with open(out_gtf, "w") as w:
        # original annotations
        with open(transcript_gtf) as rf:
            w.write(rf.read().rstrip("\n") + "\n")
        # mito annotations
        with open(mito_gtf) as mf:
            w.write(mf.read())

    # --- 5) write new genomic FASTA
    base_fa = os.path.splitext(genomic_fasta)[0]
    out_fa = f"{base_fa}_mito.fasta"
    with open(out_fa, "w") as w:
        # original genome
        with open(genomic_fasta) as gf:
            w.write(gf.read().rstrip("\n") + "\n")
        # append mito-chromosome
        with open(mito_fasta) as mf:
            w.write(mf.read())

    return out_fa, out_gtf

def main():
    parser = argparse.ArgumentParser(
        description="Download mito‐chromosome, run MitoFinder, and append "
                    "mitochondrial annotations to a genome FASTA and GTF."
    )
    parser.add_argument(
        '-a', '--accession',
        required=True,
        help="NCBI GenBank accession for the mitochondrial chromosome (e.g. NC_012920.1)"
    )
    parser.add_argument(
        '-r', '--ref_gbk',
        required=True,
        help="Path to the reference GenBank .gbk file for the mitochondrion"
    )
    parser.add_argument(
        '-g', '--genome_fasta',
        required=True,
        help="Path to the genomic FASTA to which the mito sequence will be appended"
    )
    parser.add_argument(
        '-t', '--transcript_gtf',
        required=True,
        help="Path to the transcriptome GTF to which the mito GTF will be appended"
    )
    parser.add_argument(
        '-n', '--spec_name',
        required=True,
        help="Sample/species name (used by MitoFinder to name its output folder)"
    )
    parser.add_argument(
        '-e', '--email',
        help="NCBI Entrez email address (for EFetch)"
    )
    parser.add_argument(
        '-k', '--api_key',
        help="NCBI Entrez API key (optional, raises rate limits)"
    )
    parser.add_argument(
        '-m','--mitofinder_opts',
        nargs=argparse.REMAINDER,
        default=[
            "-o", "2",
            "--blast-identity-prot", "30",
            "--blast-size", "20",
            "--blast-identity-nucl", "30"
        ],
        help=(
            "Extra args passed to MitoFinder (default: "
            "-o 2 --blast-identity-prot 30 --blast-size 20 --blast-identity-nucl 30). "
            "Specify your own by appending them after this flag."
        )
    )
    args = parser.parse_args()

    print("Running MitoFinder and appending mito sequence/annotations…")
    try:
        new_fasta, new_gtf = run_mitofinder_and_append(
            mito_accession=args.accession,
            ref_gbk=args.ref_gbk,
            genomic_fasta=args.genome_fasta,
            transcript_gtf=args.transcript_gtf,
            spec_name=args.spec_name,
            email=args.email,
            api_key=args.api_key,
            mf_opts=args.mitofinder_opts
        )
    except Exception as e:
        print(f"✖ Error: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Combined FASTA written to: {new_fasta}")
    print(f"Combined GTF written to:   {new_gtf}")

if __name__ == "__main__":
    main()