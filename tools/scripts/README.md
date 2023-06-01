# WARP Scripts

This directory contains the following scripts:

* `breakoutSnap.py` extracts the data in a snap file as csv files
* `create-merged-npz-output.py` takes a barcode.tsv, feature.tsv and matrix.mtx from STAR alignment outputs and creates 2 npy files and an npz file for row_index, col_index and the matrix. These files are required in the empty_drop step.
* `create_snss2_counts_csv.py` creates a csv file containing intron and exon counts from the Single Nucleus Smart-Seq2 pipeline
* `loomCompare.py` compares differences between loom files
* `ss2_loom_merge.py` creates a single loom file from multiple single sample loom files
* `makeCompliantBAM.py` make a BAM file with cellular barcodes in the read names compliant by moving them to the CB tag

The following scripts create a loom file from counts, metadata, and metrics from each pipeline:
* `create_loom_optimus.py` for Optimus pipeline
* `create_loom_snss2.py` for Single Nucleus Smart-Seq2 pipeline
* `create_snrna_optimus.py` for Optimus in `sn_rna` mode with `count_exons=false`
* `create_snrna_optimus_counts.py` for Optimus in `sn_rna` mode with `count_exons=true`

