## Building and running manually (outside of Docker)

The fastqpreprocessing suite comprises four programs:
* `fastqprocess`, repartitions a batch of FASTQ files to keep all reads from a
  given cell in the same file. Input must have the read structure of the 3'
  transciptomic assay produced by 10x Genomics. Output can be chosen to remain
  FASTQ, or be converted to BAM.
* `fastq_slideseq`, same as fastqprocess, but with the added flexibility of a
  user-specified read structure. Added to handle the slideseq assay.
* `fastq_metrics`, summarizes total counts of UMIs and cell barcodes, and
  position weight matrices of UMIs and cell barcodes, from all reads.
* `samplefastq`, a filter, keeping just the reads matching a user-specified cell
  barcode whitelist. Requires read structure of 8C18X6C9M1X with a fixed spacer
  sequence.

In production, this code is expected to be built and deployed by the
`docker_build.sh` script one directory up from here. To compile and run locally,
run `./fetch_and_make_dep_libs.sh && make`.

Both `fastqprocess` and `fastq_slideseq` can take multiple copies of the R1, R2,
and I1 input files, specified by multiple instances of the `--R1`, `--R2`,
`--I1` flags (as below). There must be as many R1s as R2s. The R1s and R2s must
be specified in the same order. I1s are optional, but if present, there must be
as many as there are R1s.

Examples:

```
./bin/fastqprocess --verbose \
  --bam-size 0.001 \
  --barcode-length 16 \
  --umi-length 10 \
  --sample-id L8TX \
  --white-list data/L8TX/737K-august-2016.txt \
  --I1 data/L8TX/A_I1.fastq.gz \
  --R1 data/L8TX/A_R1.fastq.gz \
  --R2 data/L8TX/A_R2.fastq.gz \
  --I1 data/L8TX/B_I1.fastq.gz \
  --R1 data/L8TX/B_R1.fastq.gz \
  --R2 data/L8TX/B_R2.fastq.gz
```

```
./bin/fastq_slideseq  \
  --bam-size 30.0 \
  --white-list WhiteList.txt \
  --read-structure 11C22M \
  --sample-id EXAMPLEID \
  --output-format FASTQ \
  --I1 data/EXAMPLEID/A_I1.fastq.gz \
  --R1 data/EXAMPLEID/A_R1.fastq.gz \
  --R2 data/EXAMPLEID/A_R2.fastq.gz \
  --I1 data/EXAMPLEID/B_I1.fastq.gz \
  --R1 data/EXAMPLEID/B_R1.fastq.gz \
  --R2 data/EXAMPLEID/B_R2.fastq.gz
```

```
./bin/fastq_metrics \
  --white-list WhiteList.txt \
  --read-structure 11C22M \
  --sample-id EXAMPLEID \
  --R1 data/EXAMPLEID/A_R1.fastq.gz \
  --R1 data/EXAMPLEID/B_R1.fastq.gz
```

## Unit tests

To run the unit tests, run `./fetch_gtest.sh && make test`, then run each of the
resulting test executables that get produced in the `bin/` directory.

To add a new file full of tests "foo_test.cpp", add bin/foo_test to the
Makefile's `test:` target, then add foo_test.cpp to the `test/` directory.
Note that the name must end in _test for `make` to know how to handle it.
