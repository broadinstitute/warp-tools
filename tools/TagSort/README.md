## Building and running manually (outside of Docker)

In production, this code is expected to be built and deployed by the
`docker_build.sh` script one directory up from here. To compile and run locally,
run `./fetch_and_make_dep_libs.sh && make`.

Here are some working examples of invoking TagSort, both for gathering metrics
and for actually sorting the input. They also show how the barcode/UMI/gene tag
order is determined by the order of the command line flags.

First, here are examples of the three different types of metrics you can ask for
(cell, gene, UMI). In these examples, TagSort is not actually writing the sorted
output. (The sorting is still needed to efficiently gather those metrics, which
is why the metric gathering isn't a separate program).

* Cell metrics. Tag order of input BAM is expected to be (barcode, umi, gene)
  because `--barcode-tag` `--umi-tag` `--gene-tag` are specified in that order.
  ```
  ./bin/TagSort --bam-input example_input.bam \
      --gtf-file example_species_genome_annotation.gtf \
      --metric-output example_output.cell-metrics.csv \
      --compute-metric \
      --metric-type cell \
      --barcode-tag CB \
      --umi-tag UB \
      --gene-tag GX \
      --temp-folder /tmp/ \
      --alignments-per-thread 1000000 \
      --nthreads 2 \
      --mitochondrial-gene-names-filename example_mito_genes.txt
  ```

* Gene metrics. Tag order here is (gene, barcode, umi) because the flags come in
  the order `--gene-tag`, `--barcode-tag`, `--umi-tag`.
  ```
  ./bin/TagSort --bam-input example_input.bam \
      --metric-output example_output.gene-metrics.csv \
      --compute-metric \
      --metric-type gene \
      --gene-tag GX \
      --barcode-tag CB \
      --umi-tag UB \
      --temp-folder /tmp/ \
      --alignments-per-thread 1000000 \
      --nthreads 2 \
      --mitochondrial-gene-names-filename example_mito_genes.txt
  ```

* UMI metrics. Tag order is again (gene, barcode, umi).
  ```
  ./bin/TagSort --bam-input example_input.bam \
      --metric-output example_output.umi-metrics.csv \
      --compute-metric \
      --metric-type umi \
      --gene-tag GX \
      --barcode-tag CB \
      --umi-tag UB \
      --temp-folder /tmp/ \
      --alignments-per-thread 1000000 \
      --nthreads 2 \
      --mitochondrial-gene-names-filename example_mito_genes.txt
  ```

Next, here's how to use TagSort to actually receive a sorted version of the BAM
input file. Currently the input must be BAM, while the sorted output will be
tab-separated ASCII SAM.
```
  ./bin/TagSort --bam-input example_input.bam \
      --output-sorted-info
      --sorted-output example_sorted_output.sam
      --gene-tag GX \
      --barcode-tag CB \
      --umi-tag UB \
      --temp-folder /tmp/ \
      --alignments-per-thread 1000000 \
      --nthreads 2 \
      --mitochondrial-gene-names-filename example_mito_genes.txt
  ```

Final notes:
* `--nthreads` is set to 2 because these examples are meant for
  playing around on a laptop or small VM; in production it should be as many
  cores as the VM has.
* TagSort can only compute one metric type per run.
* TagSort can also compute a metric in the same run that it writes the sorted
  output; simply specify both `--compute-metric` and `--output-sorted-info`.
* `--mitochondrial-gene-names-filename` is optional. If you don't provide it,
  TagSort will only consider genes that that with mt- or MT-.

## Unit tests

To run the unit tests, run `./fetch_gtest.sh && make test`, then run each of the
resulting test executables that get produced in the `bin/` directory.

To add a new file full of tests "foo_test.cpp", add bin/foo_test to the
Makefile's `test:` target, then add foo_test.cpp to the `test/` directory.
Note that the name must end in _test for `make` to know how to handle it.
