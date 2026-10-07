# emptyDropsWrapper

Reads in an rds file containing the count matrix in dgCMatrix format in gene x droplet orientation and produces a CSV file indicating cell filtering and other cell metadat as returned by the emptyDrops function.

## Requirements

This script requires the following R packages to be installed: 

* optparse
* DropletUtils

DropletUtils are available from bioconductor and require R version 3.5

emptyDrops p-values come from a Monte-Carlo simulation, so output varies between runs unless `--seed` is given.

https://bioconductor.org/packages/release/bioc/html/DropletUtils.html

## Testing

The test checks that two runs with the same `--seed` give identical output (CI runs it on every image build):

```
cd /tools/emptyDropsWrapper/test/ 
./test_emptyDropsWrapper.sh
```
