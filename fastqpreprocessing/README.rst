Single Cell Tools
#################

.. image:: https://img.shields.io/circleci/project/github/HumanCellAtlas/sctools.svg?label=Unit%20Test%20on%20Circle%20CI%20&style=flat-square&logo=circleci
  :target: https://circleci.com/gh/HumanCellAtlas/sctools/tree/master
  :alt: Unit Test Status

.. image:: https://img.shields.io/codecov/c/github/HumanCellAtlas/sctools/master.svg?label=Test%20Coverage&logo=codecov&style=flat-square
  :target: https://codecov.io/gh/HumanCellAtlas/sctools
  :alt: Test Coverage on Codecov

.. image:: https://img.shields.io/readthedocs/sctools/latest.svg?label=ReadtheDocs%3A%20Latest&logo=Read%20the%20Docs&style=flat-square
  :target: http://sctools.readthedocs.io/en/latest/?badge=latest
  :alt: Documentation Status

.. image:: https://img.shields.io/snyk/vulnerabilities/github/HumanCellAtlas/sctools/requirements.txt.svg?label=Snyk%20Vulnerabilities&logo=Snyk
  :target: https://snyk.io/test/github/HumanCellAtlas/sctools/?targetFile=requirements.txt
  :alt: Snyk Vulnerabilities for GitHub Repo (Specific Manifest)

.. image:: https://img.shields.io/github/release/HumanCellAtlas/sctools.svg?label=Latest%20Release&style=flat-square&colorB=green
  :target: https://github.com/HumanCellAtlas/sctools/releases
  :alt: Latest Release

.. image:: https://img.shields.io/github/license/HumanCellAtlas/sctools.svg?style=flat-square
  :target: https://img.shields.io/github/license/HumanCellAtlas/sctools.svg?style=flat-square
  :alt: License

.. image:: https://img.shields.io/badge/python-3.6-green.svg?style=flat-square&logo=python&colorB=blue
  :target: https://img.shields.io/badge/python-3.6-green.svg?style=flat-square&logo=python&colorB=blue
  :alt: Language

.. image:: https://img.shields.io/badge/Code%20Style-black-000000.svg?style=flat-square
  :target: https://github.com/ambv/black
  :alt: Code Style

Single Cell Tools provides utilities for manipulating sequence data formats suitable for use in
distributed systems analyzing large biological datasets.

Download and Installation
=========================

.. code bash
   git clone https://github.com/humancellatlas/sctools.git
   cd sctools
   pip3 install .
   pytest  # verify installation; run tests

sctools Package
===============

The sctools package provides both command line utilities and classes designed for use in python
programs.

Command Line Utilities
======================

1. Attach10XBarcodes: Attached barcodes stored in fastq files to reads in an unaligned bam file
2. SplitBam: Split a bam file into chunks, guaranteeing that cells are contained in 1 chunk
3. CalculateGeneMetrics: Calculate information about genes in an experiment or chunk
4. CalculateCellMetrics: Calculate information about cells in an experiment or chunk
5. MergeGeneMetrics: Merge gene metrics calculated from different chunks of an experiment
6. MergeCellMetrics Merge cell metrics calculated from different chunks of an experiment

Main Package Classes
====================

1. **Platform**: an abstract class that defines a common data structure for different 3' sequencing
   formats. All algorithms and methods in this package that are designed to work on 3' sequencing data
   speak to this common data structure. Currently 10X_v2 is defined.

2. **Reader**: a general iterator over arbitrarily zipped file(s) that is extended to work with common
   sequence formats like fastq (fastq.Reader) and gtf (gtf.Reader). We recommend using the pysam
   package for reading sam and bam files.

3. **TwoBit & ThreeBit** DNA encoders that store DNA in 2- and 3-bit form. 2-bit is smaller but
   randomizes "N" nucleotides. Both classes support fastq operations over common sequence tasks such
   as the calculation of GC content.

4. **ObservedBarcodeSet & PriorBarcodeSet**: classes for analysis and comparison of sets of barcodes
   such as the cell barcodes used by 10X genomics. Supports operations like summarizing hamming
   distances and comparing observed sequence diversity to expected (normally uniform) diversity.

5. **gtf.Reader & gtf.Record** GTF iterator and GTF record class that exposes the gtf
   fields as a lightweight, lazy-parsed python object.

6. **fastq.Reader & fastq.Record** fastq reader and fastq record class that exposes the fastq fields
   as a lightweight, lazy-parsed python object.

7. **Metrics** calculate information about the genes and cells of an experiment

8. **Bam** Split bam files into chunks and attach barcodes as tags


Viewing Test Results and Coverage
=================================

To calculate and view test coverage cd to the ``sctools`` directory and
type the following two commands to generate the report and open it in your web browser:

.. code:: bash

   pytest --cov-report html:cov_html --cov=sctools
   open cov_html/index.html

Definitions
===========

Several definitions are helpful to understand how sequence data is analyzed.

1. **Cell**: an individual cell, the target of single-cell RNA-seq experiments and the entity that we
wish to characterize

2. **Capture Primer**: A DNA oligonucleotide containing amplification machinery, a fixed cell barcode,
a random molecule barcode, and an oligo-dT tail to capture poly-adenylated RNA

3. **Molecule**: A molecule refers to a single mRNA molecule that is captured by an oligo-dT capture
primer in a single-cell sequencing experiment

4. **Molecule Barcode**: A molecule barcode (alias: UMI, RMT) is a short, random DNA barcode attached
to the capture primer that has adequate length to be probabilistically unique across the experiment.
Therefore, when multiple molecules of the same gene are captured in the same cell, they can be
differentiated through having different molecule barcodes. The proposed GA4GH standard tag for a
molecule barcode is UB and molecule barcode qualities is UY

5. **Cell Barcode**: A short DNA barcode that is typically selected from a whitelist of barcodes that
will be used in an experiment. All capture primers for a given cell will contain the same cell
barcode. The proposed GA4GH standard tag for a cell barcode is CB and cell barcode qualities is CY

6. **Fragment**: During library construction, mRNA molecules captured on capture primers are amplified,
and the resulting amplified oligonucleotides are fragmented. In 3' experiments, only the fragment
that contains the 3' end is retained, but the break point will be random, which means fragments
often have different lengths. Once sequenced, different fragments can be identified as unique
combinations of cell barcode, molecule barcode, the chromosome the sequence aligns to, and the
position it aligns to on that chromosome, after correcting for clipping that the aligner may add

7. **Bam/Sam file**: The GA4GH standard file type for the storage of aligned sequencing reads.
Unless specified, our Single Cell Tools will operate over bam files containing either aligned or
unaligned reads

Development
===========

Code Style
----------
The sctools code base is complying with the PEP-8 and using `Black <https://github.com/ambv/black>`_ to
format our code, in order to avoid "nitpicky" comments during the code review process so we spend more time discussing about the logic, 
not code styles.

In order to enable the auto-formatting in the development process, you have to spend a few seconds setting 
up the ``pre-commit`` the first time you clone the repo:

1. Install ``pre-commit`` by running: ``pip install pre-commit`` (or simply run ``pip install -r requirements.txt``).
2. Run `pre-commit install` to install the git hook.

Once you successfully install the ``pre-commit`` hook to this repo, the Black linter/formatter will be automatically triggered and run on this repo. Please make sure you followed the above steps, otherwise your commits might fail at the linting test!

If you really want to manually trigger the linters and formatters on your code, make sure ``Black`` and ``flake8`` are installed in your Python environment and run ``flake8 DIR1 DIR2`` and ``black DIR1 DIR2 --skip-string-normalization`` respectively.
