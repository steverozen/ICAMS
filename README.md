
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ICAMS

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/steverozen/ICAMS.svg?branch=master)](https://travis-ci.org/steverozen/ICAMS)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/ICAMS)](https://cran.r-project.org/package=ICAMS)

<!-- badges: end -->

In-depth Characterization and Analysis of Mutational Signatures
(‘ICAMS’)

## Purpose

Analysis and visualization of experimentally elucidated mutational
signatures – the kind of analysis and visualization in Boot et al.,
“In-depth characterization of the cisplatin mutational signature in
human cell lines and in esophageal and liver tumors”, Genome Research
2018, <https://doi.org/10.1101/gr.230219.117>. ‘ICAMS’ stands for
In-depth Characterization and Analysis of Mutational Signatures. ‘ICAMS’
has functions to read in variant call files (VCFs) and to collate the
corresponding catalogs of mutational spectra and to analyze and plot
catalogs of mutational spectra and signatures. Handles both
“counts-based” and “density-based” catalogs of mutational spectra or
signatures.

## Installation

Install the stable version of ICAMS from
[CRAN](https://cran.r-project.org/) with the R command line:

``` r
install.packages("ICAMS")
```

After that, install the necessary dependency package from
[Bioconductor](https://www.bioconductor.org/) in order to successfully
load ICAMS:

``` r
install.packages("BiocManager")
BiocManager::install("BSgenome")
```

### Development version

To use new features, you can install ICAMS from the master branch on
[GitHub](https://github.com/), which may not be stable:

``` r
install.packages("devtools")
devtools::install_github("steverozen/ICAMS")
```

## Reference manual

<https://github.com/steverozen/ICAMS/blob/master/data-raw/ICAMS_2.0.9.9010.pdf>
