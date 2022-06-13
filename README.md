
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ICAMS

<!-- badges: start -->

![R build
status](https://github.com/steverozen/ICAMS/workflows/R-CMD-check/badge.svg)
![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/steverozen/ICAMS?branch=master&svg=true)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ICAMS)](https://cran.r-project.org/package=ICAMS)
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)

<!-- badges: end -->

In-depth Characterization and Analysis of Mutational Signatures
(‘ICAMS’)

## Purpose

Analysis and visualization of experimentally elucidated mutational
signatures – the kind of analysis and visualization in Boot et al.,
“In-depth characterization of the cisplatin mutational signature in
human cell lines and in esophageal and liver tumors”, Genome Research
2018, <https://doi.org/10.1101/gr.230219.117> and “Characterization of
colibactin-associated mutational signature in an Asian oral squamous
cell carcinoma and in other mucosal tumor types”, Genome Research 2020
<https://doi.org/10.1101/gr.255620.119>. ‘ICAMS’ stands for In-depth
Characterization and Analysis of Mutational Signatures. ‘ICAMS’ has
functions to read in variant call files (VCFs) and to collate the
corresponding catalogs of mutational spectra and to analyze and plot
catalogs of mutational spectra and signatures. Handles both
“counts-based” and “density-based” (i.e. representation as mutuations
per megabase) mutational spectra or signatures.

## Installation

### To install the stable version from [CRAN](https://CRAN.R-project.org/package=ICAMS):

*IMPORTANT* Install the [Bioconductor](https://www.bioconductor.org/)
dependencies first:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("BSgenome")
```

This may be slow; please be patient.

Afterwards, install the latest release version of ICAMS with the R
command line:

``` r
install.packages("ICAMS")
```

### To install the latest stable version from [GitHub](https://github.com/steverozen/ICAMS):

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github(repo = "steverozen/ICAMS")
```

### Get the development version

To use new features in the development version, you can install ICAMS
from the master branch on
[GitHub](https://github.com/steverozen/ICAMS/tree/master), which may not
be stable:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github(repo = "steverozen/ICAMS", ref = "master")
```

## Reference manual

<https://github.com/steverozen/ICAMS/blob/master/data-raw/ICAMS_3.0.6.pdf>

## Frequently asked questions

#### How to do normalization for “counts-based” catalogs of mutational spectra or signatures to account for differing abundances of the source sequence of the mutations?

You can use *exported* function `TransformCatalog` in ICAMS to normalize
the data. Please refer to the documentation and example of
`TransformCatalog` for more details.

## Citing ICAMS

If you use ICAMS in your work, please cite:

> Rozen SG, Jiang NH, Boot A, Liu M, Wu Y, Huang MN, Chang JG (2021).
> ICAMS:In-depth Characterization and Analysis of Mutational Signatures.
> R package version 3.0.6, <https://CRAN.R-project.org/package=ICAMS>.
