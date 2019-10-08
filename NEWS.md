# ICAMS 2.0.9.000x
* Added public functions AnnotateSBSVCF, AnnotateDBSVCF and
  AnnotateIDVCF.
* Changed column name 'gene.name' to 'gene.symbol' in trans.ranges.
* In FindDelMH, cryptic repeats (i.e. un-normalized deletions in a repeat 
  such as GAGG deleted from CCCAGGGAGGGTCCC, which should be normalized
  to a deletion of AGGG) are now ignored with a warning rather than
  causing a stop().
* Fixed a rarely encountered bug in FindDelMH, which previously did not flag the
  cryptic repeat in what is now the second example.
* Updated plotting functions to make y axis labels more informative.
* Made the return from VCFsToIDCatalogs a list; 1st element is the spectrum catalog 
  (previously the only return); 2nd element is a list of further annotated VCFs.

# ICAMS 2.0.9
* as.catalog supports creation of the catalog from a vector (interpreted
  as a 1-column matrix) and optionally infers the class from the
  number of rows in the input.
* Exported some functions that are used in categorizing deletions as 
  a way to make documentation of the underlying logic more visible.

# ICAMS 2.0.8
* Moved two Bioconductor packages(BSgenome.Hsapiens.1000genomes.hs37d5, 
  BSgenome.Hsapiens.UCSC.hg38) from Imports to Suggests so that they can
  be used conditionally in "ICAMS" as requested by CRAN.
* Added ability to calculate VAF from newer version Mutect VCF files.
* Documentation updates.

# ICAMS 2.0.7
* Additional corrections for submission to CRAN
  (putting on.exit(par(opar)) immediately after
  opar <- par(...)
  in every case.

# ICAMS 2.0.6
* Additional changes for submission to CRAN (re-set
  graphics parameters to original values after plotting
  calls; documentation updates).

# ICAMS 2.0.5
* Changes in response to initial submission to CRAN.
* Added support for BSgenome.Mmusculus.UCSC.mm10. 

# ICAMS 2.0
* First stable release. Mutational catalogs are S3 classes.

# ICAMS 1.0
* Only available on GitHub; interfaces not compatible with ICAMS 2.0.

