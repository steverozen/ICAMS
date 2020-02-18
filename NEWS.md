# ICAMS 2.1.0.900x
* Slight changes to behavior of PlotCatalog() and PlotCatalogToPdf() for class
  SBS96Catalog: (New) Allow setting ylim and cex.
  (New) For PlotCatalog (not PlotCatalogToPdf), allow plotting of a 96 x 2 catalog,
  in which case behavior is a stacked
  bar chart. (New) Plot x axis tick marks if xlabels is not TRUE; set
  \code{par(tck = 0)} to suppress. 
  ? Minor changes in default font sizes.
* Added a new dependency package "zip" in ICAMS, to be used in three new exported  
  functions MutectVCFFilesToZipFile(), StrelkaSBSVCFFilesToZipFile() and 
  StrelkaIDVCFFilesToZipFile(). 
* Renamed internal function creating custom abundance to GetCustomKmerCounts().
* Fixed a bug in function TransformCatalog() and updated its documentation
  for parameter "target.abundance".
* Deprecated counts -> counts transformation in TransformCatalog(); see documentation
  for rationale.
* Added code to handle multiple alternate alleles in Strelka VCFs, e.g cases 
  where REF = T and ALT = A,C which indicates the alternate allele in some reads
  is A and in some reads C.
* Added an internal function to infer trans ranges from input ref.genome.
  Updated the documentation of the exported functions which have argument "trans.ranges"
  to make it optional.
* Updated documentation for optional arguments in the functions reading VCF files.
  Used @inheritParams to reduce repetitive documentation.
* Updated function PlotTransBiasGeneExpToPdf() so that ymax on the plot will be changed 
based on plot.type
* Fixed bugs in internal function CheckAndFixChrNames() and updated the automated tests.
* Corrected error message in TransformCatalog.
* Fixed a bug in exported function GetMutectVAF() and updated the warning
  message to make it more informative. 
* Fixed bugs in functions generating zip archive from VCF files and added 
  more automated test. 
* Added an additional argument "name.of.VCF" in internal functions 
  ReadStrelkaSBSVCF(), ReadStrelkaIDVCF() and exported function GetStrelkaVAF().
* Added an optional argument "flag.mismatches" in functions VCFsToIDCatalogs(),
  MutectVCFFilesToCatalog(), MutectVCFFilesToCatalogAndPlotToPdf(),
  MutectVCFFilesToZipFile(), StrelkaIDVCFFilesToCatalog(),
  StrelkaIDVCFFilesToCatalogAndPlotToPdf() and StrelkaIDVCFFilesToZipFile().

# ICAMS 2.0.10
* (Non backward-compatible.) 
  Made the return from functions VCFsToIDCatalogs(), StrelkaIDVCFFilesToCatalog()
  and StrelkaIDVCFFilesToCatalogAndPlotToPdf() a list; 1st element is the
  spectrum catalog (previously the only return); 2nd element is a list of
  VCFs with additional annotations.
* (Non backward-compatible.) 
  Made the return from function PlotCatalog() a list. The first element is 
  a logical value indicating whether the plot is successful. The second element 
  is a numeric vector giving the coordinates of all the bar midpoints drawn,
  useful for adding to the graph (currently only implemented for SBS96Catalog).
* (Slightly non-backward compatible.) Changed the handling of the output.file argument in
  MutectVCFFilesToCatalogAndPlotToPdf(), StrelkaSBSVCFFilesToCatalogAndPlotToPdf(), 
  and StrelkaIDVCFFilesToCatalogAndPlotToPdf()
  so that an indicator of the catalog type plus ".pdf" is simply
  appended to the base output.file name. Also made this argument
  optional with sensible default behavior.
* Updated methods so that they are compatible with new versions of
  R, in which the class of a matrix is c("matrix", "array").
* Added public functions AnnotateSBSVCF, AnnotateDBSVCF and
  AnnotateIDVCF.
* Changed column name 'gene.name' to 'gene.symbol' and added Ensembl gene ID 
  information in package variable trans.ranges. 
* In FindDelMH, cryptic repeats (i.e. un-normalized deletions in a repeat 
  such as GAGG deleted from CCCAGGGAGGGTCCC, which should be normalized
  to a deletion of AGGG) are now ignored with a warning rather than
  causing a stop().
* Fixed a rarely encountered bug in FindDelMH, which previously did not flag the
  cryptic repeat in what is now the second example in the function documentation.
* Updated plotting functions to make y axis labels more informative.
* Updated splitting functions for Mutect VCF to detect and move complex
  indels to "other" category.
* Added two exported functions PlotTransBiasGeneExp() and  PlotTransBiasGeneExpToPdf().
* Added an additional optional argument "names.of.VCFs" in functions
  ReadAndSplitMutectVCFs(), ReadAndSplitStrelkaSBSVCFs(), ReadStrelkaIDVCFs(),
  MutectVCFFilesToCatalog(), MutectVCFFilesToCatalogAndPlotToPdf(),
  StrelkaIDVCFFilesToCatalog(), StrelkaIDVCFFilesToCatalogAndPlotToPdf(),
  StrelkaSBSVCFFilesToCatalog and StrelkaSBSVCFFilesToCatalogAndPlotToPdf()
  for users to specify the names of samples in the VCF files.
* Added error checking on the order of rownames in the input of as.catalog().
* Added code to handle the case where the #CHROM lines appear
  interspersed with rows for variants in a Mutect VCF file.
* Added two exported package variables, gene.expression.data.HepG2 and
  gene.expression.data.MCF10A.
* Added an additional optional argument "tumor.col.names" in functions
  ReadAndSplitMutectVCFs(), MutectVCFFilesToCatalog() and
  MutectVCFFilesToCatalogAndPlotToPdf() to specify the column of the VCF
  that contains sequencing statistics such as sequencing depth; this column
  is often called "unknown" in Mutect.
* Added a "Comments" section in the documentation of functions     
  MutectVCFFilesToCatalog(),
  MutectVCFFilesToCatalogAndPlotToPdf(), StrelkaSBSVCFFilesToCatalog(),
  StrelkaSBSVCFFilesToCatalogAndPlotToPdf(), VCFsToSBSCatalogs(),
  VCFsToDBSCatalogs(), ReadCatalog() informing the user how to change
  attributes of the generated catalog.

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

