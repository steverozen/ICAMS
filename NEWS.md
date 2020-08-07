# ICAMS 2.1.2.90xx
## Added
* Uploaded latest Windows and macOS binary package to GitHub so that users can download to their own computer and install ICAMS locally. See [README.md](https://github.com/steverozen/ICAMS/blob/master/README.md) for more details.
* Added documentation to *exported* functions `ReadAndSplitStrelkaSBSVCFs`,
`StrelkaSBSVCFFilesToCatalog`, `StrelkaSBSVCFFilesToCatalogAndPlotToPdf` and `StrelkaSBSVCFFilesToZipFile` informing the user that the function will find and merge adjacent SBS pairs into DBS if their VAFs are very similar. The default threshold value for VAF is 0.02.
* Added new exported data of catalog row order for SBS96, SBS1536 and DBS78 in SigProfiler format to `catalog.row.order.sp`.
* Created a new internal function `ConvertICAMSCatalogToSigProSBS96`.
* Added a new *exported* function `GetFreebayesVAF` for calculating variant allel frequencies from Freebayes VCF.
* Created two internal function `ReadVCF` and `ReadVCFs` with argument `variant.caller` that allows user to specify the variant caller that produced the VCFs. These two functions are able to read mixed VCFs.
* Added test data for Strelka mixed VCF (GRCh37).
* Added time zone information to file "run-information.txt" when calling
functions `MutectVCFFilesToZipFile`, `StrelkaSBSVCFFilesToZipFile` and
`StrelkaIDVCFFilesToZipFile`.
* Enabled "counts" -> "counts.signature" catalog transformation when the source catalog
has NULL abundance. But the target catalog should have the same `abundance`, `ref.genome`
and `region` attribute as the source catalog, otherwise the function will raise an error.
* Added an extra section "ID classification" in documentation of functions generating ID catalogs, `FindDelMH` and `FindMaxRepeatDel`.

## Changed
* **(Non backward-compatible)** Changed the return of functions
`MutectVCFFilesToCatalog`, `MutectVCFFilesToCatalogAndPlotToPdf` and
`MutectVCFFilesToZipFile`. The object `catID` in the returned list has been
changed from a matrix catalog to a **list**: 1st element is a matrix of the ID
(small insertion and deletion) catalog. 2nd element is a list of further
annotated VCFs (three additional columns `seq.context.width` , `seq.context` and `ID.class` are added to the original VCF).
* Minor changes to documentation of functions `PlotCatalog`, `PlotCatalogToPdf`, `StrelkaSBSVCFFilesToZipFile`, `StrelkaIDVCFFilesToZipFile` and `MutectVCFFilesToZipFile`.
* Updated documentation for the return value of functions `StrelkaIDVCFFilesToCatalog`, `StrelkaIDVCFFilesToCatalogAndPlotToPdf`, `StrelkaIDVCFFilesToZipFile` and `VCFsToIDCatalogs` to make it clearer to the user.
* Changed the return of functions `PlotCatalog` and `PlotCatalogToPdf` to an **invisible** list.

## Fixed
* Fixed bugs in `if` statement in internal functions
`GetCustomKmerCounts`、`GetStrandedKmerCounts` and `GetGenomeKmerCounts`.
* Fixed a bug in internal function `CreateOneColIDMatrix` when there is NA ID category.
* Fixed a bug in *exported* function `GetMutectVAF` for checking whether the
VCF is indeed a Mutect VCF.
* Fixed a bug in internal function `CreateOneColDBSMatrix` when the DBS VCF does not have any variant on the transcribed region.

<br/>

# ICAMS 2.1.2
## Added
* Created an internal function `MakeDataFrameFromVCF` to read in data lines of a VCF.
* New argument `name.of.VCF` in internal function `CheckAndFixChrNames` to make the error message more informative.
* New argument `name.of.VCF` in *exported* function `AnnotateIDVCF` to make the error message more informative.

## Changed
* Updated internal function `ReadStrelkaIDVCF` to make the error message more informative.
* **(Non backward-compatible)** Changed the return of *exported* function `AnnotateIDVCF` to a list. The first element `annotated.vcf` contains the annotated VCF. If there are rows that are discarded, the function will generate a warning and
a second element `discarded.variants` will be included in the returned list.

## Deprecated
* Argument `flag.mismatches` deprecated in *exported* function `AnnotateIDVCF`. If there are mismatches to references, the
function will automatically discard these rows. User can refer to the
element `discarded.variants` in the return value for the discarded variants.

## Fixed
* Fixed a bug in internal function `SplitStrelkaSBSVCF` when there are no non.SBS mutations in the input.
* Fixed a bug in internal function `MakeDataFrameFromMutectVCF` when a Mutect VCF has no meta-information lines.
* Fixed a bug in internal function `CreateOneColSBSMatrix` when an annotated SBS VCF has variants on transcribed regions that **all** fall on transcripts on **both** strand.
* Fixed a bug in internal function `CreateOneColDBSMatrix` when an annotated DBS VCF has variants on transcribed regions that **all** fall on transcripts on **both** strand.

<br/>

# ICAMS 2.1.1
## Added
* **(New)** Added columns of VAF (variant allel frequency) and read depth information
  to the split DBS.vcfs from merged SBSs when calling function
  `ReadAndSplitStrelkaSBSVCFs`.
* Added a new dependency package "zip" in ICAMS, to be used in three new *exported*  
  functions `MutectVCFFilesToZipFile`, `StrelkaSBSVCFFilesToZipFile` and 
  `StrelkaIDVCFFilesToZipFile`. 
* Added code to handle multiple alternate alleles in Strelka VCFs, e.g cases 
  where REF = T and ALT = A,C which indicates the alternate allele in some reads
  is A and in some reads C.
* Added an internal function to infer trans ranges from input ref.genome.
  Updated the documentation of the *exported* functions which have argument `trans.ranges`
  to make it optional.
* Added an additional argument `name.of.VCF` in internal functions 
  `ReadStrelkaSBSVCF`, `ReadStrelkaIDVCF` and *exported* function `GetStrelkaVAF`.
* Added an optional argument `flag.mismatches` in functions `VCFsToIDCatalogs`,
  `MutectVCFFilesToCatalog`, `MutectVCFFilesToCatalogAndPlotToPdf`,
  `MutectVCFFilesToZipFile`, `StrelkaIDVCFFilesToCatalog`,
  `StrelkaIDVCFFilesToCatalogAndPlotToPdf` and `StrelkaIDVCFFilesToZipFile`.
  
## Changed
* **(Non backward-compatible)** Changed the return of function `GetStrelkaVAF` and  
  `GetMutectVAF` to a data frame which contains the VAF and read depth information.
* **(Non backward-compatible)** Made the return from function `PlotCatalogToPdf` a
  list. The first element is a logical value indicating whether the plot is
  successful. The second element is a list containing the strand bias statistics 
  (only for SBS192Catalog with "counts" catalog.type
  and non-NULL abundance and argument `plot.SBS12` = TRUE).  
* Slight changes to behavior of `PlotCatalog` and `PlotCatalogToPdf`:
  For class SBS96Catalog: 
  **(New)** Allow setting ylim and cex.
  **(New)** For `PlotCatalog` (not `PlotCatalogToPdf`), allow plotting of a 96 x 2 catalog,
  in which case behavior is a stacked bar chart. 
  **(New)** Plot x axis tick marks if `xlabels` is not TRUE; set `par(tck = 0)` to suppress. 
  For class IndelCatalog:
  **(New)** Allow setting ylim.
* Renamed internal function creating custom abundance to `GetCustomKmerCounts`.
* Updated documentation for optional arguments in the functions reading VCF files.
  Used @inheritParams to reduce repetitive documentation.
* Updated function `PlotTransBiasGeneExpToPdf` so that ymax on the plot will be changed 
  based on `plot.type`.
* Changed the class of values in internal data `flat.abundance` from "numeric" to
  "integer".

## Deprecated
* Deprecated counts -> counts transformation in `TransformCatalog`; see documentation
  for rationale.

## Fixed
* Fixed a bug in function `TransformCatalog` and updated its documentation
  for parameter `target.abundance`.
* Fixed bugs in internal function `CheckAndFixChrNames` and updated the automated tests.
* Corrected error message in `TransformCatalog`.
* Fixed a bug in *exported* function `GetMutectVAF` and updated the warning
  message to make it more informative. 
* Fixed bugs in functions generating zip archive from VCF files and added 
  more automated test. 
* Updated methods of different types of catalogs for the generic function `cbind`
  to check the attributes of the incoming catalogs and assign attributes accordingly.
* Updated function `TransformCatalog` to check the attributes of the catalog to be 
  transformed in the first place.
* Fixed bugs in code removing complex indels.

<br/>

# ICAMS 2.0.10
## Added
* Added *exported* functions `AnnotateSBSVCF`, `AnnotateDBSVCF` and
  `AnnotateIDVCF`.
* Added two *exported* functions `PlotTransBiasGeneExp` and  `PlotTransBiasGeneExpToPdf`.
* Added an additional optional argument `names.of.VCFs` in functions
  `ReadAndSplitMutectVCFs`, `ReadAndSplitStrelkaSBSVCFs`, `ReadStrelkaIDVCFs`,
  `MutectVCFFilesToCatalog`, `MutectVCFFilesToCatalogAndPlotToPdf`,
  `StrelkaIDVCFFilesToCatalog`, `StrelkaIDVCFFilesToCatalogAndPlotToPdf`,
  `StrelkaSBSVCFFilesToCatalog` and `StrelkaSBSVCFFilesToCatalogAndPlotToPdf`
  for users to specify the names of samples in the VCF files.
* Added error checking on the order of rownames in the input of `as.catalog`.
* Added code to handle the case where the #CHROM lines appear
  interspersed with rows for variants in a Mutect VCF file.
* Added two *exported* package variables, `gene.expression.data.HepG2` and
  `gene.expression.data.MCF10A`.
* Added an additional optional argument `tumor.col.names` in functions
  `ReadAndSplitMutectVCFs`, `MutectVCFFilesToCatalog` and
  `MutectVCFFilesToCatalogAndPlotToPdf` to specify the column of the VCF
  that contains sequencing statistics such as sequencing depth; this column
  is often called "unknown" in Mutect.
* Added a "Comments" section in the documentation of functions     
  `MutectVCFFilesToCatalog`,
  `MutectVCFFilesToCatalogAndPlotToPdf`, `StrelkaSBSVCFFilesToCatalog`,
  `StrelkaSBSVCFFilesToCatalogAndPlotToPdf`, `VCFsToSBSCatalogs`,
  `VCFsToDBSCatalogs`, `ReadCatalog` informing the user how to change
  attributes of the generated catalog.
  
## Changed
* **(Non backward-compatible)**
  Made the return from functions `VCFsToIDCatalogs`, `StrelkaIDVCFFilesToCatalog`
  and `StrelkaIDVCFFilesToCatalogAndPlotToPdf` a list; 1st element is the
  spectrum catalog (previously the only return); 2nd element is a list of
  VCFs with additional annotations.
* **(Non backward-compatible)** 
  Made the return from function `PlotCatalog` a list. The first element is 
  a logical value indicating whether the plot is successful. The second element 
  is a numeric vector giving the coordinates of all the bar midpoints drawn,
  useful for adding to the graph (only implemented for SBS96Catalog).
* **(Slightly Non backward-compatible)** Changed the handling of the `output.file` argument in
  `MutectVCFFilesToCatalogAndPlotToPdf`, `StrelkaSBSVCFFilesToCatalogAndPlotToPdf`, 
  and `StrelkaIDVCFFilesToCatalogAndPlotToPdf`
  so that an indicator of the catalog type plus ".pdf" is simply
  appended to the base `output.file` name. Also made this argument
  optional with sensible default behavior.
* Changed column name 'gene.name' to 'gene.symbol' and added Ensembl gene ID 
  information in *exported* package variables `trans.ranges.GRCh37`, `trans.ranges.GRCh38`   and `trans.ranges.GRCm38`.
* In `FindDelMH`, cryptic repeats (i.e. un-normalized deletions in a repeat 
  such as GAGG deleted from CCCAGGGAGGGTCCC, which should be normalized
  to a deletion of AGGG) are now ignored with a warning rather than
  causing a `stop`.
* Updated plotting functions to make y axis labels more informative.
* Updated splitting functions for Mutect VCF to detect and move complex
  indels to "other" category.
  
## Fixed
* Updated methods so that they are compatible with new versions of
  R, in which the class of a matrix is c("matrix", "array").  
* Fixed a rarely encountered bug in `FindDelMH`, which previously did not flag the
  cryptic repeat in what is now the second example in the function documentation.

<br/>

# ICAMS 2.0.9
## Added
* `as.catalog` supports creation of the catalog from a vector (interpreted
  as a 1-column matrix) and optionally infers the class from the
  number of rows in the input.
* Exported some functions that are used in categorizing deletions as 
  a way to make documentation of the underlying logic more visible.

<br/>

# ICAMS 2.0.8
## Added
* Added ability to calculate VAF from newer version Mutect VCF files.

## Changed
* Moved two Bioconductor packages(BSgenome.Hsapiens.1000genomes.hs37d5, 
  BSgenome.Hsapiens.UCSC.hg38) from Imports to Suggests so that they can
  be used conditionally in "ICAMS" as requested by CRAN.
* Documentation updates.

<br/>

# ICAMS 2.0.7
## Fixed
* Additional corrections for submission to CRAN
  (putting on.exit(par(opar)) immediately after
  opar <- par(...)
  in every case.

<br/>

# ICAMS 2.0.6
## Fixed
* Additional changes for submission to CRAN (re-set
  graphics parameters to original values after plotting
  calls; documentation updates).

<br/>

# ICAMS 2.0.5
## Added
* Added support for BSgenome.Mmusculus.UCSC.mm10. 

## Changed
* Changes in response to initial submission to CRAN.

<br/>

# ICAMS 2.0
* First stable release. Mutational catalogs are S3 classes.

<br/>

# ICAMS 1.0
* Only available on GitHub; interfaces not compatible with ICAMS 2.0.
