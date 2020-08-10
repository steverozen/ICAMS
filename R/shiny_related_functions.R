#' Create a zip file which contains catalogs and plot PDFs from Strelka SBS VCF files
#'
#' Create 3 SBS catalogs (96, 192, 1536), 3 DBS catalogs (78, 136, 144) from the
#' Strelka SBS VCFs specified by \code{dir}, save the catalogs as CSV files,
#' plot them to PDF and generate a zip archive of all the output files. The
#' function will find and merge adjacent SBS pairs into DBS if their VAFs are
#' very similar. The default threshold value for VAF is 0.02.
#'
#' This function calls \code{\link{StrelkaSBSVCFFilesToCatalog}},
#' \code{\link{PlotCatalogToPdf}}, \code{\link{WriteCatalog}} and
#' \code{\link[zip]{zipr}}.
#'
#' @param dir Pathname of the directory which contains \strong{only} the Strelka
#'   SBS VCF files. Each Strelka SBS VCF \strong{must} have a file extension
#'   ".vcf" (case insensitive) and share the \strong{same} \code{ref.genome} and
#'   \code{region}.
#'   
#' @inheritParams MutectVCFFilesToZipFile
#'
#' @importFrom utils glob2rx
#' 
#' @importFrom zip zipr 
#'
#' @inheritSection StrelkaSBSVCFFilesToCatalogAndPlotToPdf Value
#'
#' @inheritSection StrelkaSBSVCFFilesToCatalogAndPlotToPdf Note
#' 
#' @inheritSection MutectVCFFilesToZipFile Comments
#' 
#' @export
#' 
#' @examples 
#' dir <- c(system.file("extdata/Strelka-SBS-vcf",
#'                      package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <- 
#'     StrelkaSBSVCFFilesToZipFile(dir, 
#'                                 zipfile = file.path(tempdir(), "test.zip"),
#'                                 ref.genome = "hg19", 
#'                                 trans.ranges = trans.ranges.GRCh37,
#'                                 region = "genome",
#'                                 base.filename = "Strelka-SBS")
#'   unlink(file.path(tempdir(), "test.zip"))}
StrelkaSBSVCFFilesToZipFile <- 
  function(dir,
           zipfile, 
           ref.genome, 
           trans.ranges = NULL, 
           region = "unknown", 
           names.of.VCFs = NULL, 
           base.filename = "",
           return.annotated.vcfs = FALSE,
           suppress.discarded.variants.warnings = TRUE) {
  files <- list.files(path = dir, pattern = "\\.vcf$", 
                      full.names = TRUE, ignore.case = TRUE)
  vcf.names <- basename(files)
  catalogs0 <- 
    StrelkaSBSVCFFilesToCatalog(files, ref.genome, trans.ranges, 
                                region, names.of.VCFs, 
                                return.annotated.vcfs, 
                                suppress.discarded.variants.warnings)
  mutation.loads <- GetMutationLoadsFromStrelkaSBSVCFs(catalogs0)
  strand.bias.statistics<- NULL
  
  # Retrieve the catalog matrix from catalogs0
  catalogs <- catalogs0
  catalogs$discarded.variants <- catalogs$annotated.vcfs <- NULL
  
  output.file <- ifelse(base.filename == "",
                        paste0(tempdir(), .Platform$file.sep),
                        file.path(tempdir(), paste0(base.filename, ".")))
  
  for (name in names(catalogs)) {
    WriteCatalog(catalogs[[name]],
                 file = paste0(output.file, name, ".csv"))
  }
  
  for (name in names(catalogs)) {
    PlotCatalogToPdf(catalogs[[name]],
                     file = paste0(output.file, name, ".pdf"))
    
    if (name == "catSBS192") {
      list <- PlotCatalogToPdf(catalogs[[name]],
                               file = paste0(output.file, "SBS12.pdf"),
                               plot.SBS12 = TRUE)
      strand.bias.statistics<- c(strand.bias.statistics, 
                                 list$strand.bias.statistics)
    }
  }
  
  zipfile.name <- basename(zipfile)
  AddRunInformation(files, vcf.names, zipfile.name, vcftype = "strelka.sbs",
                    ref.genome, region, mutation.loads, strand.bias.statistics)
  file.names <- list.files(path = tempdir(), pattern = "\\.(pdf|csv|txt)$", 
                           full.names = TRUE)
  zip::zipr(zipfile = zipfile, files = file.names)
  unlink(file.names)
  invisible(catalogs0)
}

#' Create a zip file which contains ID (small insertion and deletion) catalog
#' and plot PDF from Strelka ID VCF files
#' 
#' Create ID (small insertion and deletion) catalog from the Strelka ID VCFs
#' specified by \code{dir}, save the catalog as CSV file, plot it to PDF and
#' generate a zip archive of all the output files.
#'
#' This function calls \code{\link{StrelkaIDVCFFilesToCatalog}},
#' \code{\link{PlotCatalogToPdf}}, \code{\link{WriteCatalog}} and
#' \code{\link[zip]{zipr}}.
#' 
#' @inheritParams MutectVCFFilesToZipFile
#' 
#' @param dir Pathname of the directory which contains \strong{only} the Strelka
#'   ID VCF files. Each Strelka ID VCF \strong{must} have a file extension
#'   ".vcf" (case insensitive) and share the \strong{same} \code{ref.genome} and
#'   \code{region}.
#'   
#' @param base.filename Optional. The base name of the CSV and PDF file to be
#'   produced; the file is ending in \code{catID.csv} and \code{catID.pdf}
#'   respectively.
#'
#' @importFrom utils glob2rx
#' 
#' @importFrom zip zipr 
#' 
#' @inheritSection StrelkaIDVCFFilesToCatalog Value 
#' 
#' @inheritSection VCFsToIDCatalogs ID classification
#' 
#' @inheritSection VCFsToIDCatalogs Note
#'
#' @export
#' 
#' @examples  
#' dir <- c(system.file("extdata/Strelka-ID-vcf",
#'                      package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <- 
#'     StrelkaIDVCFFilesToZipFile(dir, 
#'                                zipfile = file.path(tempdir(), "test.zip"),
#'                                ref.genome = "hg19", 
#'                                region = "genome",
#'                                base.filename = "Strelka-ID")
#'   unlink(file.path(tempdir(), "test.zip"))}
StrelkaIDVCFFilesToZipFile <- 
  function(dir,
           zipfile, 
           ref.genome, 
           region = "unknown", 
           names.of.VCFs = NULL, 
           base.filename = "",
           flag.mismatches = 0,
           return.annotated.vcfs = FALSE,
           suppress.discarded.variants.warnings = TRUE) {
    files <- list.files(path = dir, pattern = "\\.vcf$", 
                        full.names = TRUE, ignore.case = TRUE)
    vcf.names <- basename(files)
    catalogs0 <-
      StrelkaIDVCFFilesToCatalog(files, ref.genome, region, names.of.VCFs,
                                 flag.mismatches, return.annotated.vcfs,
                                 suppress.discarded.variants.warnings)
    mutation.loads <- GetMutationLoadsFromStrelkaIDVCFs(catalogs0)
    strand.bias.statistics<- NULL
    
    # Retrieve the catalog matrix from catalogs0
    catalogs <- catalogs0
    catalogs$discarded.variants <- catalogs$annotated.vcfs <- NULL
    
    output.file <- ifelse(base.filename == "",
                          paste0(tempdir(), .Platform$file.sep),
                          file.path(tempdir(), paste0(base.filename, ".")))
    
    WriteCatalog(catalogs$catalog, 
                 file = paste0(output.file, "catID.csv"))
    
    PlotCatalogToPdf(catalogs$catalog, 
                     file = paste0(output.file, "catID.pdf"))
    
    zipfile.name <- basename(zipfile)
    AddRunInformation(files, vcf.names, zipfile.name, vcftype = "strelka.id",
                      ref.genome, region, mutation.loads, strand.bias.statistics)
    file.names <- list.files(path = tempdir(), pattern = "\\.(pdf|csv|txt)$", 
                             full.names = TRUE)
    zip::zipr(zipfile = zipfile, files = file.names)
    unlink(file.names)
    invisible(catalogs0)
  }

#' Create a zip file which contains catalogs and plot PDFs from Mutect VCF files
#'
#' Create 3 SBS catalogs (96, 192, 1536), 3 DBS catalogs (78, 136, 144) and
#' Indel catalog from the Mutect VCFs specified by \code{dir}, save the catalogs
#' as CSV files, plot them to PDF and generate a zip archive of all the output files.
#'
#' This function calls \code{\link{MutectVCFFilesToCatalog}},
#' \code{\link{PlotCatalogToPdf}}, \code{\link{WriteCatalog}} and
#' \code{\link[zip]{zipr}}.
#'
#' @param dir Pathname of the directory which contains \strong{only} the Mutect
#'   VCF files. Each Mutect VCF \strong{must} have a file extension ".vcf" (case
#'   insensitive) and share the \strong{same} \code{ref.genome} and
#'   \code{region}.
#'   
#' @param zipfile Pathname of the zip file to be created.    
#'   
#' @param ref.genome  A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}. 
#'   
#' @param trans.ranges Optional. If \code{ref.genome} specifies one of the
#'   \code{\link{BSgenome}} object 
#'   \enumerate{
#'     \item \code{\link[BSgenome.Hsapiens.1000genomes.hs37d5]{BSgenome.Hsapiens.1000genomes.hs37d5}}
#'     \item \code{\link[BSgenome.Hsapiens.UCSC.hg38]{BSgenome.Hsapiens.UCSC.hg38}}
#'     \item \code{\link[BSgenome.Mmusculus.UCSC.mm10]{BSgenome.Mmusculus.UCSC.mm10}}
#'   }
#'   then the function will infer \code{trans.ranges} automatically. Otherwise,
#'   user will need to provide the necessary \code{trans.ranges}. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.
#'   If \code{is.null(trans.ranges)} do not add transcript range
#'   information.
#'   
#' @param region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'  
#' @param names.of.VCFs Optional. Character vector of names of the VCF files.
#'   The order of names in \code{names.of.VCFs} should match the order of VCFs
#'   listed in \code{dir}. If \code{NULL}(default), this function will remove
#'   all of the path up to and including the last path separator (if any) in
#'   \code{dir} and file paths without extensions (and the leading dot) will be
#'   used as the names of the VCF files.
#'   
#' @param tumor.col.names Optional. Character vector of column names in VCFs which contain
#'   the tumor sample information. The order of names in \code{tumor.col.names}
#'   should match the order of VCFs listed in \code{dir}. If
#'   \code{tumor.col.names} is equal to \code{NA}(default), this function will
#'   use the 10th column in all the VCFs to calculate VAFs.
#'   See \code{\link{GetMutectVAF}} for more details.
#'   
#' @param base.filename Optional. The base name of the CSV and PDF files to be
#'   produced; multiple files will be generated, each ending in
#'   \eqn{x}\code{.csv} or \eqn{x}\code{.pdf}, where \eqn{x} indicates the type
#'   of catalog.
#'   
#' @param flag.mismatches Deprecated. If there are ID variants whose \code{REF}
#'   do not match the extracted sequence from \code{ref.genome}, the function
#'   will automatically discard these variants and an element
#'   \code{discarded.variants} will appear in the return value. See
#'   \code{\link{AnnotateIDVCF}} for more details.
#'   
#' @param return.annotated.vcfs Logical. Whether to return the annotated VCFs
#'   with additional columns showing mutation class for each variant. Default is
#'   FALSE.
#'   
#' @param suppress.discarded.variants.warnings Logical. Whether to suppress
#'   warning messages showing information about the discarded variants. Default
#'   is TRUE.
#'   
#' @importFrom utils glob2rx 
#' 
#' @importFrom zip zipr 
#'
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf Value
#' 
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf ID classification
#'
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf Note
#'   
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf Comments
#' 
#' @export
#' 
#' @examples 
#' dir <- c(system.file("extdata/Mutect-vcf",
#'                      package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <- 
#'     MutectVCFFilesToZipFile(dir, 
#'                             zipfile = file.path(tempdir(), "test.zip"),
#'                             ref.genome = "hg19", 
#'                             trans.ranges = trans.ranges.GRCh37,
#'                             region = "genome",
#'                             base.filename = "Mutect")
#'   unlink(file.path(tempdir(), "test.zip"))}
MutectVCFFilesToZipFile <- 
  function(dir,
           zipfile, 
           ref.genome, 
           trans.ranges = NULL, 
           region = "unknown", 
           names.of.VCFs = NULL, 
           tumor.col.names = NA,
           base.filename = "",
           flag.mismatches = 0,
           return.annotated.vcfs = FALSE,
           suppress.discarded.variants.warnings = TRUE) {
    files <- list.files(path = dir, pattern = "\\.vcf$", 
                        full.names = TRUE, ignore.case = TRUE)
    vcf.names <- basename(files)
    catalogs0 <- MutectVCFFilesToCatalog(files, ref.genome, trans.ranges, 
                                         region, names.of.VCFs, tumor.col.names,
                                         flag.mismatches, return.annotated.vcfs,
                                         suppress.discarded.variants.warnings)
    mutation.loads <- GetMutationLoadsFromMutectVCFs(catalogs0)
    strand.bias.statistics <- NULL
    
    # Retrieve the catalog matrix from catalogs0
    catalogs <- catalogs0
    catalogs$discarded.variants <- catalogs$annotated.vcfs <- NULL
    
    output.file <- ifelse(base.filename == "",
                          paste0(tempdir(), .Platform$file.sep),
                          file.path(tempdir(), paste0(base.filename, ".")))
    
    for (name in names(catalogs)) {
        WriteCatalog(catalogs[[name]],
                     file = paste0(output.file, name, ".csv"))
    }
    
    for (name in names(catalogs)) {
      PlotCatalogToPdf(catalogs[[name]],
                       file = paste0(output.file, name, ".pdf"))
      
      if (name == "catSBS192") {
        list <- PlotCatalogToPdf(catalogs[[name]],
                                 file = paste0(output.file, "SBS12.pdf"),
                                 plot.SBS12 = TRUE)
        strand.bias.statistics <- 
          c(strand.bias.statistics, list$strand.bias.statistics)
      }
    }
    
    zipfile.name <- basename(zipfile)
    AddRunInformation(files, vcf.names, zipfile.name, vcftype = "mutect", 
                      ref.genome, region, mutation.loads, strand.bias.statistics)
    file.names <- list.files(path = tempdir(), pattern = "\\.(pdf|csv|txt)$", 
                             full.names = TRUE)
    zip::zipr(zipfile = zipfile, files = file.names)
    unlink(file.names)
    invisible(catalogs0)
  }

#' @keywords internal
CombineAndReturnCatalogsForStrelkaSBSVCFs <-
  function(split.vcfs.list, SBS.list, DBS.list) {
    discarded.variants.list <- 
      list(SBS = SBS.list$discarded.variants,
           DBS = DBS.list$discarded.variants,
           ThreePlus = split.vcfs.list$ThreePlus,
           multiple.alt = split.vcfs.list$multiple.alt,
           not.analyzed = split.vcfs.list$not.analyzed)
    annotated.vcfs.list <- list(SBS = SBS.list$annotated.vcfs,
                                DBS = DBS.list$annotated.vcfs)
    # Remove NULL elements from the list
    discarded.variants.list2 <- Filter(Negate(is.null), discarded.variants.list)
    if (length(discarded.variants.list2) == 0) {
      discarded.variants.list2 <- NULL
    }
    annotated.vcfs.list2 <- Filter(Negate(is.null), annotated.vcfs.list)
    if (length(annotated.vcfs.list2) == 0) {
      annotated.vcfs.list2 <- NULL
    }
    
    combined.list <- list(catSBS96 = SBS.list$catSBS96,
                          catSBS192 = SBS.list$catSBS192,
                          catSBS1536 = SBS.list$catSBS1536,
                          catDBS78 = DBS.list$catDBS78,
                          catDBS136 = DBS.list$catDBS136,
                          catDBS144 = DBS.list$catDBS144,
                          discarded.variants = discarded.variants.list2,
                          annotated.vcfs = annotated.vcfs.list2)
    # Remove NULL elements from the list
    combined.list2 <- Filter(Negate(is.null), combined.list)
    return(combined.list2)
  }

#' Create SBS and DBS catalogs from Strelka SBS VCF files
#'
#' Create 3 SBS catalogs (96, 192, 1536) and 3 DBS catalogs (78, 136, 144) from
#' the Strelka SBS VCFs specified by \code{files}. The function will find and
#' merge adjacent SBS pairs into DBS if their VAFs are very similar. The default
#' threshold value for VAF is 0.02.
#' 
#' This function calls \code{\link{VCFsToSBSCatalogs}} and
#' \code{\link{VCFsToDBSCatalogs}}.
#'
#' @param files Character vector of file paths to the Strelka SBS VCF files.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#' 
#' @inheritSection StrelkaSBSVCFFilesToCatalogAndPlotToPdf Value
#'
#' @inheritSection StrelkaSBSVCFFilesToCatalogAndPlotToPdf Note
#' 
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf Comments
#' 
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-SBS-vcf",
#'                       "Strelka.SBS.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <- StrelkaSBSVCFFilesToCatalog(file, ref.genome = "hg19",
#'                                           trans.ranges = trans.ranges.GRCh37,
#'                                           region = "genome")}
StrelkaSBSVCFFilesToCatalog <-
  function(files, ref.genome, trans.ranges = NULL, region = "unknown", 
           names.of.VCFs = NULL, return.annotated.vcfs = FALSE, 
           suppress.discarded.variants.warnings = TRUE) {
    split.vcfs <- 
      ReadAndSplitStrelkaSBSVCFs(files, names.of.VCFs, 
                                 suppress.discarded.variants.warnings)
    SBS.list <- VCFsToSBSCatalogs(split.vcfs$SBS.vcfs, ref.genome, 
                                  trans.ranges, region, return.annotated.vcfs,
                                  suppress.discarded.variants.warnings)
    DBS.list <- VCFsToDBSCatalogs(split.vcfs$DBS.vcfs, ref.genome, 
                                  trans.ranges, region, return.annotated.vcfs,
                                  suppress.discarded.variants.warnings)
    CombineAndReturnCatalogsForStrelkaSBSVCFs(split.vcfs, SBS.list, DBS.list)
  }

#' Create ID (small insertion and deletion) catalog from Strelka ID VCF files
#'
#' Create ID (small insertion and deletion) catalog from the Strelka ID VCFs
#' specified by \code{files}
#' 
#' This function calls \code{\link{VCFsToIDCatalogs}}
#'
#' @param files Character vector of file paths to the Strelka ID VCF files.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#' 
#' @section Value:
#' A \strong{list} of elements:
#'   * \code{catalog}: The ID (small insertion and deletion) catalog with
#'   attributes added. See \code{\link{as.catalog}} for more details.
#' 
#' * \code{discarded.variants}: 
#' \strong{Only appearing when} there are variants that were excluded in the
#' analysis. 
#' A list of elements:
#'     + \code{ID}: ID variants discarded that can belong to the following
#'       categories:
#'         - Variants which have empty REF or ALT allels.
#'         - Variants whose REF allels do not match the extracted sequence from
#'          \code{ref.genome}.
#'         - Variants which cannot be categorized according to the canonical
#'         representation. See catalog.row.order$ID for the canonical
#'         representation.
#'     + \code{not.analyzed}: Variants discarded immediately after reading in
#'     the VCFs:
#'         - Duplicated "CHROM" and "POS" values.
#'         - Chromosome names that contain "#".
#'         - Chromosome names that contain "GL".
#'         - Chromosome names that contain "KI".
#'         - Chromosome names that contain "random".
#'         - Chromosome names that contain "Hs".
#'         - Chromosome names that contain "M".
#' 
#'   * \code{annotated.vcfs}: 
#' \strong{Only appearing when} \code{return.annotated.vcfs} = TRUE. A list of
#' data frames which contain the original VCF's ID mutation rows with three
#' additional columns \code{seq.context.width}, \code{seq.context} and
#' \code{ID.class} added. The category assignment of each ID mutation in VCF can
#' be obtained from \code{ID.class} column.
#' @md
#' 
#' @inheritSection VCFsToIDCatalogs ID classification
#' 
#' @inheritSection VCFsToIDCatalogs Note
#' 
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-ID-vcf",
#'                       "Strelka.ID.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catID <- StrelkaIDVCFFilesToCatalog(file, ref.genome = "hg19", 
#'                                       region = "genome")}
StrelkaIDVCFFilesToCatalog <- 
  function(files, ref.genome, region = "unknown", names.of.VCFs = NULL,
           flag.mismatches = 0, return.annotated.vcfs = FALSE,
           suppress.discarded.variants.warnings = TRUE) {
    vcfs <- ReadStrelkaIDVCFs(files, names.of.VCFs, 
                              suppress.discarded.variants.warnings)
    ID.list <- VCFsToIDCatalogs(vcfs$ID.vcfs, ref.genome, region, 
                                flag.mismatches, return.annotated.vcfs,
                                suppress.discarded.variants.warnings)
    
    discarded.variants.list <- 
      list(ID = ID.list$discarded.variants,
           not.analyzed = vcfs$discarded.variants)
    annotated.vcfs.list <- ID.list$annotated.vcfs
    # Remove NULL elements from the list
    discarded.variants.list2 <- Filter(Negate(is.null), discarded.variants.list)
    if (length(discarded.variants.list2) == 0) {
      discarded.variants.list2 <- NULL
    }
    annotated.vcfs.list2 <- Filter(Negate(is.null), annotated.vcfs.list)
    if (length(annotated.vcfs.list2) == 0) {
      annotated.vcfs.list2 <- NULL
    }
    
    combined.list <- list(catalog = ID.list$catalog,
                          discarded.variants = discarded.variants.list2,
                          annotated.vcfs = annotated.vcfs.list2)
    # Remove NULL elements from the list
    combined.list2 <- Filter(Negate(is.null), combined.list)
    return(combined.list2)
  }

#' @keywords internal
CombineAndReturnCatalogsForMutectVCFs <-
  function(split.vcfs.list, SBS.list, DBS.list, ID.list) {
    discarded.variants.list <- 
      list(SBS = SBS.list$discarded.variants,
           DBS = DBS.list$discarded.variants,
           ID = ID.list$discarded.variants,
           other.subs = split.vcfs.list$other.subs,
           multiple.alt = split.vcfs.list$multiple.alt,
           not.analyzed = split.vcfs.list$not.analyzed)
    annotated.vcfs.list <- list(SBS = SBS.list$annotated.vcfs,
                                DBS = DBS.list$annotated.vcfs,
                                ID = ID.list$annotated.vcfs)
    # Remove NULL elements from the list
    discarded.variants.list2 <- Filter(Negate(is.null), discarded.variants.list)
    if (length(discarded.variants.list2) == 0) {
      discarded.variants.list2 <- NULL
    }
    annotated.vcfs.list2 <- Filter(Negate(is.null), annotated.vcfs.list)
    if (length(annotated.vcfs.list2) == 0) {
      annotated.vcfs.list2 <- NULL
    }
    
    combined.list <- list(catSBS96 = SBS.list$catSBS96,
                          catSBS192 = SBS.list$catSBS192,
                          catSBS1536 = SBS.list$catSBS1536,
                          catDBS78 = DBS.list$catDBS78,
                          catDBS136 = DBS.list$catDBS136,
                          catDBS144 = DBS.list$catDBS144,
                          catID = ID.list$catalog,
                          discarded.variants = discarded.variants.list2,
                          annotated.vcfs = annotated.vcfs.list2)
    # Remove NULL elements from the list
    combined.list2 <- Filter(Negate(is.null), combined.list)
    return(combined.list2)
  }

#' Create SBS, DBS and Indel catalogs from Mutect VCF files
#'
#' Create 3 SBS catalogs (96, 192, 1536), 3 DBS catalogs (78, 136, 144) and
#' Indel catalog from the Mutect VCFs specified by \code{files}
#'
#' This function calls \code{\link{VCFsToSBSCatalogs}},
#' \code{\link{VCFsToDBSCatalogs}} and \code{\link{VCFsToIDCatalogs}}
#'
#' @param files Character vector of file paths to the Mutect VCF files.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#' 
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf Value
#' 
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf ID classification
#'
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf Note
#'   
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf Comments
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <- MutectVCFFilesToCatalog(file, ref.genome = "hg19", 
#'                                       trans.ranges = trans.ranges.GRCh37,
#'                                       region = "genome")}
MutectVCFFilesToCatalog <-
  function(files, ref.genome, trans.ranges = NULL, region = "unknown", 
           names.of.VCFs = NULL, tumor.col.names = NA, flag.mismatches = 0,
           return.annotated.vcfs = FALSE, 
           suppress.discarded.variants.warnings = TRUE) {
    split.vcfs <- 
      ReadAndSplitMutectVCFs(files, names.of.VCFs, tumor.col.names,
                             suppress.discarded.variants.warnings)
    
    SBS.list <- VCFsToSBSCatalogs(split.vcfs$SBS, ref.genome, 
                               trans.ranges, region, return.annotated.vcfs,
                               suppress.discarded.variants.warnings)
    DBS.list <- VCFsToDBSCatalogs(split.vcfs$DBS, ref.genome, 
                               trans.ranges, region, return.annotated.vcfs,
                               suppress.discarded.variants.warnings)
    ID.list <- VCFsToIDCatalogs(split.vcfs$ID, ref.genome, 
                                region, flag.mismatches, return.annotated.vcfs,
                                suppress.discarded.variants.warnings)
    CombineAndReturnCatalogsForMutectVCFs(split.vcfs, SBS.list, 
                                          DBS.list, ID.list)
  }

#' Read and split Strelka SBS VCF files
#' 
#' The function will find and merge adjacent SBS pairs into DBS if their VAFs
#' are very similar. The default threshold value for VAF is 0.02.
#'
#' @param files Character vector of file paths to the Strelka SBS VCF files.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#'   
#' @section Value: 
#' A list of elements as follows:
#'
#'    * \code{SBS.vcfs}: List of data.frames of pure SBS mutations -- no DBS or
#'    3+BS mutations.
#'
#'    * \code{DBS.vcfs}: List of data.frames of pure DBS mutations -- no SBS or
#'    3+BS mutations.
#'
#'    * \code{ThreePlus}: \strong{Only appearing when} there is list of
#'    data.tables with the key CHROM, LOW.POS, HIGH.POS which contain rows in
#'    the input that did not represent SBSs or DBSs.
#'
#'   * \code{multiple.alt}: \strong{Only appearing when} there is list of VCF
#'   like data.frames with rows for variants with multiple alternative alleles,
#'   for example ACA mutated to both AGA and ATA at the same position.
#'
#'   * \code{not.analyzed}: \strong{Only appearing when} there is list of VCF
#'   like data.frames with rows for variants that were discarded immediately
#'   after reading in the VCFs. The variants not analyzed can belong to the
#'   following categories:
#'       + Duplicated "CHROM" and "POS" values.
#'       + Chromosome names that contain "#".
#'       + Chromosome names that contain "GL".
#'       + Chromosome names that contain "KI".
#'       + Chromosome names that contain "random".
#'       * Chromosome names that contain "Hs".
#'       * Chromosome names that contain "M".
#' @md
#'
#' @seealso \code{\link{StrelkaSBSVCFFilesToCatalog}}
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-SBS-vcf",
#'                       "Strelka.SBS.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs(file)
ReadAndSplitStrelkaSBSVCFs <- 
  function(files, names.of.VCFs = NULL,
           suppress.discarded.variants.warnings = TRUE) {
    vcfs <- ReadStrelkaSBSVCFs(files, names.of.VCFs, 
                               suppress.discarded.variants.warnings)
    split.vcfs <- 
      SplitListOfStrelkaSBSVCFs(vcfs, 
                                suppress.discarded.variants.warnings)
    return(split.vcfs)
  }

#' Read Strelka ID (small insertion and deletion) VCF files
#'
#' @inheritParams ReadMutectVCFs
#'
#' @section Value: 
#' A list of elements:
#'    * \code{ID.vcfs}: List of data.frames containing ID mutations.
#'
#'   * \code{discarded.variants}: \strong{Only appearing when} there is list of VCF
#'   like data.frames with rows for variants that were discarded immediately
#'   after reading in the VCFs. The variants not analyzed can belong to the
#'   following categories:
#'       + Duplicated "CHROM" and "POS" values.
#'       + Chromosome names that contain "#".
#'       + Chromosome names that contain "GL".
#'       + Chromosome names that contain "KI".
#'       + Chromosome names that contain "random".
#'       * Chromosome names that contain "Hs".
#'       * Chromosome names that contain "M".
#' @md
#'
#' @inheritSection VCFsToIDCatalogs Note
#'
#' @seealso \code{\link{StrelkaIDVCFFilesToCatalog}}
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-ID-vcf",
#'                       "Strelka.ID.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadStrelkaIDVCFs(file)
ReadStrelkaIDVCFs <- function(files, names.of.VCFs = NULL, 
                              suppress.discarded.variants.warnings = TRUE) {
  vcfs <- 
    lapply(files, FUN = ReadStrelkaIDVCF, name.of.VCF = names.of.VCFs,
           suppress.discarded.variants.warnings = 
             suppress.discarded.variants.warnings)
  if (is.null(names.of.VCFs)) {
    names(vcfs) <- tools::file_path_sans_ext(basename(files))
  } else {
    # Check whether the number of VCFs match the number of names
    # in names.of.VCFs
    CheckNamesOfVCFs(files, names.of.VCFs)
    names(vcfs) <- names.of.VCFs
  }
  
  list.of.ID.vcfs <- lapply(vcfs, function(f1) f1$df)
  list.of.discarded.variants <- 
    lapply(vcfs, function(f2) f2$discarded.variants)
  # Remove NULL elements from the list
  list.of.discarded.variants2 <- Filter(Negate(is.null), list.of.discarded.variants)
  
  if (length(list.of.discarded.variants2) == 0) {
    return(list(ID.vcfs = list.of.ID.vcfs))
  } else {
    return(list(ID.vcfs = list.of.ID.vcfs, 
                discarded.variants = list.of.discarded.variants2))
  }
}

#' Read and split Mutect VCF files
#'
#' @param files Character vector of file paths to the Mutect VCF files.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#'   
#' @section Value: A list containing the following objects:
#'
#'   * \code{SBS}: List of VCFs with only single base substitutions.
#'
#'
#'   * \code{DBS}: List of VCFs with only doublet base substitutions as called
#'   by Mutect.
#'
#'   * \code{ID}: List of VCFs with only small insertions and deletions.
#'
#'   * \code{other.subs}: \strong{Only appearing when} there is list of VCF
#'   like data.frames with rows for coordinate substitutions involving 3 or more
#'   nucleotides (e.g. ACT > TGA or AACT > GGTA) and rows for complex indels.
#'
#'   * \code{multiple.alt}: \strong{Only appearing when} there is list of VCF
#'   like data.frames with rows for variants with multiple alternative alleles,
#'   for example ACA mutated to both AGA and ATA at the same position.
#'
#'   * \code{not.analyzed}: \strong{Only appearing when} there is list of VCF
#'   like data.frames with rows for variants that were discarded immediately
#'   after reading in the VCFs. The variants not analyzed can belong to the
#'   following categories:
#'       + Duplicated "CHROM" and "POS" values.
#'       + Chromosome names that contain "#".
#'       + Chromosome names that contain "GL".
#'       + Chromosome names that contain "KI".
#'       + Chromosome names that contain "random".
#'       * Chromosome names that contain "Hs".
#'       * Chromosome names that contain "M".
#' @md
#' 
#' @seealso \code{\link{MutectVCFFilesToCatalog}}
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitMutectVCFs(file)
ReadAndSplitMutectVCFs <- 
  function(files, names.of.VCFs = NULL, tumor.col.names = NA,
           suppress.discarded.variants.warnings = TRUE) {
    vcfs <- ReadMutectVCFs(files, names.of.VCFs, tumor.col.names,
                           suppress.discarded.variants.warnings)
    split.vcfs <- 
      SplitListOfMutectVCFs(vcfs, suppress.discarded.variants.warnings)
    return(split.vcfs)
  }

#' Check and return SBS catalogs
#'
#' @param catSBS96 An SBS96 catalog.
#' 
#' @param catSBS1536 An SBS1536 catalog.
#' 
#' @param catSBS192 An SBS192 catalog.
#' 
#' @param discarded.variants A list of discarded variants.
#' 
#' @param annotated.vcfs A list of annotated VCFs.
#'
#' @return A list of SBS catalogs. Also return the discarded variants and
#'   annotated VCFs if they exit.
#' 
#' @keywords internal
CheckAndReturnSBSCatalogs <- 
  function(catSBS96, catSBS1536, catSBS192 = NULL, discarded.variants,
           annotated.vcfs) {
    if (is.null(catSBS192)) {
      if (length(discarded.variants) == 0) {
        if (length(annotated.vcfs) == 0) {
          return(list(catSBS96 = catSBS96, catSBS1536 = catSBS1536))
        } else {
          return(list(catSBS96 = catSBS96, catSBS1536 = catSBS1536,
                      annotated.vcfs = annotated.vcfs))
        }
      } else {
        if (length(annotated.vcfs) == 0) {
          return(list(catSBS96 = catSBS96, catSBS1536 = catSBS1536,
                      discarded.variants = discarded.variants))
        } else {
          return(list(catSBS96 = catSBS96, catSBS1536 = catSBS1536,
                      discarded.variants = discarded.variants,
                      annotated.vcfs = annotated.vcfs))
        }
      }
    } else {
      if (length(discarded.variants) == 0) {
        if (length(annotated.vcfs) == 0) {
          return(list(catSBS96 = catSBS96, catSBS192 = catSBS192,
                      catSBS1536 = catSBS1536))
        } else {
          return(list(catSBS96 = catSBS96, catSBS192 = catSBS192,
                      catSBS1536 = catSBS1536,
                      annotated.vcfs = annotated.vcfs))
        }
      } else {
        if (length(annotated.vcfs) == 0) {
          return(list(catSBS96 = catSBS96, catSBS192 = catSBS192,
                      catSBS1536 = catSBS1536,
                      discarded.variants = discarded.variants))
        } else {
          return(list(catSBS96 = catSBS96, catSBS192 = catSBS192,
                      catSBS1536 = catSBS1536,
                      discarded.variants = discarded.variants,
                      annotated.vcfs = annotated.vcfs))
        }
      }
    }
  }

#' Create SBS catalogs from SBS VCFs
#'
#' Create a list of 3 catalogs (one each for 96, 192, 1536)
#' out of the contents in list.of.SBS.vcfs. The SBS VCFs must not contain
#' DBSs, indels, or other types of mutations.
#'
#' @param list.of.SBS.vcfs List of in-memory data frames of pure SBS mutations
#'   -- no DBS or 3+BS mutations. The list names will be the sample ids in the
#'   output catalog.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#' 
#' @section Value:  
#' A list containing the following objects:
#' 
#' * \code{catSBS96}, \code{catSBS192}, \code{catSBS1536}: Matrix of 
#' 3 SBS catalogs (one each for 96, 192, and 1536).
#' 
#' * \code{discarded.variants}: \strong{Only appearing when} there are SBS
#' variants whose pentanucleotide context contains "N".
#' 
#' * \code{annotated.vcfs}: 
#' \strong{Only appearing when} \code{return.annotated.vcfs} = TRUE.
#'     SBS VCF annotated by \code{\link{AnnotateSBSVCF}} with
#'     three new columns \code{SBS96.class}, \code{SBS192.class} and
#'     \code{SBS1536.class} showing the mutation class for each SBS variant.
#' 
#' If \code{trans.ranges} is not provided by user and cannot be inferred by
#' ICAMS, SBS 192 catalog will not be generated. Each catalog has attributes
#' added. See \code{\link{as.catalog}} for more details.
#' @md
#'
#' @note SBS 192 catalogs only contain mutations in transcribed regions.
#'
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf Comments
#' 
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' list.of.SBS.vcfs <- ReadAndSplitMutectVCFs(file)$SBS
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs.SBS <- VCFsToSBSCatalogs(list.of.SBS.vcfs, ref.genome = "hg19",
#'                                     trans.ranges = trans.ranges.GRCh37,
#'                                     region = "genome")}
VCFsToSBSCatalogs <- function(list.of.SBS.vcfs, ref.genome, 
                              trans.ranges = NULL, region = "unknown",
                              return.annotated.vcfs = FALSE,
                              suppress.discarded.variants.warnings = TRUE) {
  ncol <- length(list.of.SBS.vcfs)
  
  catSBS96 <- empty.cats$catSBS96
  catSBS192 <- empty.cats$catSBS192
  catSBS1536 <- empty.cats$catSBS1536
  trans.ranges <- InferTransRanges(ref.genome, trans.ranges)
  
  annotated.vcfs <- discarded.variants <- list()
  vcf.names <- character()
  
  for (i in 1:ncol) {
    SBS.vcf <- list.of.SBS.vcfs[[i]]
    sample.id <- names(list.of.SBS.vcfs)[i]
    annotated.SBS.vcf <- AnnotateSBSVCF(SBS.vcf, ref.genome, trans.ranges)
    if (suppress.discarded.variants.warnings == TRUE) {
      SBS.cat <- 
        suppressWarnings(CreateOneColSBSMatrix(annotated.SBS.vcf, sample.id,
                                               return.annotated.vcfs))  
    } else {
      SBS.cat <- CreateOneColSBSMatrix(annotated.SBS.vcf, sample.id,
                                       return.annotated.vcfs)
    }
    catSBS96 <- cbind(catSBS96, SBS.cat$catSBS96)
    if (!is.null(trans.ranges)) {
      catSBS192 <- cbind(catSBS192, SBS.cat$catSBS192)
    }
    catSBS1536 <- cbind(catSBS1536, SBS.cat$catSBS1536)
    if (return.annotated.vcfs == TRUE) {
      annotated.vcfs <- c(annotated.vcfs, list(SBS.cat$annotated.vcf))
    }
    if (!is.null(SBS.cat$discarded.variants)) {
      discarded.variants <- 
        c(discarded.variants, list(SBS.cat$discarded.variants))
      vcf.names <- c(vcf.names, sample.id)
    }
  }
  
  colnames(catSBS96) <- colnames(catSBS1536) <- names(list.of.SBS.vcfs)
  catSBS96 <- as.catalog(catSBS96, ref.genome = ref.genome,
                         region = region, catalog.type = "counts")
  catSBS1536 <- as.catalog(catSBS1536, ref.genome = ref.genome,
                           region = region, catalog.type = "counts")
  if (return.annotated.vcfs == TRUE) {
    names(annotated.vcfs) <- names(list.of.SBS.vcfs)
  }
  names(discarded.variants) <- vcf.names
  
  if (is.null(trans.ranges)) {
    CheckAndReturnSBSCatalogs(catSBS96, catSBS1536, discarded.variants,
                              annotated.vcfs)
  }
  
  colnames(catSBS192) <- names(list.of.SBS.vcfs)
  in.transcript.region <- ifelse(region == "genome", "transcript", region)
  catSBS192 <- as.catalog(catSBS192, ref.genome = ref.genome,
                          region = in.transcript.region, 
                          catalog.type = "counts")
  
  CheckAndReturnSBSCatalogs(catSBS96, catSBS1536, catSBS192, 
                            discarded.variants, annotated.vcfs)
}

#' Check and return DBS catalogs
#'
#' @param catDBS78 An DBS78 catalog.
#' 
#' @param catDBS136 An DBS136 catalog.
#' 
#' @param catDBS144 An DBS144 catalog.
#' 
#' @param discarded.variants A list of discarded variants.
#' 
#' @param annotated.vcfs A list of annotated VCFs.
#'
#' @return A list of DBS catalogs. Also return the discarded variants and
#'   annotated VCFs if they exit.
#' 
#' @keywords internal
CheckAndReturnDBSCatalogs <- 
  function(catDBS78, catDBS136, catDBS144 = NULL, discarded.variants,
           annotated.vcfs) {
    if (is.null(catDBS144)) {
      if (length(discarded.variants) == 0) {
        if (length(annotated.vcfs) == 0) {
          return(list(catDBS78 = catDBS78, catDBS136 = catDBS136))
        } else {
          return(list(catDBS78 = catDBS78, catDBS136 = catDBS136,
                      annotated.vcfs = annotated.vcfs))
        }
      } else {
        if (length(annotated.vcfs) == 0) {
          return(list(catDBS78 = catDBS78, catDBS136 = catDBS136,
                      discarded.variants = discarded.variants))
        } else {
          return(list(catDBS78 = catDBS78, catDBS136 = catDBS136,
                      discarded.variants = discarded.variants,
                      annotated.vcfs = annotated.vcfs))
        }
      }
    } else {
      if (length(discarded.variants) == 0) {
        if (length(annotated.vcfs) == 0) {
          return(list(catDBS78 = catDBS78, catDBS136 = catDBS136,
                      catDBS144 = catDBS144))
        } else {
          return(list(catDBS78 = catDBS78, catDBS136 = catDBS136,
                      catDBS144 = catDBS144,
                      annotated.vcfs = annotated.vcfs))
        }
      } else {
        if (length(annotated.vcfs) == 0) {
          return(list(catDBS78 = catDBS78, catDBS136 = catDBS136,
                      catDBS144 = catDBS144,
                      discarded.variants = discarded.variants))
        } else {
          return(list(catDBS78 = catDBS78, catDBS136 = catDBS136,
                      catDBS144 = catDBS144,
                      discarded.variants = discarded.variants,
                      annotated.vcfs = annotated.vcfs))
        }
      }
    }
  }

#' Create DBS catalogs from VCFs
#'
#' Create a list of 3 catalogs (one each for DBS78, DBS144 and DBS136)
#' out of the contents in list.of.DBS.vcfs. The VCFs must not contain
#' any type of mutation other then DBSs.
#'
#' @param list.of.DBS.vcfs List of in-memory data frames of pure DBS mutations
#'   -- no SBS or 3+BS mutations. The list names will be the sample ids in the
#'   output catalog.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#'
#' @section Value:  
#' A list containing the following objects:
#' 
#' * \code{catDBS78}, \code{catDBS136}, \code{catDBS144}: Matrix of
#' 3 DBS catalogs (one each for 78, 136, and 144).
#'
#' * \code{discarded.variants}: \strong{Only appearing when} there are DBS
#' variants whose tetranucleotide context contains "N".
#'   
#' * \code{annotated.vcfs}: \strong{Only appearing when}
#' \code{return.annotated.vcfs} = TRUE. DBS VCF annotated by
#' \code{\link{AnnotateDBSVCF}} with three new columns \code{DBS78.class},
#' \code{DBS136.class} and \code{DBS144.class} showing the mutation class for
#' each DBS variant.
#' 
#' If \code{trans.ranges} is not provided by user and cannot be inferred by
#' ICAMS, DBS 144 catalog will not be generated. Each catalog has
#' attributes added. See \code{\link{as.catalog}} for more details.
#' @md
#'
#' @note DBS 144 catalog only contains mutations in transcribed regions.
#'
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf Comments
#' 
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' list.of.DBS.vcfs <- ReadAndSplitMutectVCFs(file)$DBS
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs.DBS <- VCFsToDBSCatalogs(list.of.DBS.vcfs, ref.genome = "hg19",
#'                                     trans.ranges = trans.ranges.GRCh37,
#'                                     region = "genome")}
VCFsToDBSCatalogs <- function(list.of.DBS.vcfs, ref.genome, 
                              trans.ranges = NULL, region = "unknown",
                              return.annotated.vcfs = FALSE,
                              suppress.discarded.variants.warnings = TRUE) {
  ncol <- length(list.of.DBS.vcfs)
  
  catDBS78 <- empty.cats$catDBS78
  catDBS136 <- empty.cats$catDBS136
  catDBS144 <- empty.cats$catDBS144
  trans.ranges <- InferTransRanges(ref.genome, trans.ranges)
  
  annotated.vcfs <- discarded.variants <- list()
  vcf.names <- character()
  
  for (i in 1 : ncol) {
    DBS.vcf <- list.of.DBS.vcfs[[i]]
    sample.id <- names(list.of.DBS.vcfs)[i]
    annotated.DBS.vcf <- AnnotateDBSVCF(DBS.vcf, ref.genome, trans.ranges)
    if (suppress.discarded.variants.warnings == TRUE) {
      DBS.cat <- 
        suppressWarnings(CreateOneColDBSMatrix(annotated.DBS.vcf, sample.id,
                                               return.annotated.vcfs))
    } else {
      DBS.cat <- CreateOneColDBSMatrix(annotated.DBS.vcf, sample.id,
                                       return.annotated.vcfs)
    }
    catDBS78 <- cbind(catDBS78, DBS.cat$catDBS78)
    catDBS136 <- cbind(catDBS136, DBS.cat$catDBS136)
    if (!is.null(trans.ranges)) {
      catDBS144 <- cbind(catDBS144, DBS.cat$catDBS144)
    }
    if (return.annotated.vcfs == TRUE) {
      annotated.vcfs <- c(annotated.vcfs, list(DBS.cat$annotated.vcf))
    }
    if (!is.null(DBS.cat$discarded.variants)) {
      discarded.variants <- 
        c(discarded.variants, list(DBS.cat$discarded.variants))
      vcf.names <- c(vcf.names, sample.id)
    }
  }
  
  colnames(catDBS78) <- colnames(catDBS136) <- names(list.of.DBS.vcfs)
  catDBS78 <- as.catalog(catDBS78, ref.genome = ref.genome,
                         region = region, catalog.type = "counts")
  catDBS136 <- as.catalog(catDBS136, ref.genome = ref.genome,
                          region = region, catalog.type = "counts")
  if (return.annotated.vcfs == TRUE) {
    names(annotated.vcfs) <- names(list.of.DBS.vcfs)
  }
  names(discarded.variants) <- vcf.names
  
  if (is.null(trans.ranges)) {
    CheckAndReturnDBSCatalogs(catDBS78, catDBS136, discarded.variants,
                              annotated.vcfs)
  }
  
  colnames(catDBS144) <- names(list.of.DBS.vcfs)
  in.transcript.region <- ifelse(region == "genome", "transcript", region)
  catDBS144 <- as.catalog(catDBS144, ref.genome = ref.genome,
               region = in.transcript.region, catalog.type = "counts")
  
  CheckAndReturnDBSCatalogs(catDBS78, catDBS136, catDBS144,
                            discarded.variants, annotated.vcfs)
}

#' Check and return ID catalog 
#'
#' @param catID An ID catalog.
#' 
#' @param discarded.variants A list of discarded variants.
#' 
#' @param annotated.vcfs A list of annotated VCFs.
#'
#' @return A list of ID catalog. Also return the discarded variants and
#'   annotated VCFs if they exit.
#'   
#' @keywords internal
CheckAndReturnIDCatalog <- 
  function(catID, discarded.variants, annotated.vcfs) {
    if (length(discarded.variants) == 0) {
      if (length(annotated.vcfs) == 0) {
        return(list(catalog = catID))
      } else {
        return(list(catalog = catID, annotated.vcfs = annotated.vcfs))
      }
    } else {
      if (length(annotated.vcfs) == 0) {
        return(list(catalog = catID, discarded.variants = discarded.variants))
      } else {
        return(list(catalog = catID, discarded.variants = discarded.variants,
                    annotated.vcfs = annotated.vcfs))
      }
    }
  }

#' Create ID (small insertion and deletion) catalog from ID VCFs
#'
#' @param list.of.vcfs List of in-memory ID VCFs. The list names will be
#' the sample ids in the output catalog.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region A character string acting as a region identifier, one of
#' "genome", "exome".
#' 
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#'   
#' @section Value:
#' A \strong{list} of elements:
#'   * \code{catalog}: The ID (small insertion and deletion) catalog with
#'   attributes added. See \code{\link{as.catalog}} for more details.
#' 
#'   * \code{discarded.variants}: 
#' \strong{Only appearing when} there are ID variants that were discarded. A
#' list of data frames which contain the discarded variants from the original
#' VCF. The discarded variants can belong to the following categories:
#'       + Variants which have the same number of bases for REF and ALT alleles.
#'       + Variants which have empty REF or ALT allels.
#'       + Complex indels.
#'       + Variants whose \code{REF} do not match the extracted sequence from
#'         \code{ref.genome}.
#'       + Variants which cannot be categorized according to the canonical
#'       representation. See catalog.row.order$ID for the canonical
#'       representation.
#' 
#'   * \code{annotated.vcfs}: 
#' \strong{Only appearing when} \code{return.annotated.vcfs} = TRUE. A list of
#' data frames which contain the original VCF's ID mutation rows with three
#' additional columns \code{seq.context.width}, \code{seq.context} and
#' \code{ID.class} added. The category assignment of each ID mutation in VCF can
#' be obtained from \code{ID.class} column.
#' @md
#' 
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf ID classification
#'   
#' @section Note:
#'  In ID (small insertion and deletion) catalogs, deletion repeat sizes range
#'  from 0 to 5+, but for plotting and end-user documentation deletion repeat
#'  sizes range from 1 to 6+.
#'   
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-ID-vcf/",
#'                       "Strelka.ID.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' list.of.ID.vcfs <- ReadStrelkaIDVCFs(file)$ID.vcfs                      
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5",
#'  quietly = TRUE)) {
#'   catID <- VCFsToIDCatalogs(list.of.ID.vcfs, ref.genome = "hg19",
#'                             region = "genome")}
VCFsToIDCatalogs <- function(list.of.vcfs, ref.genome, region = "unknown",
                             flag.mismatches = 0,
                             return.annotated.vcfs = FALSE,
                             suppress.discarded.variants.warnings = TRUE) {
  ncol <- length(list.of.vcfs)
  
  # Create a 0-column matrix with the correct row labels.
  catID <- matrix(0, nrow = length(ICAMS::catalog.row.order$ID), ncol = 0)
  rownames(catID) <- ICAMS::catalog.row.order$ID
  annotated.vcfs <- discarded.variants <- list()
  discarded.variants.vcf.name <- character()
  
  vcf.names <- names(list.of.vcfs)
  for (i in 1:ncol) {
    ID <- list.of.vcfs[[i]]
    sample.id <- names(list.of.vcfs)[i]
    if (suppress.discarded.variants.warnings == TRUE) {
      list <- 
        suppressWarnings(AnnotateIDVCF(ID, ref.genome = ref.genome,
                                       flag.mismatches = flag.mismatches,
                                       name.of.VCF = vcf.names[i]))
    } else {
      list <- AnnotateIDVCF(ID, ref.genome = ref.genome,
                            flag.mismatches = flag.mismatches,
                            name.of.VCF = vcf.names[i])
    }

    # Create an empty data frame for discarded variants
    df <- ID[0, ]
    
    if (!is.null(list$discarded.variants)) {
      df <- dplyr::bind_rows(df, list$discarded.variants)
    }
    # Unlike the case for SBS and DBS, we do not add transcript information.
    if (suppress.discarded.variants.warnings == TRUE) {
      tmp <- suppressWarnings({
        CreateOneColIDMatrix(list$annotated.vcf, 
                             sample.id = sample.id, 
                             return.annotated.vcf = return.annotated.vcfs)
      })
    } else {
      tmp <- CreateOneColIDMatrix(list$annotated.vcf, sample.id = sample.id, 
                                  return.annotated.vcf = return.annotated.vcfs)
    }
    one.ID.column <- tmp$catalog
    rm(ID)
    catID <- cbind(catID, one.ID.column)
    
    if (return.annotated.vcfs == TRUE) {
      annotated.vcfs <- c(annotated.vcfs, list(tmp$annotated.vcf))
    }
    if (!is.null(tmp$discarded.variants)) {
      df <- dplyr::bind_rows(df, tmp$discarded.variants)
    }
    if (nrow(df) != 0) {
      discarded.variants <- c(discarded.variants, list(df))
      discarded.variants.vcf.name <- 
        c(discarded.variants.vcf.name, vcf.names[i])
    }
  }
  
  colnames(catID) <- names(list.of.vcfs)
  catID <- as.catalog(catID, ref.genome = ref.genome,
                      region = region, catalog.type = "counts")
  if (return.annotated.vcfs == TRUE) {
    names(annotated.vcfs) <- names(list.of.vcfs)
  }
  names(discarded.variants) <- discarded.variants.vcf.name
  CheckAndReturnIDCatalog(catID, discarded.variants, annotated.vcfs)
}

#' Calculate the number of space needed to add strand bias statistics to
#' the run-information.txt file.
#' 
#' @param list A list containing strand bias statistics.
#'
#' @return A matrix containing the space information.
#' 
#' @keywords internal
CalculateNumberOfSpace <- function(list) {
  mat <- matrix(0, nrow = 6, ncol = 3)
  rownames(mat) <- rownames(list[[1]])
  colnames(mat) <- c("space.counts", "space.q.value", "space.total")
  
  GetInformation <- function(df, info) {
    if (info == "counts") {
      df1 <- df[, c("transcribed", "untranscribed")]
      df1$max.counts <- apply(df1, MARGIN = 1, FUN = max)
      df1$counts.space <- nchar(df1$max.counts)
      return(df1[, "counts.space", drop = FALSE])
    } else if (info == "q.value") {
      df$q.values.sci <- formatC(df$q.values, format = "e", digits = 2)
      df$q.values.space <- nchar(df$q.values.sci)
      return(df[, "q.values.space", drop = FALSE])
    }
  }
  list.counts.space <- lapply(list, FUN = GetInformation, info = "counts")
  list.q.values.space <- lapply(list, FUN = GetInformation, info = "q.value")
  df.counts.space <- do.call("cbind", list.counts.space)
  df.q.values.space <- do.call("cbind", list.q.values.space)
  
  mat[, "space.counts"] <- 
    pmax(apply(df.counts.space, MARGIN = 1, FUN = max), nchar("counts"))
  mat[, "space.q.value"] <- 
    pmax(apply(df.q.values.space, MARGIN = 1, FUN = max), nchar("p-value"))
  mat[, "space.total"] <- rowSums(mat) + 3
  return(mat)
}

#' Create a run information text file from generating zip archive from VCF
#' files.
#' 
#' @importFrom  stringi stri_pad
#' 
#' @importFrom tools md5sum
#' 
#' @importFrom utils packageVersion
#' 
#' @keywords internal
AddRunInformation <- 
  function(files, vcf.names, zipfile.name, vcftype, ref.genome, 
           region, mutation.loads, strand.bias.statistics) {
    
    run.info <- 
      file(description = file.path(tempdir(), "run-information.txt"), open = "w")
    
    # Add the header information
    time.info <- strftime(Sys.time(), usetz = TRUE) # Get time zone information
    time.info1 <- 
      gsub(pattern = "+", replacement = "UTC+", x = time.info, fixed = TRUE)
    header <- paste0("run-information.txt file for ", zipfile.name, 
                     " created on ", time.info1)
    char.length <- nchar(header)
    writeLines(paste(rep("-", char.length), collapse = ""), run.info)
    writeLines(header, run.info)
    writeLines(paste(rep("-", char.length), collapse = ""), run.info)
    
    # Add section on purpose of ICAMS software
    writeLines("", run.info)
    writeLines("--- About ICAMS ---", run.info)
    writeLines(c("Analysis and visualization of experimentally elucidated mutational",
                 "signatures - the kind of analysis and visualization in Boot et al.,",
                 "'In-depth characterization of the cisplatin mutational signature in",
                 "human cell lines and in esophageal and liver tumors', ", 
                 "Genome Research 2018, https://doi.org/10.1101/gr.230219.117.",
                 "'ICAMS' stands for In-depth Characterization and Analysis of",
                 "Mutational Signatures. 'ICAMS' has functions to read in variant",
                 "call files (VCFs) and to collate the corresponding catalogs of",
                 "mutational spectra and to analyze and plot catalogs of mutational",
                 'spectra and signatures. Handles both "counts-based" and ', 
                 '"density-based" catalogs of mutational spectra or signatures.'), 
               run.info)
    writeLines("", run.info)
    writeLines(c("For complete documentation of ICAMS, please refer to ",
                 "https://cran.rstudio.com/web/packages/ICAMS/index.html"), run.info)
    
    # Add ICAMS and R version used
    writeLines("", run.info)
    writeLines("--- Version of the software ---", run.info)
    writeLines(paste0("ICAMS version: ", packageVersion("ICAMS")), run.info)
    writeLines(paste0("R version:     ", getRversion()), run.info)
    
    # Add input parameters specified by the user
    writeLines("", run.info)
    writeLines("--- Input parameters ---", run.info)
    if (vcftype == "strelka.sbs") {
      vcftype <- "Strelka SBS VCF"
    } else if (vcftype == "strelka.id") {
      vcftype <- "Strelka ID VCF"
    } else if (vcftype == "mutect") {
      vcftype <- "Mutect VCF"
    }
    
    if (ref.genome == "hg19") {
      ref.genome <- "Human GRCh37/hg19"
    } else if (ref.genome == "hg38") {
      ref.genome <- "Human GRCh38/hg38"
    } else if (ref.genome == "mm10") {
      ref.genome <- "Mouse GRCm38/mm10"
    }
    writeLines(paste0("Type of VCF:      ", vcftype), run.info)
    writeLines(paste0("Reference genome: ", ref.genome), run.info)
    writeLines(paste0("Region:           ", region), run.info)
    
    # Add input files information
    writeLines("", run.info)
    writeLines("--- Input files ---", run.info)
    max.num.of.char <- max(nchar(vcf.names))
    # Add a description of the information listed for input files
    writeLines(paste0(stri_pad("Name", width = max.num.of.char,
                               side = "right"), "  ", 
                      "# of data lines", "  ",
                      stri_pad("MD5", width = 32,
                               side = "right"), "  ",
                      "# of SBS", "  ",
                      "# of DBS", "  ",
                      "# of ID", "  ",
                      "# of discarded variants*", "  "),
               run.info)
    
    num.of.file <- length(files)
    
    for (i in 1:num.of.file) {
      writeLines(paste0(stri_pad(vcf.names[i], 
                                 width = max.num.of.char,
                                 side = "right"), "  ",
                        stri_pad(mutation.loads$total.variants[i], 
                                 width = 15, side = "right"), "  ",
                        tools::md5sum(files[i]), "  ",
                        stri_pad(mutation.loads$SBS[i], width = 8,
                                 side = "right"), "  ",
                        stri_pad(mutation.loads$DBS[i], width = 8,
                                 side = "right"), "  ",
                        stri_pad(mutation.loads$ID[i], width = 7,
                                 side = "right"), "  ",
                        stri_pad(mutation.loads$discarded.variants[i], 
                                 width = 23, side = "right")), 
                 run.info)
      
    }
    # Add a disclaimer about discarded variants in the analysis
    writeLines("", run.info)
    writeLines(paste0("* For details about the various types of discarded ", 
                      "variants, please refer to element discarded.variants ",
                      "in the return value."), run.info)
    
    # Add strand bias statistics for SBS12 plot
    if (!is.null(strand.bias.statistics)) {
      writeLines("", run.info)
      writeLines("--- Transcription strand bias statistics ---", run.info)
      list0 <- strand.bias.statistics
      num.of.sample <- length(names(list0))
      space.mat <- CalculateNumberOfSpace(list0)
      
      for (i in 1:num.of.sample) {
        transcribed.counts <- list0[[i]][, "transcribed"]
        untranscribed.counts <- list0[[i]][, "untranscribed"]
        q.values <- list0[[i]][, "q.values"]
        q.values.symbol <- lapply(q.values, FUN = AssignNumberOfAsterisks)
        q.values.sci <- formatC(q.values, format = "e", digits = 2)
        
        transcribed.info <- character(0)
        untranscribed.info <- character(0)
        header1 <- header2 <- character(0)
        mutation.class <- rownames(list0[[1]])
        
        for (j in 1:6) {
          header1 <- paste0(header1, stri_pad(mutation.class[j], 
                                     width = space.mat[j, "space.total"], 
                                     side = "both"), "|")
          
          header2 <- 
            paste0(header2, " ", 
                   stri_pad("counts", 
                            width = space.mat[j, "space.counts"], 
                            side = "right"), " ",
                   stri_pad("Q-value", 
                            width = space.mat[j, "space.q.value"], 
                            side = "right"), " ", "|")
          
          transcribed.info <- 
            paste0(transcribed.info, " ",
                   stri_pad(transcribed.counts[j], 
                            width = space.mat[j, "space.counts"], 
                            side = "right"), " ", 
                   stri_pad(q.values.sci[j], 
                            width = space.mat[j, "space.q.value"], 
                            side = "right"), " ", "|")
          
          untranscribed.info <- 
            paste0(untranscribed.info, " ",
                   stri_pad(untranscribed.counts[j], 
                            width = space.mat[j, "space.counts"], 
                            side = "right"), " ", 
                   stri_pad(ifelse(is.null(q.values.symbol[[j]]), 
                                   "", q.values.symbol[[j]]), 
                            width = space.mat[j, "space.q.value"], 
                            side = "right"), " ", "|")
        }
        
        # Add description lines of the information listed for strand bias statistics
        writeLines(paste0(stri_pad("", width = 13), " |", header1), run.info)
        writeLines(paste0(stri_pad("Strand", width = 13, side = "right"), " |",
                          header2, "Sample name"), run.info)
        
        # Write the transcription strand bias statistics
        writeLines(paste0(stri_pad("transcribed", width = 13, side = "right"),
                          " |", transcribed.info, names(list0)[i]), run.info)
        writeLines(paste0(stri_pad("untranscribed", width = 13, side = "right"),
                          " |", untranscribed.info, names(list0)[i]), run.info)
        
        writeLines("", run.info)
      }
      
      # Add a description about the symbol denoting p-value
      writeLines(
        paste0("Legend: *Q<0.05, **Q<0.01, ***Q<0.001 (Benjamini-Hochberg ",
               "false discovery rates based on two-tailed binomial tests)"), run.info)
      
      # Add a note about direction of strand bias
      writeLines(paste0("Direction of strand bias: Fewer mutations on ",
                        "transcribed strand indicates that DNA damage occurred on ",
                        "pyrimidines,"), run.info) 
      writeLines(paste0("                          Fewer mutations on ", 
                        "untranscribed strand indicates that DNA damage occurred on ", 
                        "purines."), run.info)
    }
    close(run.info)
  }

#' Get mutation loads information from Mutect VCF files.
#'
#' @param catalogs A list generated by calling function
#'   \code{\link{MutectVCFFilesToCatalog}} to Mutect VCF files.
#'
#' @return A list containing mutation loads information from Mutect VCF files:
#'   
#' \enumerate{
#'  \item \code{total.variants} Total number of mutations.
#'  
#'  \item \code{SBS} Number of single base substitutions.
#'
#'  \item \code{DBS} Number of double base substitutions.
#'
#'  \item \code{ID} Number of small insertions and deletions.
#'
#'  \item \code{discarded.variants} Number of other types of mutations which are
#'  excluded in the analysis in \code{\link{ICAMS}}.
#'
#' }
#'   
#' @keywords internal
GetMutationLoadsFromMutectVCFs <- function(catalogs) {
  catSBS96 <- catalogs$catSBS96
  catDBS78 <- catalogs$catDBS78
  catID <- catalogs$catID
  vcf.names <- colnames(catSBS96)
  num.of.SBS <- colSums(catSBS96)
  num.of.DBS <- colSums(catDBS78)
  num.of.ID <- colSums(catID)
  
  if (is.null(catalogs$discarded.variants)) {
    num.of.total.variants <- num.of.SBS + num.of.DBS + num.of.ID
    num.of.discarded.variants <- 
      stats::setNames(rep(0, length(vcf.names)), vcf.names)
  } else {
    discarded.variants <- catalogs$discarded.variants
    mat <- matrix(data = 0, nrow = length(vcf.names),
                  ncol = length(discarded.variants), 
                  dimnames = list(vcf.names, names(discarded.variants)))
    
    for (type in names(discarded.variants)) {
      tmp <- discarded.variants[[type]]
      num.of.variants <- sapply(tmp, FUN = nrow)
      mat[names(num.of.variants), type] <- num.of.variants
    }
    num.of.discarded.variants <- rowSums(mat)
    num.of.total.variants <- 
      num.of.SBS + num.of.DBS + num.of.ID + num.of.discarded.variants
  }
  return(list(total.variants = num.of.total.variants,
              SBS = num.of.SBS, DBS = num.of.DBS, ID = num.of.ID,
              discarded.variants = num.of.discarded.variants))
}

#' Get mutation loads information from Strelka SBS VCF files.
#'
#' @param split.vcfs A list generated by calling function
#'   \code{\link{StrelkaSBSVCFFilesToCatalog}} to Strelka SBS VCF files.
#'
#' @return A list containing mutation loads information from Strelka SBS VCF
#'   files:
#'   
#' \enumerate{
#'  \item \code{total.variants} Total number of mutations.
#'  
#'  \item \code{SBS} Number of single base substitutions.
#'
#'  \item \code{DBS} Number of double base substitutions.
#'
#'  \item \code{ID} Number of small insertions and deletions.
#'
#'  \item \code{discarded.variants} Number of other types of mutations which are
#'  excluded in the analysis in \code{\link{ICAMS}}.
#'
#' }
#'   
#' @keywords internal
GetMutationLoadsFromStrelkaSBSVCFs <- function(catalogs) {
  catSBS96 <- catalogs$catSBS96
  catDBS78 <- catalogs$catDBS78
  vcf.names <- colnames(catSBS96)
  num.of.SBS <- colSums(catSBS96)
  num.of.DBS <- colSums(catDBS78)
  num.of.ID <- stats::setNames(rep(0, length(vcf.names)), vcf.names)
  
  if (is.null(catalogs$discarded.variants)) {
    num.of.total.variants <- num.of.SBS + num.of.DBS + num.of.ID
    num.of.discarded.variants <- 
      stats::setNames(rep(0, length(vcf.names)), vcf.names)
  } else {
    discarded.variants <- catalogs$discarded.variants
    mat <- matrix(data = 0, nrow = length(vcf.names),
                  ncol = length(discarded.variants), 
                  dimnames = list(vcf.names, names(discarded.variants)))
    
    for (type in names(discarded.variants)) {
      tmp <- discarded.variants[[type]]
      num.of.variants <- sapply(tmp, FUN = nrow)
      mat[names(num.of.variants), type] <- num.of.variants
    }
    num.of.discarded.variants <- rowSums(mat)
    num.of.total.variants <- 
      num.of.SBS + num.of.DBS + num.of.ID + num.of.discarded.variants
  }
  return(list(total.variants = num.of.total.variants,
              SBS = num.of.SBS, DBS = num.of.DBS, ID = num.of.ID,
              discarded.variants = num.of.discarded.variants))
}

#' Get mutation loads information from Strelka ID VCF files.
#'
#' @param list.of.vcfs A list generated by calling function
#'   \code{\link{ReadStrelkaIDVCFs}} to Strelka ID VCF files.
#'
#' @return A list containing mutation loads information from Strelka ID VCF
#'   files:
#'   
#' \enumerate{
#'  \item \code{total.variants} Total number of mutations.
#'  
#'  \item \code{SBS} Number of single base substitutions.
#'
#'  \item \code{DBS} Number of double base substitutions.
#'
#'  \item \code{ID} Number of small insertions and deletions.
#'
#'  \item \code{excluded.variants} Number of other types of mutations which are
#'  excluded in the analysis in \code{\link{ICAMS}}.
#'
#' }
#'   
#' @keywords internal
GetMutationLoadsFromStrelkaIDVCFs <- function(catalogs) {
  catID <- catalogs$catalog
  vcf.names <- colnames(catID)
  num.of.ID <- colSums(catID)
  num.of.SBS <- stats::setNames(rep(0, length(vcf.names)), vcf.names)
  num.of.DBS <- stats::setNames(rep(0, length(vcf.names)), vcf.names)
  
  if (is.null(catalogs$discarded.variants)) {
    num.of.total.variants <- num.of.SBS + num.of.DBS + num.of.ID
    num.of.discarded.variants <- 
      stats::setNames(rep(0, length(vcf.names)), vcf.names)
  } else {
    discarded.variants <- catalogs$discarded.variants
    mat <- matrix(data = 0, nrow = length(vcf.names),
                  ncol = length(discarded.variants), 
                  dimnames = list(vcf.names, names(discarded.variants)))
    
    for (type in names(discarded.variants)) {
      tmp <- discarded.variants[[type]]
      num.of.variants <- sapply(tmp, FUN = nrow)
      mat[names(num.of.variants), type] <- num.of.variants
    }
    num.of.discarded.variants <- rowSums(mat)
    num.of.total.variants <- 
      num.of.SBS + num.of.DBS + num.of.ID + num.of.discarded.variants
  }
  return(list(total.variants = num.of.total.variants,
              SBS = num.of.SBS, DBS = num.of.DBS, ID = num.of.ID,
              discarded.variants = num.of.discarded.variants))
}
