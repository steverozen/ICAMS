#' Create a zip file which contains catalogs and plot PDFs from Strelka SBS VCF files
#'
#' Create 3 SBS catalogs (96, 192, 1536), 3 DBS catalogs (78, 136, 144) from the
#' Strelka SBS VCFs specified by \code{dir}, save the catalogs as CSV files,
#' plot them to PDF and generate a zip archive of all the output files.
#'
#' This function calls \code{\link{StrelkaSBSVCFFilesToCatalog}},
#' \code{\link{PlotCatalogToPdf}}, \code{\link{WriteCatalog}} and
#' \code{\link[zip]{zipr}}.
#'
#' @param dir Pathname of the directory which contains the Strelka SBS VCF
#'   files. Each Strelka SBS VCF \strong{must} have a file extension ".vcf"
#'   (case insensitive) and share the \strong{same} \code{ref.genome} and
#'   \code{region}.
#'   
#' @inheritParams MutectVCFFilesToZipFile
#'
#' @importFrom utils glob2rx
#' 
#' @importFrom zip zipr 
#'
#' @return  A list of 3 SBS catalogs (one each for 96, 192, and 1536) and 3 DBS
#'   catalogs (one each for 78, 136, and 144). If trans.ranges = NULL, SBS 192
#'   and DBS 144 catalog will not be generated and plotted. Each catalog has
#'   attributes added. See \code{\link{as.catalog}} for more details.
#'
#' @note SBS 192 and DBS 144 catalogs include only mutations in transcribed
#'   regions. 
#' 
#' @inheritSection MutectVCFFilesToZipFile Comments
#' 
#' @export
#' 
#' @examples 
#' dir <- c(system.file("extdata/Strelka-SBS-vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <- 
#'     StrelkaSBSVCFFilesToZipFile(dir, 
#'                                 zipfile = file.path(tempdir(), "test.zip"),
#'                                 ref.genome = "hg19", 
#'                                 trans.ranges = trans.ranges.GRCh37,
#'                                 region = "genome",
#'                                 base.filename = "StrelkaSBS")
#'   unlink(file.path(tempdir(), "test.zip"))}
StrelkaSBSVCFFilesToZipFile <- function(dir,
                                        zipfile, 
                                        ref.genome, 
                                        trans.ranges = NULL, 
                                        region = "unknown", 
                                        names.of.VCFs = NULL, 
                                        base.filename = "",
                                        updateProgress = NULL) {
  .StrelkaSBSVCFFilesToZipFile(dir, zipfile, ref.genome, trans.ranges, region,
                               names.of.VCFs, base.filename, updateProgress)
}

#' The argument updateProgress is to be used in ICAMS.shiny package.
#' @keywords internal
.StrelkaSBSVCFFilesToZipFile <- function(dir,
                                         zipfile, 
                                         ref.genome, 
                                         trans.ranges = NULL, 
                                         region = "unknown", 
                                         names.of.VCFs = NULL, 
                                         base.filename = "",
                                         updateProgress = NULL) {
  files <- list.files(path = dir, pattern = "\\.vcf$", 
                      full.names = TRUE, ignore.case = TRUE)
  
  catalogs <-
    .StrelkaSBSVCFFilesToCatalog(files, ref.genome, trans.ranges, 
                                 region, names.of.VCFs, updateProgress)
  
  if (base.filename != "") {
    output.file <- file.path(tempdir(), paste0(base.filename, "."))
  } else {
    output.file <- file.path(tempdir(), base.filename)
  }
  
  for (name in names(catalogs)) {
    WriteCatalog(catalogs[[name]],
                 file = paste0(output.file, name, ".csv"))
  }
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "wrote catalogs to CSV files")
  }
  
  for (name in names(catalogs)) {
    PlotCatalogToPdf(catalogs[[name]],
                     file = paste0(output.file, name, ".pdf"))
    if (name == "catSBS192") {
      PlotCatalogToPdf(catalogs[[name]],
                       file = paste0(output.file, "SBS12.pdf"),
                       plot.SBS12 = TRUE)
    }
  }
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "plotted catalogs to PDF files")
  }
  
  file.names <- list.files(path = tempdir(), pattern = glob2rx("*.csv|pdf"), 
                           full.names = TRUE)
  zip::zipr(zipfile = zipfile, files = file.names)
  unlink(file.names)
  invisible(catalogs)
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
#' @param dir Pathname of the directory which contains the Strelka ID VCF files.
#'   Each Strelka ID VCF \strong{must} have a file extension ".vcf" (case
#'   insensitive) and share the \strong{same} \code{ref.genome} and
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
#' @return  A list of two elements. 1st element is an S3 object containing an ID
#'   (small insertion and deletion) catalog with class "IndelCatalog". See
#'   \code{\link{as.catalog}} for more details. 2nd element is a list of further
#'   annotated VCFs.
#'
#' @note In ID (small insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation deletion
#'   repeat sizes range from 1 to 6+.
#' 
#' @export
#' 
#' @examples 
#' dir <- c(system.file("extdata/Strelka-ID-vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <- 
#'     StrelkaIDVCFFilesToZipFile(dir, 
#'                                zipfile = file.path(tempdir(), "test.zip"),
#'                                ref.genome = "hg19", 
#'                                region = "genome",
#'                                base.filename = "StrelkaID")
#'   unlink(file.path(tempdir(), "test.zip"))}
StrelkaIDVCFFilesToZipFile <- function(dir,
                                       zipfile, 
                                       ref.genome, 
                                       region = "unknown", 
                                       names.of.VCFs = NULL, 
                                       base.filename = "",
                                       updateProgress = NULL){
  .StrelkaIDVCFFilesToZipFile(dir, zipfile, ref.genome, region,
                              names.of.VCFs, base.filename, updateProgress)
}

#' The argument updateProgress is to be used in ICAMS.shiny package.
#' @keywords internal
.StrelkaIDVCFFilesToZipFile <- function(dir,
                                        zipfile, 
                                        ref.genome, 
                                        region = "unknown", 
                                        names.of.VCFs = NULL, 
                                        base.filename = "",
                                        updateProgress = NULL){
  files <- list.files(path = dir, pattern = "\\.vcf$", 
                      full.names = TRUE, ignore.case = TRUE)
  list <-
    .StrelkaIDVCFFilesToCatalog(files, ref.genome, region, names.of.VCFs,
                                updateProgress)
  
  if (base.filename != "") {
    output.file <- file.path(tempdir(), paste0(base.filename, "."))
  } else {
    output.file <- file.path(tempdir(), base.filename)
  }
  
  WriteCatalog(list$catalog, file = paste0(output.file, "catID", ".csv"))
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "wrote catalog to CSV file")
  }
  
  PlotCatalogToPdf(list$catalog, file = paste0(output.file, "catID", ".pdf"))
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "plotted catalog to PDF file")
  }
  
  file.names <- list.files(path = tempdir(), pattern = glob2rx("*.csv|pdf"), 
                           full.names = TRUE)
  zip::zipr(zipfile = zipfile, files = file.names)
  unlink(file.names)
  invisible(list)
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
#' @param dir Pathname of the directory which contains the Mutect VCF files.
#'   Each Mutect VCF \strong{must} have a file extension ".vcf" (case
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
#' @param updateProgress Optional. Currently only used in ICAMS.shiny package to
#'   update the progress indicator.
#'
#' @importFrom utils glob2rx 
#' 
#' @importFrom zip zipr 
#'
#' @return  A list of 3 SBS catalogs (one each for 96, 192, and 1536), 3 DBS
#'   catalogs (one each for 78, 136, and 144) and Indel catalog. If trans.ranges
#'   = NULL, SBS 192 and DBS 144 catalog will not be generated and plotted. Each
#'   catalog has attributes added. See \code{\link{as.catalog}} for more
#'   details.
#'
#' @note SBS 192 and DBS 144 catalogs include only mutations in transcribed
#'   regions. In ID (small insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation deletion
#'   repeat sizes range from 1 to 6+.
#' 
#' @section Comments:
#' To add or change attributes of the catalog, you can use function \code{\link[base]{attr}}. \cr
#' For example, \code{attr(catalog, "abundance") <- custom.abundance}.
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
MutectVCFFilesToZipFile <- function(dir,
                                    zipfile, 
                                    ref.genome, 
                                    trans.ranges = NULL, 
                                    region = "unknown", 
                                    names.of.VCFs = NULL, 
                                    tumor.col.names = NA,
                                    base.filename = "",
                                    updateProgress = NULL){
  .MutectVCFFilesToZipFile(dir, zipfile, ref.genome, trans.ranges, region, 
                           names.of.VCFs, tumor.col.names, base.filename, 
                           updateProgress)
}

#' The argument updateProgress is to be used in ICAMS.shiny package.
#' @keywords internal
.MutectVCFFilesToZipFile <- function(dir,
                                     zipfile, 
                                     ref.genome, 
                                     trans.ranges = NULL, 
                                     region = "unknown", 
                                     names.of.VCFs = NULL, 
                                     tumor.col.names = NA,
                                     base.filename = "",
                                     updateProgress = NULL){
  files <- list.files(path = dir, pattern = "\\.vcf$", 
                      full.names = TRUE, ignore.case = TRUE)
  catalogs <-
    .MutectVCFFilesToCatalog(files, ref.genome, trans.ranges, 
                             region, names.of.VCFs, tumor.col.names,
                             updateProgress)
  
  if (base.filename != "") {
    output.file <- file.path(tempdir(), paste0(base.filename, "."))
  } else {
    output.file <- file.path(tempdir(), base.filename)
  }
  
  for (name in names(catalogs)) {
    WriteCatalog(catalogs[[name]],
                 file = paste0(output.file, name, ".csv"))
  }
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "wrote catalogs to CSV files")
  }
  
  for (name in names(catalogs)) {
    PlotCatalogToPdf(catalogs[[name]],
                     file = paste0(output.file, name, ".pdf"))
    if (name == "catSBS192") {
      PlotCatalogToPdf(catalogs[[name]],
                       file = paste0(output.file, "SBS12.pdf"),
                       plot.SBS12 = TRUE)
    }
  }
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "plotted catalogs to PDF files")
  }
  
  file.names <- list.files(path = tempdir(), pattern = glob2rx("*.csv|pdf"), 
                           full.names = TRUE)
  zip::zipr(zipfile = zipfile, files = file.names)
  unlink(file.names)
  invisible(catalogs)
}

#' Create SBS and DBS catalogs from Strelka SBS VCF files.
#'
#' Create 3 SBS catalogs (96, 192, 1536) and 3 DBS catalogs (78, 136, 144)
#' from the Strelka SBS VCFs specified by \code{files}
#'
#' This function calls \code{\link{VCFsToSBSCatalogs}} and
#' \code{\link{VCFsToDBSCatalogs}}.
#'
#' @param files Character vector of file paths to the Strelka SBS VCF files.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#' 
#' @return  A list of 3 SBS catalogs (one each for 96, 192, and 1536) and 3 DBS
#'   catalogs (one each for 78, 136, and 144). If trans.ranges = NULL, SBS 192
#'   and DBS 144 catalog will not be generated. Each catalog has attributes
#'   added. See \code{\link{as.catalog}} for more details.
#'
#' @note SBS 192 and DBS 144 catalog only contains mutations in transcribed regions.
#' 
#' @inheritSection MutectVCFFilesToCatalog Comments
#' 
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-SBS-vcf",
#'                       "Strelka.SBS.GRCh37.vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <- StrelkaSBSVCFFilesToCatalog(file, ref.genome = "hg19",
#'                                           trans.ranges = trans.ranges.GRCh37,
#'                                           region = "genome")}
StrelkaSBSVCFFilesToCatalog <-
  function(files, ref.genome, trans.ranges = NULL, 
           region = "unknown", names.of.VCFs = NULL) {
    .StrelkaSBSVCFFilesToCatalog(files, ref.genome, trans.ranges, 
                                 region, names.of.VCFs)
  }

#' The argument updateProgress is to be used in ICAMS.shiny package.
#' @keywords internal
.StrelkaSBSVCFFilesToCatalog <-
  function(files, ref.genome, trans.ranges = NULL, 
           region = "unknown", names.of.VCFs = NULL, updateProgress = NULL) {
    split.vcfs <- 
      .ReadAndSplitStrelkaSBSVCFs(files, names.of.VCFs, updateProgress)
    return(c(.VCFsToSBSCatalogs(split.vcfs$SBS.vcfs, ref.genome, 
                                trans.ranges, region, updateProgress),
             .VCFsToDBSCatalogs(split.vcfs$DBS.vcfs, ref.genome, 
                                trans.ranges, region, updateProgress)))
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
#' @return A list of two elements. 1st element is an S3 object containing an ID
#'   (small insertion and deletion) catalog with class "IndelCatalog". See
#'   \code{\link{as.catalog}} for more details. 2nd element is a list of further
#'   annotated VCFs.
#'
#' @note In ID (small insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-ID-vcf",
#'                       "Strelka.ID.GRCh37.vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catID <- StrelkaIDVCFFilesToCatalog(file, ref.genome = "hg19", 
#'                                       region = "genome")}
StrelkaIDVCFFilesToCatalog <- 
  function(files, ref.genome, region = "unknown", names.of.VCFs = NULL) {
    .StrelkaIDVCFFilesToCatalog(files, ref.genome, region, names.of.VCFs)
  }

#' The argument updateProgress is to be used in ICAMS.shiny package.
#' @keywords internal
.StrelkaIDVCFFilesToCatalog <- 
  function(files, ref.genome, region = "unknown", names.of.VCFs = NULL,
           updateProgress = NULL) {
    vcfs <- .ReadStrelkaIDVCFs(files, names.of.VCFs, updateProgress)
    return(.VCFsToIDCatalogs(vcfs, ref.genome, region, updateProgress))
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
#' @param updateProgress Optional. Currently only used in ICAMS.shiny package to
#'   update the progress indicator.
#'
#' @return  A list of 3 SBS catalogs (one each for 96, 192, and 1536), 3 DBS
#'   catalogs (one each for 78, 136, and 144) and ID catalog. If trans.ranges =
#'   NULL, SBS 192 and DBS 144 catalog will not be generated. Each catalog has
#'   attributes added. See \code{\link{as.catalog}} for more details.
#'
#' @note SBS 192 and DBS 144 catalogs include only mutations in transcribed
#'   regions. In ID (small insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation deletion
#'   repeat sizes range from 1 to 6+.
#'   
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf Comments
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <- MutectVCFFilesToCatalog(file, ref.genome = "hg19", 
#'                                       trans.ranges = trans.ranges.GRCh37,
#'                                       region = "genome")}
MutectVCFFilesToCatalog <-
  function(files, ref.genome, trans.ranges = NULL, region = "unknown", 
           names.of.VCFs = NULL, tumor.col.names = NA, updateProgress = NULL) {
    .MutectVCFFilesToCatalog(files, ref.genome, trans.ranges, region,
                             names.of.VCFs, tumor.col.names, updateProgress)
  }

#' The argument updateProgress is to be used in ICAMS.shiny package.
#' @keywords internal
.MutectVCFFilesToCatalog <-
  function(files, ref.genome, trans.ranges = NULL, region = "unknown", 
           names.of.VCFs = NULL, tumor.col.names = NA, updateProgress = NULL) {
    split.vcfs <- 
      .ReadAndSplitMutectVCFs(files, names.of.VCFs, tumor.col.names,
                              updateProgress)
    return(c(.VCFsToSBSCatalogs(split.vcfs$SBS, ref.genome, 
                                trans.ranges, region, updateProgress),
             .VCFsToDBSCatalogs(split.vcfs$DBS, ref.genome, 
                                trans.ranges, region, updateProgress),
             list(catID = .VCFsToIDCatalogs(split.vcfs$ID, ref.genome, 
                                            region, updateProgress)[[1]])))
  }

#' Read and split Strelka SBS VCF files.
#'
#' @param files Character vector of file paths to the Strelka SBS VCF files.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#'   
#' @return A list of 3 in-memory objects as follows:
#' \enumerate{
#'    \item \code{SBS.vcfs} List of data.frames of pure SBS mutations -- no DBS or 3+BS mutations.
#'
#'    \item \code{DBS.vcfs} List of data.frames of pure DBS mutations -- no SBS or 3+BS mutations.
#'
#'    \item \code{ThreePlus} List of data.tables with the key CHROM, LOW.POS, HIGH.POS. containing
#'    rows that that in the input that did not represent SBSs or DBSs.
#'
#'    }
#'
#' @seealso \code{\link{StrelkaSBSVCFFilesToCatalog}}
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-SBS-vcf",
#'                       "Strelka.SBS.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs(file)
ReadAndSplitStrelkaSBSVCFs <- function(files, names.of.VCFs = NULL) {
  .ReadAndSplitStrelkaSBSVCFs(files, names.of.VCFs)
}

#' The argument updateProgress is to be used in ICAMS.shiny package.
#' @keywords internal
.ReadAndSplitStrelkaSBSVCFs <- function(files, names.of.VCFs = NULL, 
                                        updateProgress = NULL) {
  vcfs <- ReadStrelkaSBSVCFs(files, names.of.VCFs)
  split.vcfs <- SplitListOfStrelkaSBSVCFs(vcfs)
  if (is.function(updateProgress)) {
    updateProgress(value = 0.2, detail = "read and split VCFs")
  }
  return(split.vcfs)
}

#' Read Strelka ID (small insertion and deletion) VCF files.
#'
#' @param files Character vector of file paths to the Strelka ID VCF files.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#' 
#' @return A list of vcfs from \code{files}.
#'
#' @note In ID (small insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#'
#' @seealso \code{\link{StrelkaIDVCFFilesToCatalog}}
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-ID-vcf",
#'                       "Strelka.ID.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadStrelkaIDVCFs(file)
ReadStrelkaIDVCFs <- function(files, names.of.VCFs = NULL) {
  .ReadStrelkaIDVCFs(files, names.of.VCFs)
}

#' The argument updateProgress is to be used in ICAMS.shiny package.
#' @keywords internal
.ReadStrelkaIDVCFs <- function(files, names.of.VCFs = NULL, 
                               updateProgress = NULL) {
  vcfs <- lapply(files, FUN = ReadStrelkaIDVCF)
  if (is.null(names.of.VCFs)) {
    names(vcfs) <- tools::file_path_sans_ext(basename(files))
  } else {
    CheckNamesOfVCFs(files, names.of.VCFs)
    names(vcfs) <- names.of.VCFs
  }
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "read VCFs")
  }
  return(vcfs)
}

#' Read and split Mutect VCF files.
#'
#' @param files Character vector of file paths to the Mutect VCF files.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#'   
#' @return A list with 3 in-memory VCFs and two left-over
#' VCF-like data frames with rows that were not incorporated
#' into the first 3 VCFs, as follows:
#'
#' \enumerate{
#'
#'  \item \code{SBS} VCF with only single base substitutions.
#'
#'  \item \code{DBS} VCF with only doublet base substitutions
#'   as called by Mutect.
#'
#'  \item \code{ID} VCF with only small insertions and deletions.
#'
#'  \item \code{other.subs} VCF like data.frame with
#'  rows for coordinate substitutions involving
#'  3 or more nucleotides, e.g. ACT > TGA or AACT > GGTA.
#'
#'  \item \code{multiple.alternative.alleles} VCF like data.frame with
#'  rows for variants with multiple alternative alleles, for example
#'  ACT mutated to both AGT and ACT at the same position.
#'
#' }
#'
#' @seealso \code{\link{MutectVCFFilesToCatalog}}
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitMutectVCFs(file)
ReadAndSplitMutectVCFs <- 
  function(files, names.of.VCFs = NULL, tumor.col.names = NA) {
    .ReadAndSplitMutectVCFs(files, names.of.VCFs, tumor.col.names)
  }

#' The argument updateProgress is to be used in ICAMS.shiny package.
#' @keywords internal
.ReadAndSplitMutectVCFs <- 
  function(files, names.of.VCFs = NULL, tumor.col.names = NA, 
           updateProgress = NULL) {
    vcfs <- ReadMutectVCFs(files, names.of.VCFs, tumor.col.names)
    split.vcfs <- SplitListOfMutectVCFs(vcfs)
    
    if (is.function(updateProgress)) {
      updateProgress(value = 0.1, detail = "read and split VCFs")
    }
    
    return(split.vcfs)
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
#' @return A list of 3 SBS catalogs, one each for 96, 192, 1536: catSBS96
#'   catSBS192 catSBS1536. If trans.ranges = NULL, SBS 192 catalog will not be
#'   generated. Each catalog has attributes added. See \code{\link{as.catalog}}
#'   for more details.
#'
#' @note SBS 192 catalogs only contain mutations in transcribed regions.
#'
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf Comments
#' 
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.SBS.vcfs <- ReadAndSplitMutectVCFs(file)$SBS
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs.SBS <- VCFsToSBSCatalogs(list.of.SBS.vcfs, ref.genome = "hg19",
#'                                     trans.ranges = trans.ranges.GRCh37,
#'                                     region = "genome")}
VCFsToSBSCatalogs <- function(list.of.SBS.vcfs, ref.genome, 
                              trans.ranges = NULL, region = "unknown") {
  .VCFsToSBSCatalogs(list.of.SBS.vcfs, ref.genome, trans.ranges, region)
}

#' The argument updateProgress is to be used in ICAMS.shiny package.
#' @keywords internal
.VCFsToSBSCatalogs <- function(list.of.SBS.vcfs, ref.genome, 
                               trans.ranges = NULL, region = "unknown",
                               updateProgress = NULL) {
  ncol <- length(list.of.SBS.vcfs)
  
  catSBS96 <- empty.cats$catSBS96
  catSBS192 <- empty.cats$catSBS192
  catSBS1536 <- empty.cats$catSBS1536
  trans.ranges <- InferTransRanges(ref.genome, trans.ranges)
  
  for (i in 1:ncol) {
    SBS.vcf <- list.of.SBS.vcfs[[i]]
    
    annotated.SBS.vcf <- AnnotateSBSVCF(SBS.vcf, ref.genome, trans.ranges)
    
    SBS.cat <- CreateOneColSBSMatrix(annotated.SBS.vcf)
    catSBS96 <- cbind(catSBS96, SBS.cat$catSBS96)
    if (!is.null(trans.ranges)) {
      catSBS192 <- cbind(catSBS192, SBS.cat$catSBS192)
    }
    catSBS1536 <- cbind(catSBS1536, SBS.cat$catSBS1536)
  }
  
  colnames(catSBS96) <- names(list.of.SBS.vcfs)
  colnames(catSBS1536) <- names(list.of.SBS.vcfs)
  
  catSBS96 <-
    as.catalog(catSBS96, ref.genome = ref.genome,
               region = region, catalog.type = "counts")
  
  catSBS1536 <-
    as.catalog(catSBS1536, ref.genome = ref.genome,
               region = region, catalog.type = "counts",
               abundance = NULL)
  if (is.null(trans.ranges)) {
    return(list(catSBS96 = catSBS96, catSBS1536 = catSBS1536))
  }
  
  colnames(catSBS192) <- names(list.of.SBS.vcfs)
  in.transcript.region <- ifelse(region == "genome", "transcript", region)
  catSBS192 <-
    as.catalog(catSBS192, ref.genome = ref.genome,
               region = in.transcript.region, 
               catalog.type = "counts",
               abundance = NULL)
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.3, detail = "generated SBS catalogs")
  }
  
  return(list(catSBS96 = catSBS96, catSBS192 = catSBS192, 
              catSBS1536 = catSBS1536))
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
#' @return A list of 3 DBS catalogs, one each for 78, 144, 136: catDBS78
#'   catDBS144 catDBS136. If trans.ranges = NULL, DBS 144 catalog will not be
#'   generated. Each catalog has attributes added. See \code{\link{as.catalog}}
#'   for more details.
#'
#' @note DBS 144 catalog only contains mutations in transcribed regions.
#'
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf Comments
#' 
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.DBS.vcfs <- ReadAndSplitMutectVCFs(file)$DBS
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs.DBS <- VCFsToDBSCatalogs(list.of.DBS.vcfs, ref.genome = "hg19",
#'                                     trans.ranges = trans.ranges.GRCh37,
#'                                     region = "genome")}
VCFsToDBSCatalogs <- function(list.of.DBS.vcfs, ref.genome, 
                              trans.ranges = NULL, region = "unknown") {
  .VCFsToDBSCatalogs(list.of.DBS.vcfs, ref.genome, trans.ranges, region)
}

#' The argument updateProgress is to be used in ICAMS.shiny package.
#' @keywords internal
.VCFsToDBSCatalogs <- function(list.of.DBS.vcfs, ref.genome, 
                               trans.ranges = NULL, region = "unknown",
                               updateProgress = NULL) {
  ncol <- length(list.of.DBS.vcfs)
  
  catDBS78 <- empty.cats$catDBS78
  catDBS136 <- empty.cats$catDBS136
  catDBS144 <- empty.cats$catDBS144
  trans.ranges <- InferTransRanges(ref.genome, trans.ranges)
  
  for (i in 1 : ncol) {
    DBS.vcf <- list.of.DBS.vcfs[[i]]
    
    annotated.DBS.vcf <- AnnotateDBSVCF(DBS.vcf, ref.genome, trans.ranges)
    
    DBS.cat <- CreateOneColDBSMatrix(annotated.DBS.vcf)
    catDBS78 <- cbind(catDBS78, DBS.cat$catDBS78)
    catDBS136 <- cbind(catDBS136, DBS.cat$catDBS136)
    if (!is.null(trans.ranges)) {
      catDBS144 <- cbind(catDBS144, DBS.cat$catDBS144)
    }
  }
  
  colnames(catDBS78) <- names(list.of.DBS.vcfs)
  colnames(catDBS136) <- names(list.of.DBS.vcfs)
  
  catDBS78 <-
    as.catalog(catDBS78, ref.genome = ref.genome,
               region = region, catalog.type = "counts")
  
  catDBS136 <-
    as.catalog(catDBS136, 
               ref.genome = ref.genome,
               region = region,
               catalog.type = "counts",
               abundance = NULL)
  
  if (is.null(trans.ranges)) {
    return(list(catDBS78 = catDBS78, catDBS136 = catDBS136))
  }
  colnames(catDBS144) <- names(list.of.DBS.vcfs)
  in.transcript.region <- ifelse(region == "genome", "transcript", region)
  catDBS144 <-
    as.catalog(catDBS144, 
               ref.genome = ref.genome,
               region = in.transcript.region, 
               catalog.type = "counts",
               abundance = NULL)
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.3, detail = "generated DBS catalogs")
  }
  
  return(list(catDBS78 = catDBS78, catDBS136 = catDBS136, 
              catDBS144 = catDBS144))
}

#' Create ID (small insertion and deletion) catalog from ID VCFs
#'
#' @param list.of.vcfs List of in-memory VCFs. The list names will be
#' the sample ids in the output catalog.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region A character string acting as a region identifier, one of
#' "genome", "exome".
#'
#' @return A list of two elements. 1st element is an S3 object containing an ID
#'   (small insertion and deletion) catalog with class "IndelCatalog". See
#'   \code{\link{as.catalog}} for more details. 2nd element is a list of further
#'   annotated VCFs.
#'   
#' @note In ID (small insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#'   
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-ID-vcf/",
#'                       "Strelka.ID.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.ID.vcfs <- ReadStrelkaIDVCFs(file)
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5",
#'  quietly = TRUE)) {
#'   catID <- VCFsToIDCatalogs(list.of.ID.vcfs, ref.genome = "hg19",
#'                             region = "genome")}
VCFsToIDCatalogs <- function(list.of.vcfs, ref.genome, region = "unknown") {
  .VCFsToIDCatalogs(list.of.vcfs, ref.genome, region)
}

#' The argument updateProgress is to be used in ICAMS.shiny package.
#' @keywords internal
.VCFsToIDCatalogs <- function(list.of.vcfs, ref.genome, region = "unknown",
                              updateProgress = NULL) {
  ncol <- length(list.of.vcfs)
  
  # Create a 0-column matrix with the correct row labels.
  catID <- matrix(0, nrow = length(ICAMS::catalog.row.order$ID), ncol = 0)
  rownames(catID) <- ICAMS::catalog.row.order$ID
  out.list.of.vcfs <- list()
  
  for (i in 1:ncol) {
    ID <- list.of.vcfs[[i]]
    ID <- AnnotateIDVCF(ID, ref.genome = ref.genome)
    # Unlike the case for SBS and DBS, we do not
    # add transcript information.
    tmp <- CreateOneColIDMatrix(ID)
    one.ID.column <- tmp[[1]]
    out.list.of.vcfs <- c(out.list.of.vcfs, list(tmp[[2]]))
    rm(ID)
    catID <- cbind(catID, one.ID.column)
  }
  
  colnames(catID) <- names(list.of.vcfs)
  names(out.list.of.vcfs) <- names(list.of.vcfs)
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "generated ID catalogs")
  }
  
  return(list(catalog = 
                as.catalog(catID, ref.genome = ref.genome,
                           region = region, catalog.type = "counts"),
              annotated.vcfs = out.list.of.vcfs))
}