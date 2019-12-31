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
#'                                 zipfile = paste0(tempdir(), "/test.zip"),
#'                                 ref.genome = "hg19", 
#'                                 trans.ranges = trans.ranges.GRCh37,
#'                                 region = "genome",
#'                                 base.filename = "StrelkaSBS")
#'   unlink(paste0(tempdir(), "/test.zip"))}
StrelkaSBSVCFFilesToZipFile <- function(dir,
                                        zipfile, 
                                        ref.genome, 
                                        trans.ranges = NULL, 
                                        region = "unknown", 
                                        names.of.VCFs = NULL, 
                                        base.filename = "") {
  .StrelkaSBSVCFFilesToZipFile(dir, zipfile, ref.genome, trans.ranges,
                               region, names.of.VCFs, base.filename)
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
    output.file <- paste0(tempdir(), "\\", base.filename, ".")
  } else {
    output.file <- paste0(tempdir(), "\\", base.filename)
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
    updateProgress(value = 0.1, detail = "read and split VCFs")
  }
  return(split.vcfs)
}