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
    updateProgress(value = 0.2, detail = "read VCFs")
  }
  return(vcfs)
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
    updateProgress(value = 0.4, detail = "generated SBS catalogs")
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
    updateProgress(value = 0.6, detail = "generated ID catalogs")
  }
  
  return(list(catalog = 
                as.catalog(catID, ref.genome = ref.genome,
                           region = region, catalog.type = "counts"),
              annotated.vcfs = out.list.of.vcfs))
}