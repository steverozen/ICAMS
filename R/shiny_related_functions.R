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
StrelkaSBSVCFFilesToZipFile <- function(dir,
                                        zipfile, 
                                        ref.genome, 
                                        trans.ranges = NULL, 
                                        region = "unknown", 
                                        names.of.VCFs = NULL, 
                                        base.filename = "") {
  files <- list.files(path = dir, pattern = "\\.vcf$", 
                      full.names = TRUE, ignore.case = TRUE)
  vcf.names <- basename(files)
  split.vcfs <- ReadAndSplitStrelkaSBSVCFs(files, names.of.VCFs)
  mutation.loads <- GetMutationLoadsFromStrelkaSBSVCFs(split.vcfs)
  SBS.catalogs <- VCFsToSBSCatalogs(split.vcfs$SBS.vcfs, ref.genome, 
                                    trans.ranges, region)
  DBS.catalogs <- VCFsToDBSCatalogs(split.vcfs$DBS.vcfs, ref.genome, 
                                    trans.ranges, region)
  catalogs <- c(SBS.catalogs, DBS.catalogs)
  
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
      PlotCatalogToPdf(catalogs[[name]],
                       file = paste0(output.file, "SBS12.pdf"),
                       plot.SBS12 = TRUE)
    }
  }
  
  zipfile.name <- basename(zipfile)
  AddRunInformation(files, vcf.names, zipfile.name, vcftype = "strelka.sbs",
                    ref.genome, region, mutation.loads)
  file.names <- list.files(path = tempdir(), pattern = glob2rx("*.csv|pdf|txt"), 
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
StrelkaIDVCFFilesToZipFile <- function(dir,
                                       zipfile, 
                                       ref.genome, 
                                       region = "unknown", 
                                       names.of.VCFs = NULL, 
                                       base.filename = "",
                                       flag.mismatches = 0) {
  files <- list.files(path = dir, pattern = "\\.vcf$", 
                      full.names = TRUE, ignore.case = TRUE)
  vcf.names <- basename(files)
  list.of.vcfs <- ReadStrelkaIDVCFs(files, names.of.VCFs)
  mutation.loads <- GetMutationLoadsFromStrelkaIDVCFs(list.of.vcfs)
  list <- VCFsToIDCatalogs(list.of.vcfs, ref.genome, region, flag.mismatches)
  
  output.file <- ifelse(base.filename == "",
                        paste0(tempdir(), .Platform$file.sep),
                        file.path(tempdir(), paste0(base.filename, ".")))
  
  WriteCatalog(list$catalog, 
               file = paste0(output.file, "catID.csv"))
  
  PlotCatalogToPdf(list$catalog, 
                   file = paste0(output.file, "catID.pdf"))
  
  zipfile.name <- basename(zipfile)
  AddRunInformation(files, vcf.names, zipfile.name, vcftype = "strelka.id",
                    ref.genome, region, mutation.loads)
  file.names <- list.files(path = tempdir(), pattern = glob2rx("*.csv|pdf|txt"), 
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
#' @param flag.mismatches Optional. If > 0, then if there are mismatches to
#'   references in the ID (insertion/deletion) VCF, generate messages showing
#'   the mismatched rows and continue. Otherwise \code{stop} if there are
#'   mismatched rows. See \code{\link{AnnotateIDVCF}} for more details.
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
                                    flag.mismatches = 0){
  files <- list.files(path = dir, pattern = "\\.vcf$", 
                      full.names = TRUE, ignore.case = TRUE)
  vcf.names <- basename(files)
  split.vcfs <- ReadAndSplitMutectVCFs(files, names.of.VCFs, tumor.col.names)
  mutation.loads <- GetMutationLoadsFromMutectVCFs(split.vcfs)
  SBS.catalogs <- VCFsToSBSCatalogs(split.vcfs$SBS, ref.genome, 
                                    trans.ranges, region)
  DBS.catalogs <- VCFsToDBSCatalogs(split.vcfs$DBS, ref.genome, 
                                    trans.ranges, region)
  ID.catalog <- VCFsToIDCatalogs(split.vcfs$ID, ref.genome, 
                                 region, flag.mismatches)[[1]]
  catalogs <- c(SBS.catalogs, DBS.catalogs, list(catID = ID.catalog))
  
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
      PlotCatalogToPdf(catalogs[[name]],
                       file = paste0(output.file, "SBS12.pdf"),
                       plot.SBS12 = TRUE)
    }
  }
  
  zipfile.name <- basename(zipfile)
  AddRunInformation(files, vcf.names, zipfile.name, vcftype = "mutect", 
                    ref.genome, region, mutation.loads)
  file.names <- list.files(path = tempdir(), pattern = glob2rx("*.csv|pdf|txt"), 
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
  function(files, ref.genome, trans.ranges = NULL, region = "unknown", 
           names.of.VCFs = NULL) {
    split.vcfs <- ReadAndSplitStrelkaSBSVCFs(files, names.of.VCFs)
    return(c(VCFsToSBSCatalogs(split.vcfs$SBS.vcfs, ref.genome, 
                               trans.ranges, region),
             VCFsToDBSCatalogs(split.vcfs$DBS.vcfs, ref.genome, 
                               trans.ranges, region)))
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
  function(files, ref.genome, region = "unknown", names.of.VCFs = NULL,
           flag.mismatches = 0) {
    vcfs <- ReadStrelkaIDVCFs(files, names.of.VCFs)
    return(VCFsToIDCatalogs(vcfs, ref.genome, region, flag.mismatches))
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
           names.of.VCFs = NULL, tumor.col.names = NA, flag.mismatches = 0) {
    split.vcfs <- 
      ReadAndSplitMutectVCFs(files, names.of.VCFs, tumor.col.names)
    
    return(c(VCFsToSBSCatalogs(split.vcfs$SBS, ref.genome, 
                               trans.ranges, region),
             VCFsToDBSCatalogs(split.vcfs$DBS, ref.genome, 
                               trans.ranges, region),
             list(catID = VCFsToIDCatalogs(split.vcfs$ID, ref.genome, 
                                           region, flag.mismatches)[[1]])))
  }

#' Read and split Strelka SBS VCF files.
#'
#' @param files Character vector of file paths to the Strelka SBS VCF files.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#'   
#' @return A list of 3 in-memory objects as follows:
#' 
#' \enumerate{
#' 
#'    \item \code{SBS.vcfs} List of data.frames of pure SBS mutations -- no DBS
#'    or 3+BS mutations.
#'    
#'    \item \code{DBS.vcfs} List of data.frames of pure DBS mutations -- no SBS
#'    or 3+BS mutations.
#'
#'    \item \code{ThreePlus} List of data.tables with the key CHROM, LOW.POS,
#'    HIGH.POS which contain rows in the input that did not represent SBSs or
#'    DBSs.
#'    
#'    \item \code{multiple.alt} Rows with multiple alternate alleles (removed from
#'    \code{SBS.vcfs} etc.)
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
  vcfs <- ReadStrelkaSBSVCFs(files, names.of.VCFs)
  split.vcfs <- SplitListOfStrelkaSBSVCFs(vcfs)
  return(split.vcfs)
}

#' Read Strelka ID (small insertion and deletion) VCF files.
#'
#' @param files Character vector of file paths to the Strelka ID VCF files.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#'
#' @return A list of data frames containing data lines of the VCF files.
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
  vcfs <- 
    lapply(files, FUN = ReadStrelkaIDVCF, name.of.VCF = names.of.VCFs)
  if (is.null(names.of.VCFs)) {
    names(vcfs) <- tools::file_path_sans_ext(basename(files))
  } else {
    # Check whether the number of VCFs match the number of names
    # in names.of.VCFs
    CheckNamesOfVCFs(files, names.of.VCFs)
    names(vcfs) <- names.of.VCFs
  }
  return(vcfs)
}

#' Read and split Mutect VCF files.
#'
#' @param files Character vector of file paths to the Mutect VCF files.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#'   
#' @return A list with 3 in-memory VCFs and two left-over VCF-like data frames
#'   with rows that were not incorporated into the first 3 VCFs, as follows:
#'   
#' \enumerate{
#'
#'  \item \code{SBS} VCF with only single base substitutions.
#'
#'  \item \code{DBS} VCF with only doublet base substitutions.
#'
#'  \item \code{ID} VCF with only small insertions and deletions.
#'
#'  \item \code{other.subs} VCF like data.frame with rows for coordinate
#'  substitutions involving 3 or more nucleotides (e.g. ACT > TGA or AACT >
#'  GGTA) and rows for complex indels.
#'
#'  \item \code{multiple.alt} VCF like data.frame with
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
    vcfs <- ReadMutectVCFs(files, names.of.VCFs, tumor.col.names)
    split.vcfs <- SplitListOfMutectVCFs(vcfs)
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
  ncol <- length(list.of.SBS.vcfs)
  
  catSBS96 <- empty.cats$catSBS96
  catSBS192 <- empty.cats$catSBS192
  catSBS1536 <- empty.cats$catSBS1536
  trans.ranges <- InferTransRanges(ref.genome, trans.ranges)
  
  for (i in 1:ncol) {
    SBS.vcf <- list.of.SBS.vcfs[[i]]
    sample.id <- names(list.of.SBS.vcfs)[i]
    
    annotated.SBS.vcf <- AnnotateSBSVCF(SBS.vcf, ref.genome, trans.ranges)
    
    SBS.cat <- CreateOneColSBSMatrix(annotated.SBS.vcf, sample.id)
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
#' @param flag.mismatches Optional. If > 0, then if there are mismatches to
#'   references in the ID (insertion/deletion) VCF, generate messages showing
#'   the mismatched rows and continue. Otherwise \code{stop} if there are
#'   mismatched rows. See \code{\link{AnnotateIDVCF}} for more details.
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
VCFsToIDCatalogs <- function(list.of.vcfs, ref.genome, region = "unknown",
                             flag.mismatches = 0) {
  ncol <- length(list.of.vcfs)
  
  # Create a 0-column matrix with the correct row labels.
  catID <- matrix(0, nrow = length(ICAMS::catalog.row.order$ID), ncol = 0)
  rownames(catID) <- ICAMS::catalog.row.order$ID
  out.list.of.vcfs <- list()
  
  for (i in 1:ncol) {
    ID <- list.of.vcfs[[i]]
    ID <- AnnotateIDVCF(ID, ref.genome = ref.genome,
                        flag.mismatches = flag.mismatches)
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
  
  return(list(catalog = 
                as.catalog(catID, ref.genome = ref.genome,
                           region = region, catalog.type = "counts"),
              annotated.vcfs = out.list.of.vcfs))
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
           region, mutation.loads) {
    
    run.info <- 
      file(description = file.path(tempdir(), "run-information.txt"), open = "w")
    
    # Add the header information
    header <- paste0("run-information.txt file for ", zipfile.name, 
                     " created on ", Sys.time())
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
    writeLines(paste0(stringi::stri_pad("Name", width = max.num.of.char,
                                        side = "right"), "  ", 
                      "# of data lines", "  ",
                      stringi::stri_pad("MD5", width = 32,
                                        side = "right"), "  ",
                      "# of SBS", "  ",
                      "# of DBS", "  ",
                      "# of ID", "  ",
                      "# of excluded variants", "  "),
               run.info)
    
    num.of.file <- length(files)
    
    for (i in 1:num.of.file) {
      writeLines(paste0(stringi::stri_pad(vcf.names[i], 
                                          width = max.num.of.char,
                                          side = "right"), "  ",
                        stringi::stri_pad(mutation.loads$total.variants[i], 
                                          width = 15, side = "right"), "  ",
                        tools::md5sum(files[i]), "  ",
                        stringi::stri_pad(mutation.loads$SBS[i], width = 8,
                                          side = "right"), "  ",
                        stringi::stri_pad(mutation.loads$DBS[i], width = 8,
                                          side = "right"), "  ",
                        stringi::stri_pad(mutation.loads$ID[i], width = 7,
                                          side = "right"), "  ",
                        stringi::stri_pad(mutation.loads$excluded.variants[i], 
                                          width = 22, side = "right")), 
                 run.info)
                                          
    }
    # Add a disclaimer about excluded variants in the analysis
    writeLines("", run.info)
    writeLines("Disclaimer:", run.info)
    writeLines(paste0("Triplet and above base substitutions, ", 
                      "complex indels and variants with multiple alternate ",
                      "alleles are currently excluded in the analysis."), run.info)
    close(run.info)
  }

#' Get mutation loads information from Mutect VCF files.
#'
#' @param split.vcfs A list generated by calling function
#'   \code{\link{ReadAndSplitMutectVCFs}} to Mutect VCF files.
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
#'  \item \code{excluded.variants} Number of other types of mutations which are
#'  excluded in the analysis in \code{\link{ICAMS}}.
#'
#' }
#'   
#' @keywords internal
GetMutationLoadsFromMutectVCFs <- function(split.vcfs) {
  num.of.SBS <- sapply(split.vcfs$SBS, FUN = nrow)
  num.of.DBS <- sapply(split.vcfs$DBS, FUN = nrow)
  num.of.ID <- sapply(split.vcfs$ID, FUN = nrow)
  num.of.other.subs <- sapply(split.vcfs$other.subs, FUN = nrow)
  num.of.multiple.alt <- sapply(split.vcfs$multiple.alt, FUN = nrow)
  
  num.of.excluded.variants <- num.of.other.subs + num.of.multiple.alt
  num.of.total.variants <- 
    num.of.SBS + num.of.DBS + num.of.ID + num.of.excluded.variants
  
  return(list(total.variants = num.of.total.variants,
              SBS = num.of.SBS, DBS = num.of.DBS, ID = num.of.ID,
              excluded.variants = num.of.excluded.variants))
}

#' Get mutation loads information from Strelka SBS VCF files.
#'
#' @param split.vcfs A list generated by calling function
#'   \code{\link{ReadAndSplitStrelkaSBSVCFs}} to Strelka SBS VCF files.
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
#'  \item \code{excluded.variants} Number of other types of mutations which are
#'  excluded in the analysis in \code{\link{ICAMS}}.
#'
#' }
#'   
#' @keywords internal
GetMutationLoadsFromStrelkaSBSVCFs <- function(split.vcfs) {
  num.of.SBS <- sapply(split.vcfs$SBS.vcfs, FUN = nrow)
  num.of.DBS <- sapply(split.vcfs$DBS.vcfs, FUN = nrow)
  num.of.ID <- rep(0L, length(num.of.SBS))
  names(num.of.ID) <- names(num.of.SBS)
  num.of.threeplus <- sapply(split.vcfs$ThreePlus, FUN = nrow)
  num.of.multiple.alt <- sapply(split.vcfs$multiple.alt, FUN = nrow)
  
  num.of.excluded.variants <- num.of.threeplus + num.of.multiple.alt
  num.of.total.variants <- num.of.SBS + num.of.DBS + num.of.excluded.variants
  
  return(list(total.variants = num.of.total.variants,
              SBS = num.of.SBS, DBS = num.of.DBS, ID = num.of.ID,
              excluded.variants = num.of.excluded.variants))
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
GetMutationLoadsFromStrelkaIDVCFs <- function(list.of.vcfs) {
  num.of.total.variants <- sapply(list.of.vcfs, FUN = nrow)
  CheckComplexIndel <- function(ID.df) {
    complex.indels <- 
      which(substr(ID.df$REF, 1, 1) != substr(ID.df$ALT, 1, 1))
    return(length(complex.indels))
  }
  num.of.excluded.variants <- sapply(list.of.vcfs, FUN = CheckComplexIndel)
  num.of.ID <- num.of.total.variants - num.of.excluded.variants
  num.of.SBS <- rep(0L, length(num.of.ID))
  names(num.of.SBS) <- names(num.of.ID)
  num.of.DBS <- num.of.SBS
  
  return(list(total.variants = num.of.total.variants,
              SBS = num.of.SBS, DBS = num.of.DBS, ID = num.of.ID,
              excluded.variants = num.of.excluded.variants))
}