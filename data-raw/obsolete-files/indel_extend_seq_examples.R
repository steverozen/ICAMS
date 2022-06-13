#' @examples
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' split.vcfs <- ReadAndSplitVCFs(file, variant.caller = "mutect")
#' ID.catalog <- VCFsToIDCatalogs(list.of.vcfs = split.vcfs$ID,
#'                                ref.genome = "hg19",
#'                                region = "genome",
#'                                return.annotated.vcfs = TRUE)
#' annotated.vcf <- ID.catalog$annotated.vcfs$Mutect.GRCh37.s1
#' extended.seq.contexts <-
#'   SymmetricalContextsFor1BPIndel(annotated.vcf = annotated.vcf,
#'                                  indel.class = "DEL:T:1:0")
#'
#'



#' @examples
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' split.vcfs <- ReadAndSplitVCFs(file, variant.caller = "mutect")
#' ID.catalog <- VCFsToIDCatalogs(list.of.vcfs = split.vcfs$ID,
#'                                ref.genome = "hg19",
#'                                region = "genome",
#'                                return.annotated.vcfs = TRUE)
#' annotated.vcf <- ID.catalog$annotated.vcfs$Mutect.GRCh37.s1
#' extended.seq.contexts <-
#'   SymmetricalContextsFor1BPIndel(annotated.vcf = annotated.vcf,
#'                                  indel.class = "DEL:T:1:0")
#' GeneratePlotPFMmatrix(sequences = extended.seq.contexts,
#'                       indel.class = "DEL:T:1:0",
#'                       plot.dir = file.path(tempdir(), "test.pdf"),
#'                       plot.title = "Deletion of 1T from 1T")
#'                       


                   
#' @examples
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' split.vcfs <- ReadAndSplitVCFs(file, variant.caller = "mutect")
#' ID.catalog <- VCFsToIDCatalogs(list.of.vcfs = split.vcfs$ID,
#'                                ref.genome = "hg19",
#'                                region = "genome",
#'                                return.annotated.vcfs = TRUE)
#' annotated.vcf <- ID.catalog$annotated.vcfs$Mutect.GRCh37.s1
#' extended.seq.contexts <-
#'   SymmetricalContextsFor1BPIndel(annotated.vcf = annotated.vcf,
#'                                  indel.class = "INS:T:1:4")
#' ggplot.object <- HaplotypePlot(sequences = extended.seq.contexts,
#'                                indel.class = "INS:T:1:4",
#'                                title = "Deletion of 1T from 4Ts")
#' plot(ggplot.object)  