#' Create ID (small insertions and deletions) catalog from ID VCFs
#'
#' @param list.of.vcfs List of in-memory ID VCFs. The list names will be
#' the sample ids in the output catalog. `AnnotateIDVCF`` has not been
#' on these VCFs.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param num.of.cores The number of cores to use. Not available on Windows
#'   unless \code{num.of.cores = 1}.
#'
#' @param region A character string acting as a region identifier, one of
#' "genome", "exome".
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#'
#' @section Value:
#' A \strong{list} of elements:
#'   * \code{catalog}: The ID (small insertions and deletions) catalog with
#'   attributes added. See \code{\link{as.catalog}} for details.
#'
#'   * \code{discarded.variants}: \strong{Non-NULL only if} there are variants
#'   that were excluded from the analysis. See the added extra column
#'   \code{discarded.reason} for more details.
#'
#'   * \code{annotated.vcfs}:
#' \strong{Non-NULL only if} \code{return.annotated.vcfs} = TRUE. A list of
#' data frames which contain the original VCF's ID mutation rows with three
#' additional columns \code{seq.context.width}, \code{seq.context} and
#' \code{ID.class} added. The category assignment of each ID mutation in VCF can
#' be obtained from \code{ID.class} column.
#' @md
#'
#' @inheritSection VCFsToCatalogsAndPlotToPdf ID classification
#'
#' @section Note:
#'  In ID (small insertions and deletions) catalogs, deletion repeat sizes range
#'  from 0 to 5+, but for plotting and end-user documentation deletion repeat
#'  sizes range from 1 to 6+.
#'
#' @export
#'
#' @examples
#' file <- c(system.file("extdata/Strelka-ID-vcf/",
#'                       "Strelka.ID.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' list.of.ID.vcfs <- ReadAndSplitVCFs(file, variant.caller = "strelka")$ID
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5",
#'  quietly = TRUE)) {
#'   catID <- VCFsToIDCatalogs(list.of.ID.vcfs, ref.genome = "hg19",
#'                             region = "genome")}
VCFsToIDCatalogs <- function(
  list.of.vcfs,
  ref.genome,
  num.of.cores = 1,
  trans.ranges = NULL,
  region = "unknown",
  flag.mismatches = 0,
  return.annotated.vcfs = FALSE,
  suppress.discarded.variants.warnings = TRUE
) {
  ncol <- length(list.of.vcfs)

  # Create 0-column matrices with the correct row labels.
  catID <- matrix(0, nrow = length(ICAMS::catalog.row.order$ID), ncol = 0)
  rownames(catID) <- ICAMS::catalog.row.order$ID
  catID166 <-
    matrix(0, nrow = length(ICAMS::catalog.row.order$ID166), ncol = 0)
  rownames(catID166) <- ICAMS::catalog.row.order$ID166

  annotated.vcfs <- discarded.variants <- list()

  GetIDCatalogs <- function(i, list.of.vcfs) {
    ID <- list.of.vcfs[[i]]
    sample.id <- names(list.of.vcfs)[i]

    if (suppress.discarded.variants.warnings == TRUE) {
      list <-
        suppressWarnings(AnnotateIDVCF(
          ID.vcf = ID,
          ref.genome = ref.genome,
          flag.mismatches = flag.mismatches,
          name.of.VCF = sample.id
        ))
    } else {
      list <- AnnotateIDVCF(
        ID.vcf = ID,
        ref.genome = ref.genome,
        flag.mismatches = flag.mismatches,
        name.of.VCF = sample.id
      )
    }

    # Create an empty data frame for discarded variants
    df <- ID[0, ]

    if (!is.null(list$discarded.variants)) {
      df <- dplyr::bind_rows(df, list$discarded.variants)
    }

    if (suppress.discarded.variants.warnings == TRUE) {
      tmp <- suppressWarnings({
        CreateOneColIDMatrix(
          list$annotated.vcf,
          sample.id = sample.id,
          return.annotated.vcf = return.annotated.vcfs
        )
      })
    } else {
      tmp <- CreateOneColIDMatrix(
        list$annotated.vcf,
        sample.id = sample.id,
        return.annotated.vcf = return.annotated.vcfs
      )
    }
    one.ID.column <- tmp$catalog
    one.ID166.column <- tmp$catID166
    one.ID476.column <- tmp$catID476
    rm(ID)

    if (return.annotated.vcfs == TRUE) {
      annotated.vcfs <- c(annotated.vcfs, list(tmp$annotated.vcf))
      names(annotated.vcfs) <- sample.id
    }
    if (!is.null(tmp$discarded.variants)) {
      df <- dplyr::bind_rows(df, tmp$discarded.variants)
    }
    if (nrow(df) != 0) {
      discarded.variants <- c(discarded.variants, list(df))
      names(discarded.variants) <- sample.id
    }

    return(list(
      ID.column = one.ID.column,
      ID166.column = one.ID166.column,
      ID476.column = one.ID476.column,
      discarded.variants = discarded.variants,
      annotated.vcfs = annotated.vcfs
    ))
  }

  list0 <- parallel::mclapply(
    1:ncol,
    FUN = GetIDCatalogs,
    list.of.vcfs = list.of.vcfs,
    mc.cores = num.of.cores
  )

  ID.cat <- lapply(list0, FUN = "[[", "ID.column")
  ID.cat1 <- do.call("cbind", ID.cat)
  catID <- as.catalog(
    ID.cat1,
    ref.genome = ref.genome,
    region = region,
    catalog.type = "counts"
  )

  ID166.cat <- lapply(list0, FUN = "[[", "ID166.column")
  ID166.cat1 <- do.call("cbind", ID166.cat)
  catID166 <- as.catalog(
    ID166.cat1,
    ref.genome = ref.genome,
    region = region,
    catalog.type = "counts"
  )

  ID476.cat <- lapply(list0, FUN = "[[", "ID476.column")
  # ID476.column can be NULL if Koh_476 was not in the VCF
  ID476.cat <- Filter(Negate(is.null), ID476.cat)
  catID476 <- NULL
  if (length(ID476.cat) > 0) {
    catID476 <- do.call("cbind", ID476.cat)
  }

  discarded.variants1 <- lapply(list0, FUN = "[[", "discarded.variants")
  discarded.variants2 <- do.call("c", discarded.variants1)

  annotated.vcfs1 <- lapply(list0, FUN = "[[", "annotated.vcfs")
  annotated.vcfs2 <- do.call("c", annotated.vcfs1)

  CheckAndReturnIDCatalog(
    catID = catID,
    catID166 = catID166,
    catID476 = catID476,
    discarded.variants = discarded.variants2,
    annotated.vcfs = annotated.vcfs2
  )
}
