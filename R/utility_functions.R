#' Collapse catalog functions
#'
#' Collapse a catalog matrix
#'
#' \code{Collapse192To96} Collapse a SNS 192 catalog matrix to a SNS 96 catalog matrix.
#'
#' \code{Collapse1536To96} Collapse a SNS 1536 catalog matrix to a SNS 96 catalog matrix.
#'
#' \code{Collapse144To78} Collapse a DNS 144 catalog matrix to a DNS 78 catalog matrix.
#' @param catalog A catalog matrix to be collapsed whose row names indicate the
#'   mutation types while its columns show the occurrences of each mutation
#'   type of different samples.
#' @return A canonical catalog matrix whose row names indicate the mutation
#'   types while its columns show the occurrences of each mutation type of
#'   different samples.
#' @name CollapseCatalog
NULL

#' Transform nucleotide spectra functions
#'
#' Transform count spectra from a particular organism region to an inferred
#' count spectra based on the target nucleotide abundance.
#'
#' @param catalog A matrix of mutation counts. Rownames indicate the mutation
#'   types. Each column contains the mutation counts for one sample.
#'
#' @param source.abundance An abundance matrix specified by the user, which
#'   can be created using functions \code{\link{CreateDinucAbundance}},
#'   \code{\link{CreateTrinucAbundance}}, \code{\link{CreateTetranucAbundance}},
#'   \code{\link{CreatePentanucAbundance}}.
#'   There are 6 types of predefined abundance matrix which are incorporated
#'   in this function ("GRCh37.genome", "GRCh37.exome", "GRCh38.genome",
#'   "GRCh38.exome", "GRCm38.genome", "GRCm38.exome").
#'   User can invoke a specific predefined abundance matrix by typing its name,
#'   e.g. source.abundance = "GRCh37.genome".
#'
#' @param target.abundance An abundance matrix specified by the user, which
#'   can be created using functions \code{\link{CreateDinucAbundance}},
#'   \code{\link{CreateTrinucAbundance}}, \code{\link{CreateTetranucAbundance}},
#'   \code{\link{CreatePentanucAbundance}}.
#'   There are 6 types of predefined abundance matrix which are incorporated
#'   in this function ("GRCh37.genome", "GRCh37.exome", "GRCh38.genome",
#'   "GRCh38.exome", "GRCm38.genome", "GRCm38.exome").
#'   User can invoke a specific predefined abundance matrix by typing its name,
#'   e.g. target.abundance = "GRCm38.genome".
#'
#' @return A matrix of inferred mutation counts. Rownames indicate the mutation
#'   types which are the same as those in \code{catalog}.
#'   Each column contains the inferred mutation counts for one sample based on
#'   \code{target.abundance}.
#'
#' @keywords internal
#'
#' @name TransformSpectra
NULL

#' @rdname CollapseCatalog
#' @export
Collapse192To96 <- function(catalog) {
  dt192 <- data.table(catalog)
  dt192$rn <- PyrTri(rownames(catalog))
  dt96 <- dt192[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat96 <- as.matrix(dt96[, -1])
  rownames(mat96) <- dt96$rn
  mat96 <- mat96[catalog.row.order$SNS96, , drop = FALSE]
}

#' @rdname CollapseCatalog
#' @export
Collapse1536To96 <- function(catalog) {
  dt <- data.table(catalog)
  rn <- rownames(catalog)

  # The next gsub replaces the string representing a
  # single-base mutation in pentanucleotide with the corresponding
  # sring for that mutation in a trinucleotide context.
  dt$rn <- gsub(".(...).(.)", "\\1\\2", rn, perl = TRUE)
  dt96 <- dt[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat96 <- as.matrix(dt96[, -1])
  rownames(mat96) <- dt96$rn
  mat96 <- mat96[catalog.row.order$SNS96, , drop = FALSE]
}

#' @rdname CollapseCatalog
#' @export
Collapse144To78 <- function(catalog) {
  dt144 <- data.table(catalog)
  ref <- substr(rownames(catalog), 1, 2)
  alt <- substr(rownames(catalog), 3, 4)
  dt144$rn <- CanonicalizeDNS(ref, alt)
  dt78 <- dt144[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat78 <- as.matrix(dt78[ , -1])
  rownames(mat78) <- dt78$rn
  mat78 <- mat78[catalog.row.order$DNS78, , drop = FALSE]
}

#' Handle abundance (opportunity) specications uniformly.
#'
#' @param abundance Either an abunance variable or string specifying an abundance.
#'
#' @param which.n The n for the n-mers, one of 2, 3, 4, 5 for 2-mers, 3-mers, etc.
NormalizeAbundanceArg <- function(abundance, which.n) {
  if (class(abundance) %in% c("matrix", "numeric")) {
    # stopifnot .... which.n matches the abundance vector....
    return (abundance)
  }
  if (!which.n %in% 2:5) {
    stop("Argument which.n must be in the set 2:5, got", which.n)
  }

  if (!abundance %in% c("GRCh37.genome", "GRCh37.exome",
                               "GRCh38.genome", "GRCh38.exome",
                               "GRCm38.genome", "GRCm38.exome")) {
    stop ('abundance must be either an abundance matrix created by yourself
          or one of
          ("GRCh37.genome", "GRCh37.exome", "GRCh38.genome", "GRCh38.exome",
          "GRCm38.genome", "GRCm38.exome"), got', abundance)
  }

  if (abundance == "GRCh37.genome") {
    if (which.n == 2) return(abundance.2bp.genome.GRCh37)
    if (which.n == 3) return(abundance.3bp.genome.GRCh37)
    if (which.n == 4) return(abundance.4bp.genome.GRCh37)
    if (which.n == 5) return(abundance.5bp.genome.GRCh37)
  }
  if (abundance == "GRCh37.exome"){
    if (which.n == 2) return(abundance.2bp.exome.GRCh37)
    if (which.n == 3) return(abundance.3bp.exome.GRCh37)
    if (which.n == 4) return(abundance.4bp.exome.GRCh37)
    if (which.n == 5) return(abundance.5bp.exome.GRCh37)
  }
  if (abundance == "GRCh38.genome") {
    if (which.n == 2) return(abundance.2bp.genome.GRCh38)
    if (which.n == 3) return(abundance.3bp.genome.GRCh38)
    if (which.n == 4) return(abundance.4bp.genome.GRCh38)
    if (which.n == 5) return(abundance.5bp.genome.GRCh38)
  }
  if (abundance == "GRCh38.exome"){
    if (which.n == 2) return(abundance.2bp.exome.GRCh38)
    if (which.n == 3) return(abundance.3bp.exome.GRCh38)
    if (which.n == 4) return(abundance.4bp.exome.GRCh38)
    if (which.n == 5) return(abundance.5bp.exome.GRCh38)
  }
  stop("Programming error: we should never get here")
}

# TransDinucSpectra <- function(catalog, source.abundance, target.abundance, which.n,
# type = c("signature", "density", "counts"))
#  {
# source.abundance <- NormalizeAbundanceArg(source.abundance, which.n)
# target.abundance <- NormalizeAbundanceArg(target.abundance, which.n)
# !Do some error checking, something like
# if (nrow(catalog) == 96 && which.n != 3) stop....
# if (nrow(catalog) == 192 && which.n != 3) stop...
#
# stopifnot(names(source.abudance) == names(target.abundance))
#
# out.catalog <- catalog
#
# factor <- target.abundance / source.abundance
# names(factor) <- rownames(target.abundance)
#
# # CAUTION: this function depends on how mutations are encoded in
# the row names!
# transform.n.mer <- function(source.n.mer) {
#  # First, get the rows with the given source.n.mer
#  rows <- grep(paste("^", source.n.mer, sep=''), rownames(out.catalog))
#  out.catalog[rows, ] <<- out.catalog[rows, ] * factor[source.n.mer]
# }
#
# lapply(rownames(out.op), transform.n.mer)
#
# out2 <- apply(out.catalog, MARGIN = 2, function (x) x / sum(x))
# # Each colmun in out2 sums to 1
#
# if (type = "signature") return(out2)
#
# # lazy way to get new matrix in same shape as out2
# out3's elements will be overwritten
# out3 <- out2
#
# # This going back to counts making sure that the total number
# # of counts is the same as in the input. I think this
# # is one of several(?) possible design choices. Alternatively could
# # be to keep the counts of each major mutation class (e.g. C>A, C>G, C>T...)
# # unchanged.
# for (i in 1:ncol(out2)) {
#  out3[ ,i] <- out2[ ,i] * sum(input.sig.mat[ , i])
#}
# return(out3)
#}

#' @rdname TransformSpectra
#' @export
TransDinucSpectra <- function(catalog, source.abundance, target.abundance) {
  source.abundance <- NormalizeAbundanceArg(source.abundance, 3)
  target.abundance <- NormalizeAbundanceArg(target.abundance, 3)

  stopifnot(nrow(catalog) == 78)
  stopifnot(all(rownames(catalog) %in% catalog.row.order$DNS78) == TRUE)
  n <- ncol(catalog)
  per.dinuc.freq <- matrix(0, nrow = 78, ncol = n)
  inferred.count <- matrix(0, nrow = 78, ncol = n)

  for (i in 1:78) {
    for (j in 1:n) {
      if (class(source.abundance) == "matrix") {
        per.dinuc.freq[i, j] <-
          catalog[i, j] /
          source.abundance[substr(rownames(catalog)[i], 1, 2), ]
      } else  if (source.abundance == "GRCh37.genome") {
        per.dinuc.freq[i, j] <-
          catalog[i, j] /
          abundance.2bp.genome.GRCh37[substr(rownames(catalog)[i], 1, 2), ]
      } else if (source.abundance == "GRCh37.exome") {
        per.dinuc.freq[i, j] <-
          catalog[i, j] /
          abundance.2bp.exome.GRCh37[substr(rownames(catalog)[i], 1, 2), ]
      } else if (source.abundance == "GRCh38.genome") {
        per.dinuc.freq[i, j] <-
          catalog[i, j] /
          abundance.2bp.genome.GRCh38[substr(rownames(catalog)[i], 1, 2), ]
      } else if (source.abundance == "GRCh38.exome") {
        per.dinuc.freq[i, j] <-
          catalog[i, j] /
          abundance.2bp.exome.GRCh38[substr(rownames(catalog)[i], 1, 2), ]
      } else if (source.abundance == "GRCm38.genome") {
        per.dinuc.freq[i, j] <-
          catalog[i, j] /
          abundance.2bp.genome.GRCm38[substr(rownames(catalog)[i], 1, 2), ]
      } else if (source.abundance == "GRCm38.exome") {
        per.dinuc.freq[i, j] <-
          catalog[i, j] /
          abundance.2bp.exome.GRCm38[substr(rownames(catalog)[i], 1, 2), ]
      }

      if (class(target.abundance) == "matrix") {
        inferred.count[i, j] <-
          per.dinuc.freq[i, j] *
          target.abundance[substr(rownames(catalog)[i], 1, 2), ]
      } else if (target.abundance == "GRCh37.genome") {
        inferred.count[i, j] <-
          per.dinuc.freq[i, j] *
          abundance.2bp.genome.GRCh37[substr(rownames(catalog)[i], 1, 2), ]
      } else if (target.abundance == "GRCh37.exome") {
        inferred.count[i, j] <-
          per.dinuc.freq[i, j] *
          abundance.2bp.exome.GRCh37[substr(rownames(catalog)[i], 1, 2), ]
      } else if (target.abundance == "GRCh38.genome") {
        inferred.count[i, j] <-
          per.dinuc.freq[i, j] *
          abundance.2bp.genome.GRCh38[substr(rownames(catalog)[i], 1, 2), ]
      } else if (target.abundance == "GRCh38.exome") {
        inferred.count[i, j] <-
          per.dinuc.freq[i, j] *
          abundance.2bp.exome.GRCh38[substr(rownames(catalog)[i], 1, 2), ]
      } else if (target.abundance == "GRCm38.genome") {
        inferred.count[i, j] <-
          per.dinuc.freq[i, j] *
          abundance.2bp.genome.GRCm38[substr(rownames(catalog)[i], 1, 2), ]
      } else if (target.abundance == "GRCm38.exome") {
        inferred.count[i, j] <-
          per.dinuc.freq[i, j] *
          abundance.2bp.exome.GRCm38[substr(rownames(catalog)[i], 1, 2), ]
      }
    }
  }

  mat <- round(inferred.count, 0)
  rownames(mat) <- rownames(catalog)
  colnames(mat) <- colnames(catalog)
  return(mat)
}

#' @rdname TransformSpectra
#' @export
TransTrinucSpectra <- function(catalog, source.abundance, target.abundance) {
  if (class(source.abundance) != "matrix" &&
      !source.abundance %in% c("GRCh37.genome", "GRCh37.exome",
                               "GRCh38.genome", "GRCh38.exome",
                               "GRCm38.genome", "GRCm38.exome")) {
    stop ('source.abundance must be either an abundance matrix created by yourself
          or a type from
          ("GRCh37.genome", "GRCh37.exome", "GRCh38.genome", "GRCh38.exome",
          "GRCm38.genome", "GRCm38.exome")')
  }

  if (class(target.abundance) != "matrix" &&
      !target.abundance %in% c("GRCh37.genome", "GRCh37.exome",
                               "GRCh38.genome", "GRCh38.exome",
                               "GRCm38.genome", "GRCm38.exome")) {
    stop ('target.abundance must be either an abundance matrix created by yourself
          or a type from
          ("GRCh37.genome", "GRCh37.exome", "GRCh38.genome", "GRCh38.exome",
          "GRCm38.genome", "GRCm38.exome")')
  }

  stopifnot(nrow(catalog) == 96)
  stopifnot(all(rownames(catalog) %in% catalog.row.order$SNS96) == TRUE)
  n <- ncol(catalog)
  per.trinuc.freq <- matrix(0, nrow = 96, ncol = n)
  inferred.count <- matrix(0, nrow = 96, ncol = n)

  for (i in 1:96) {
    for (j in 1:n) {
      if (class(source.abundance) == "matrix") {
        per.trinuc.freq[i, j] <-
          catalog[i, j] /
          source.abundance[substr(rownames(catalog)[i], 1, 3), ]
      } else  if (source.abundance == "GRCh37.genome") {
        per.trinuc.freq[i, j] <-
          catalog[i, j] /
          abundance.3bp.genome.GRCh37[substr(rownames(catalog)[i], 1, 3), ]
      } else if (source.abundance == "GRCh37.exome") {
        per.trinuc.freq[i, j] <-
          catalog[i, j] /
          abundance.3bp.exome.GRCh37[substr(rownames(catalog)[i], 1, 3), ]
      } else if (source.abundance == "GRCh38.genome") {
        per.trinuc.freq[i, j] <-
          catalog[i, j] /
          abundance.3bp.genome.GRCh38[substr(rownames(catalog)[i], 1, 3), ]
      } else if (source.abundance == "GRCh38.exome") {
        per.trinuc.freq[i, j] <-
          catalog[i, j] /
          abundance.3bp.exome.GRCh38[substr(rownames(catalog)[i], 1, 3), ]
      } else if (source.abundance == "GRCm38.genome") {
        per.trinuc.freq[i, j] <-
          catalog[i, j] /
          abundance.3bp.genome.GRCm38[substr(rownames(catalog)[i], 1, 3), ]
      } else if (source.abundance == "GRCm38.exome") {
        per.trinuc.freq[i, j] <-
          catalog[i, j] /
          abundance.3bp.exome.GRCm38[substr(rownames(catalog)[i], 1, 3), ]
      }

      if (class(target.abundance) == "matrix") {
        inferred.count[i, j] <-
          per.trinuc.freq[i, j] *
          target.abundance[substr(rownames(catalog)[i], 1, 3), ]
      } else if (target.abundance == "GRCh37.genome") {
        inferred.count[i, j] <-
          per.trinuc.freq[i, j] *
          abundance.3bp.genome.GRCh37[substr(rownames(catalog)[i], 1, 3), ]
      } else if (target.abundance == "GRCh37.exome") {
        inferred.count[i, j] <-
          per.trinuc.freq[i, j] *
          abundance.3bp.exome.GRCh37[substr(rownames(catalog)[i], 1, 3), ]
      } else if (target.abundance == "GRCh38.genome") {
        inferred.count[i, j] <-
          per.trinuc.freq[i, j] *
          abundance.3bp.genome.GRCh38[substr(rownames(catalog)[i], 1, 3), ]
      } else if (target.abundance == "GRCh38.exome") {
        inferred.count[i, j] <-
          per.trinuc.freq[i, j] *
          abundance.3bp.exome.GRCh38[substr(rownames(catalog)[i], 1, 3), ]
      } else if (target.abundance == "GRCm38.genome") {
        inferred.count[i, j] <-
          per.trinuc.freq[i, j] *
          abundance.3bp.genome.GRCm38[substr(rownames(catalog)[i], 1, 3), ]
      } else if (target.abundance == "GRCm38.exome") {
        inferred.count[i, j] <-
          per.trinuc.freq[i, j] *
          abundance.3bp.exome.GRCm38[substr(rownames(catalog)[i], 1, 3), ]
      }
    }
  }

  mat <- round(inferred.count, 0)
  rownames(mat) <- rownames(catalog)
  colnames(mat) <- colnames(catalog)
  return(mat)
}

#' @rdname TransformSpectra
#' @export
TransTetranucSpectra <- function(catalog, source.abundance, target.abundance) {
  if (class(source.abundance) != "matrix" &&
      !source.abundance %in% c("GRCh37.genome", "GRCh37.exome",
                               "GRCh38.genome", "GRCh38.exome",
                               "GRCm38.genome", "GRCm38.exome")) {
    stop ('source.abundance must be either an abundance matrix created by yourself
          or a type from
          ("GRCh37.genome", "GRCh37.exome", "GRCh38.genome", "GRCh38.exome",
          "GRCm38.genome", "GRCm38.exome")')
  }

  if (class(target.abundance) != "matrix" &&
      !target.abundance %in% c("GRCh37.genome", "GRCh37.exome",
                               "GRCh38.genome", "GRCh38.exome",
                               "GRCm38.genome", "GRCm38.exome")) {
    stop ('target.abundance must be either an abundance matrix created by yourself
          or a type from
          ("GRCh37.genome", "GRCh37.exome", "GRCh38.genome", "GRCh38.exome",
          "GRCm38.genome", "GRCm38.exome")')
  }

  stopifnot(nrow(catalog) == 136)
  stopifnot(all(rownames(catalog) %in% catalog.row.order$DNS136) == TRUE)
  n <- ncol(catalog)
  per.tetranuc.freq <- matrix(0, nrow = 136, ncol = n)
  inferred.count <- matrix(0, nrow = 136, ncol = n)

  for (i in 1:136) {
    for (j in 1:n) {
      if (class(source.abundance) == "matrix") {
        per.tetranuc.freq[i, j] <-
          catalog[i, j] /
          source.abundance[rownames(catalog)[i], ]
      } else  if (source.abundance == "GRCh37.genome") {
        per.tetranuc.freq[i, j] <-
          catalog[i, j] /
          abundance.4bp.genome.GRCh37[rownames(catalog)[i], ]
      } else if (source.abundance == "GRCh37.exome") {
        per.tetranuc.freq[i, j] <-
          catalog[i, j] /
          abundance.4bp.exome.GRCh37[rownames(catalog)[i], ]
      } else if (source.abundance == "GRCh38.genome") {
        per.tetranuc.freq[i, j] <-
          catalog[i, j] /
          abundance.4bp.genome.GRCh38[rownames(catalog)[i], ]
      } else if (source.abundance == "GRCh38.exome") {
        per.tetranuc.freq[i, j] <-
          catalog[i, j] /
          abundance.4bp.exome.GRCh38[rownames(catalog)[i], ]
      } else if (source.abundance == "GRCm38.genome") {
        per.tetranuc.freq[i, j] <-
          catalog[i, j] /
          abundance.4bp.genome.GRCm38[rownames(catalog)[i], ]
      } else if (source.abundance == "GRCm38.exome") {
        per.tetranuc.freq[i, j] <-
          catalog[i, j] /
          abundance.4bp.exome.GRCm38[rownames(catalog)[i], ]
      }

      if (class(target.abundance) == "matrix") {
        inferred.count[i, j] <-
          per.tetranuc.freq[i, j] *
          target.abundance[rownames(catalog)[i], ]
      } else if (target.abundance == "GRCh37.genome") {
        inferred.count[i, j] <-
          per.tetranuc.freq[i, j] *
          abundance.4bp.genome.GRCh37[rownames(catalog)[i], ]
      } else if (target.abundance == "GRCh37.exome") {
        inferred.count[i, j] <-
          per.tetranuc.freq[i, j] *
          abundance.4bp.exome.GRCh37[rownames(catalog)[i], ]
      } else if (target.abundance == "GRCh38.genome") {
        inferred.count[i, j] <-
          per.tetranuc.freq[i, j] *
          abundance.4bp.genome.GRCh38[rownames(catalog)[i], ]
      } else if (target.abundance == "GRCh38.exome") {
        inferred.count[i, j] <-
          per.tetranuc.freq[i, j] *
          abundance.4bp.exome.GRCh38[rownames(catalog)[i], ]
      } else if (target.abundance == "GRCm38.genome") {
        inferred.count[i, j] <-
          per.tetranuc.freq[i, j] *
          abundance.4bp.genome.GRCm38[rownames(catalog)[i], ]
      } else if (target.abundance == "GRCm38.exome") {
        inferred.count[i, j] <-
          per.tetranuc.freq[i, j] *
          abundance.4bp.exome.GRCm38[rownames(catalog)[i], ]
      }
    }
  }

  mat <- round(inferred.count, 0)
  rownames(mat) <- rownames(catalog)
  colnames(mat) <- colnames(catalog)
  return(mat)
}

#' @rdname TransformSpectra
#' @export
TransPentanucSpectra <- function(catalog, source.abundance, target.abundance) {
  if (class(source.abundance) != "matrix" &&
      !source.abundance %in% c("GRCh37.genome", "GRCh37.exome",
                               "GRCh38.genome", "GRCh38.exome",
                               "GRCm38.genome", "GRCm38.exome")) {
    stop ('source.abundance must be either an abundance matrix created by yourself
          or a type from
          ("GRCh37.genome", "GRCh37.exome", "GRCh38.genome", "GRCh38.exome",
          "GRCm38.genome", "GRCm38.exome")')
  }

  if (class(target.abundance) != "matrix" &&
      !target.abundance %in% c("GRCh37.genome", "GRCh37.exome",
                               "GRCh38.genome", "GRCh38.exome",
                               "GRCm38.genome", "GRCm38.exome")) {
    stop ('target.abundance must be either an abundance matrix created by yourself
          or a type from
          ("GRCh37.genome", "GRCh37.exome", "GRCh38.genome", "GRCh38.exome",
          "GRCm38.genome", "GRCm38.exome")')
  }

  stopifnot(nrow(catalog) == 1536)
  stopifnot(all(rownames(catalog) %in% catalog.row.order$SNS1536) == TRUE)
  n <- ncol(catalog)
  per.pentanuc.freq <- matrix(0, nrow = 1536, ncol = n)
  inferred.count <- matrix(0, nrow = 1536, ncol = n)

  for (i in 1:1536) {
    for (j in 1:n) {
      if (class(source.abundance) == "matrix") {
        per.pentanuc.freq[i, j] <-
          catalog[i, j] /
          source.abundance[substr(rownames(catalog)[i], 1, 5), ]
      } else  if (source.abundance == "GRCh37.genome") {
        per.pentanuc.freq[i, j] <-
          catalog[i, j] /
          abundance.5bp.genome.GRCh37[substr(rownames(catalog)[i], 1, 5), ]
      } else if (source.abundance == "GRCh37.exome") {
        per.pentanuc.freq[i, j] <-
          catalog[i, j] /
          abundance.5bp.exome.GRCh37[substr(rownames(catalog)[i], 1, 5), ]
      } else if (source.abundance == "GRCh38.genome") {
        per.pentanuc.freq[i, j] <-
          catalog[i, j] /
          abundance.5bp.genome.GRCh38[substr(rownames(catalog)[i], 1, 5), ]
      } else if (source.abundance == "GRCh38.exome") {
        per.pentanuc.freq[i, j] <-
          catalog[i, j] /
          abundance.5bp.exome.GRCh38[substr(rownames(catalog)[i], 1, 5), ]
      } else if (source.abundance == "GRCm38.genome") {
        per.pentanuc.freq[i, j] <-
          catalog[i, j] /
          abundance.5bp.genome.GRCm38[substr(rownames(catalog)[i], 1, 5), ]
      } else if (source.abundance == "GRCm38.exome") {
        per.pentanuc.freq[i, j] <-
          catalog[i, j] /
          abundance.5bp.exome.GRCm38[substr(rownames(catalog)[i], 1, 5), ]
      }

      if (class(target.abundance) == "matrix") {
        inferred.count[i, j] <-
          per.pentanuc.freq[i, j] *
          target.abundance[substr(rownames(catalog)[i], 1, 5), ]
      } else if (target.abundance == "GRCh37.genome") {
        inferred.count[i, j] <-
          per.pentanuc.freq[i, j] *
          abundance.5bp.genome.GRCh37[substr(rownames(catalog)[i], 1, 5), ]
      } else if (target.abundance == "GRCh37.exome") {
        inferred.count[i, j] <-
          per.pentanuc.freq[i, j] *
          abundance.5bp.exome.GRCh37[substr(rownames(catalog)[i], 1, 5), ]
      } else if (target.abundance == "GRCh38.genome") {
        inferred.count[i, j] <-
          per.pentanuc.freq[i, j] *
          abundance.5bp.genome.GRCh38[substr(rownames(catalog)[i], 1, 5), ]
      } else if (target.abundance == "GRCh38.exome") {
        inferred.count[i, j] <-
          per.pentanuc.freq[i, j] *
          abundance.5bp.exome.GRCh38[substr(rownames(catalog)[i], 1, 5), ]
      } else if (target.abundance == "GRCm38.genome") {
        inferred.count[i, j] <-
          per.pentanuc.freq[i, j] *
          abundance.5bp.genome.GRCm38[substr(rownames(catalog)[i], 1, 5), ]
      } else if (target.abundance == "GRCm38.exome") {
        inferred.count[i, j] <-
          per.pentanuc.freq[i, j] *
          abundance.5bp.exome.GRCm38[substr(rownames(catalog)[i], 1, 5), ]
      }
    }
  }

  mat <- round(inferred.count, 0)
  rownames(mat) <- rownames(catalog)
  colnames(mat) <- colnames(catalog)
  return(mat)
}

#' Standardize the Chromosome name annotations for a data frame
#'
#' @param df A data frame whose first column contains the Chromosome name
#'
#' @return A data frame whose Chromosome names are only in the form of 1:22, "X"
#'   and "Y".
#' @keywords internal
StandardChromName <- function(df) {
  # Is there any row in df whose Chromosome names start with "GL"?
  if (sum(grepl("^GL", df[[1]])) > 0) {
    df <- df[-grep("^GL", df[[1]]), ]
  }

  # Is there any row in df whose Chromosome names are "Hs37D5"?
  if (sum(grepl("^Hs", df[[1]])) > 0) {
    df <- df[-grep("^Hs", df[[1]]), ]
  }

  # Is there any row in df whose Chromosome names contain "M"?
  if (sum(grepl("M", df[[1]])) > 0) {
    df <- df[-grep("M", df[[1]]), ]
  }

  # Remove the "chr" character in the Chromosome's name
  df[, 1] <- sub(pattern = "chr", replacement = "", df[[1]])

  return(df)
}

#' Create a Transcript Range file from the raw GFF3 File
#'
#' @param path The name/path of the raw GFF3 File, or a complete URL.
#'
#' @return A data frame which contains chromosome name, start, end position,
#'   strand information and gene name. Only the following four gene types are
#'   kept to facilitate transcriptional strand bias analysis: protein_coding,
#'   retained_intron, processed_transcript and nonsense_mediated_decay.
#' @keywords internal
CreateTransRange <- function(path) {
  df <- read.csv(path, header = FALSE, fill = TRUE, nrows = 20)
  # Count the number of comment lines
  n <- sum(grepl("#", df[, 1]))

  # Read in the raw GFF3 File while skipping the comment lines
  dt <- fread(path, header = FALSE, sep = "\t", fill = TRUE, skip = n)

  dt1 <- dt[dt$V3 == "gene", ]

  # Select out the four gene types for transcriptional strand bias analysis
  idx <-
    grepl("protein_coding", dt1$V9) |
    grepl("retained_intron", dt1$V9) |
    grepl("processed_transcript", dt1$V9) |
    grepl("nonsense_mediated_decay", dt1$V9)
  dt2 <- dt1[idx, ]

  # Split the 9th column of dt2 according to separator ";" and get a list
  list <- stringr::str_split(dt2$V9, ";")

  # Extract the character string which contains gene name information
  names <- sapply(list, stringr::str_subset, "gene_name")

  # Remove the "gene_name" characters
  names <- sub(pattern = "gene_name.", replacement = "", names)

  # Remove the quotation marks
  names <- gsub(pattern = '\"', replacement = "", names)

  # Remove the whitespace
  dt2$V9 <- gsub(pattern = "\\s", replacement = "", names)

  return(StandardChromName(dt2[, c(1, 4, 5, 7, 9)]))
}

#' PyrTri
#'
#' @param mutstring TODO
#'
#' @return TODO
#' @keywords internal
PyrTri <- function(mutstring) {
  # TODO (Steve) document
  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  output <-
    ifelse(substr(mutstring, 2, 2) %in% c("A", "G"),
           paste0(revc(substr(mutstring, 1, 3)),
                  revc(substr(mutstring, 4, 4))),
           mutstring)
  return(output)
}

#' PyrPenta
#'
#' @param mutstring TODO
#'
#' @return TODO
#' @keywords internal
PyrPenta <- function(mutstring) {
  # TODO (Steve) document
  stopifnot(nchar(mutstring) == rep(6, length(mutstring)))
  output <-
    ifelse(substr(mutstring, 3, 3) %in% c("A", "G"),
           paste0(revc(substr(mutstring, 1, 5)),
                  revc(substr(mutstring, 6, 6))),
           mutstring)
  return(output)
}

#' Reverse complement every string in \code{string.vec}.
#'
#' @param string.vec a vector of type character.
#' @importFrom Biostrings reverseComplement DNAStringSet
#' @return A vector of type characters with the reverse complement of every
#'   string in \code{string.vec}.
#' @export
revc <- function(string.vec) {
  return(
    as.character(reverseComplement(DNAStringSet(string.vec)))
  )
}

#' @title Reverse complement strings that represent stranded SNSs
#'
#' @param mutstring A vector of 4-character strings representing
#' stranded SNSs in trinucleotide context,
#' for example "AATC" represents AAT > ACT mutations.
#'
#' @return Return the vector of
#' reverse complements of the first 3 characters
#' concatenated with the reverse complement of the
#' last character, e.g. "AATC" returns "ATTG".
#'
#' @keywords internal
RevcSNS96 <- function(mutstring) {
  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  context <- revc(substr(mutstring, 1, 3))
  target  <- revc(substr(mutstring, 4, 4))
  return(paste0(context, target))
}

#' @title Reverse complement strings that represent stranded DNSs
#'
#' @param mutstring A vector of 4-character strings representing
#' stranded DNSs, for example "AATC" represents AA > TC mutations.
#'
#' @return Return the vector of
#' reverse complements of the first 2 characters
#' concatenated with the reverse complement of the second
#' 2 characters, e.g. "AATC" returns "TTGA".
#'
#' @keywords internal
RevcDNS144 <- function(mutstring) {
  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  context <- revc(substr(mutstring, 1, 2))
  target  <- revc(substr(mutstring, 3, 4))
  return(paste0(context, target))
}

#' Read transcript ranges and strands from a gff3 format file.
#' Use this one for the new, cut down gff3 file (2018 11 24)
#'
#' @param path Path to the file with the transcript information with 1-based
#'   start end positions of genomic ranges.
#'
#' @return A data.table keyed by chrom, chromStart, and chromEnd.
#' @keywords internal
ReadTranscriptRanges <- function(path) {
  d <- utils::read.table(path)
  colnames(d) <- c("chrom", "chromStart", "chromEnd", "strand", "name")
  bed1 <- data.table(d)
  data.table::setkeyv(bed1, c("chrom", "chromStart", "chromEnd"))
  return(bed1)
}

#' Read transcript ranges and strands from a bed format file.
#'
#' This function is mostly for testing purpose, may be removed in the future.
#'
#' @param path Path to the file with the transcript information (in bed format).
#'
#' @return A data.table keyed by chrom, chromStart, and chromEnd.
#' @export
#' @keywords internal
ReadBedTranscriptRanges <- function(path) {
  names <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")
  bed <- utils::read.table(path, col.names = names, as.is = TRUE)

  # Delete duplicate entries in the BED file
  bed <- dplyr::distinct(bed, chrom, chromStart, chromEnd, strand, .keep_all = TRUE)

  # Bed file are 0 based start and 1 based end (an oversimplification).
  # We need to add 1L and not 1, otherwise the column turns to a double
  # we get a warning from data.table.
  bed$chromStart <- bed$chromStart + 1L

  bed1 <- data.table(bed)
  data.table::setkeyv(bed1, c("chrom", "chromStart", "chromEnd"))
  return(bed1)
}

#' Create trinucleotide abundance file
#'
#' @param path Path to the file with the nucleotide abundance information with 3
#'   base pairs.
#'
#' @return A matrix whose row names indicate 32 different types of 3 base pairs
#'   combinations while its column contains the occurrences of each type.
#' @keywords internal
CreateTrinucAbundance <- function(path) {
  dt <- fread(path)
  colnames(dt) <- c("3bp", "occurrences")
  dt$type <-
    ifelse(substr(dt[[1]], 2, 2) %in% c("A", "G"), revc(dt[[1]]), dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  mat <- as.matrix(dt1[, 2])
  rownames(mat) <- dt1[[1]]
  return(mat)
}

#' Create dinucleotide abundance file
#'
#' @param path Path to the file with the nucleotide abundance information with 4
#'   base pairs.
#' @import data.table
#' @return A matrix whose row names indicate 10 different types of 2 base pairs
#'   combinations while its column contains the occurrences of each type.
#' @keywords internal
CreateDinucAbundance <- function(path) {
  dt <- fread(path)
  colnames(dt) <- c("4bp", "occurrences")
  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  dt$type <-
    ifelse(substr(dt[[1]], 2, 3) %in% canonical.ref,
           substr(dt[[1]], 2, 3),
           revc(substr(dt[[1]], 2, 3)))
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  mat <- as.matrix(dt1[, 2])
  rownames(mat) <- dt1[[1]]
  return(mat)
}

#' Create tetranucleotide abundance file
#'
#' @param path Path to the file with the nucleotide abundance information with 4
#'   base pairs.
#' @import data.table
#' @return A matrix whose row names indicate 136 different types of 4 base pairs
#'   combinations while its column contains the occurrences of each type.
#' @keywords internal
CreateTetranucAbundance <- function(path) {
  dt <- fread(path)
  colnames(dt) <- c("4bp", "occurrences")
  dt$type <- CanonicalizeQUAD(dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  mat <- as.matrix(dt1[, 2])
  rownames(mat) <- dt1[[1]]
  return(mat)
}

#' Create pentanucleotide abundance file
#'
#' @param path Path to the file with the nucleotide abundance information
#'   with 5 base pairs.
#' @import data.table
#' @return A matrix whose row names indicate 512 different types of 5 base
#'   pairs combinations while its column contains the occurrences of each type.
#' @keywords internal
CreatePentanucAbundance <- function(path) {
  dt <- fread(path)
  colnames(dt) <- c("5bp", "occurrences")
  dt$type <-
    ifelse(substr(dt[[1]], 3, 3) %in% c("A", "G"), revc(dt[[1]]), dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  mat <- as.matrix(dt1[, 2])
  rownames(mat) <- dt1[[1]]
  return(mat)
}

#' Take strings representing a genome and return the \link[BSgenome]{BSgenome} object.
#'
#' @param genome Either a variable containing a \link[BSgenome]{BSgenome} object
#'   or a character string acting as a genome identifier.
#'
#' @return If \code{genome} is \link[BSgenome]{BSgenome} object, return it.
#' Otherwise return the \link[BSgenome]{BSgenome} object identified by the
#' string \code{genome}.
#'
#' @keywords internal
NormalizeGenomeArg <- function(genome) {
  if (class(genome) == "character") {
    if (genome %in% c("GRCh38", "hg38")) {
      genome <- BSgenome.Hsapiens.UCSC.hg38
    } else if (genome %in% c("GRCh37", "hg19")) {
      genome <- BSgenome.Hsapiens.1000genomes.hs37d5
    } else {
      stop("Unrecoginzed genome identifier:\n", genome,
           "\nNeed one of GRCh38, hg38, GRCh37, hg19")
    }
  }
  return(genome)
}
