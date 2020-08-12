#' "Collapse" a catalog
#' 
#' @description 
#' \enumerate{
#' \item Take a mutational spectrum or signature catalog
#' that is based on a fined-grained set of features (for example, single-nucleotide
#' substitutions in the context of the preceding and following 2 bases).
#'
#' \item Collapse it to a catalog based on a coarser-grained set of features
#' (for example, single-nucleotide substitutions in the context of the
#' immediately preceding and following bases).
#' }
#'
#' \code{Collapse192CatalogTo96} Collapse an SBS 192 catalog
#' to an SBS 96 catalog.
#'
#' \code{Collapse1536CatalogTo96} Collapse an SBS 1536 catalog
#'  to an SBS 96 catalog.
#'
#' \code{Collapse144CatalogTo78} Collapse a DBS 144 catalog
#' to a DBS 78 catalog.
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}}.
#'
#' @return A catalog as defined in \code{\link{ICAMS}}.
#'
#' @name CollapseCatalog
#' 
#' @examples 
#' # Create an SBS192 catalog and collapse it to an SBS96 catalog
#' object <- matrix(1, nrow = 192, ncol = 1, 
#'                  dimnames = list(catalog.row.order$SBS192))
#' catSBS192 <- as.catalog(object, region = "transcript")
#' catSBS96 <- Collapse192CatalogTo96(catSBS192)
NULL

#' @rdname CollapseCatalog
#' @export
Collapse192CatalogTo96 <- function(catalog) {
  dt192 <- data.table(catalog)
  dt192$rn <- PyrTri(rownames(catalog))
  dt96 <- dt192[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat96 <- as.matrix(dt96[, -1])
  rownames(mat96) <- dt96$rn
  mat96 <- mat96[ICAMS::catalog.row.order$SBS96, , drop = FALSE]
  
  abundance <- attributes(catalog)$abundance
  if (!is.null(abundance)) {
      abundance <- Collapse192AbundanceTo96(abundance)
  }
  
  cat96 <-
    as.catalog(
      object = mat96,
      ref.genome = attributes(catalog)$ref.genome,
      region = attributes(catalog)$region,
      catalog.type = attributes(catalog)$catalog.type,
      abundance = abundance
    )
  return(cat96)
}

#' @keywords internal
Collapse192AbundanceTo96 <- function(abundance192) {
  PyrTri <- function(string) {
    stopifnot(nchar(string) == rep(3, length(string)))
    output <-
      ifelse(substr(string, 2, 2) %in% c("A", "G"), revc(string), string)
    return(output)
  }
  dt <- data.table(abundance192)
  rownames(dt) <- names(abundance192)
  dt$rn <- PyrTri(rownames(dt))
  dt1 <- dt[, lapply(.SD, sum), by = rn, .SDcols = ]
  abundance96 <- unlist(dt1[, -1])
  names(abundance96) <- dt1$rn
  return(abundance96)
}

#' @rdname CollapseCatalog
#' @export
Collapse1536CatalogTo96 <- function(catalog) {
  dt <- data.table(catalog)
  rn <- rownames(catalog)

  # The next gsub replaces the string representing a
  # single-base mutation in pentanucleotide with the corresponding
  # string for that mutation in a trinucleotide context.
  dt$rn <- gsub(".(...).(.)", "\\1\\2", rn, perl = TRUE)
  dt96 <- dt[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat96 <- as.matrix(dt96[, -1])
  rownames(mat96) <- dt96$rn
  mat96 <- mat96[ICAMS::catalog.row.order$SBS96, , drop = FALSE]
  
  abundance <- attributes(catalog)$abundance
  if (!is.null(abundance)) {
    abundance <- Collapse5bpAbundanceTo3bp(abundance)
  }
  
  cat96 <-
    as.catalog(
      object = mat96,
      ref.genome = attributes(catalog)$ref.genome,
      region = attributes(catalog)$region,
      catalog.type = attributes(catalog)$catalog.type,
      abundance = abundance
    )
  return(cat96)
}

#' @keywords internal
Collapse5bpAbundanceTo3bp <- function(abundance1536) {
  dt <- data.table(abundance1536)
  rownames(dt) <- names(abundance1536)
  dt$rn <- substr(rownames(dt), 2, 4)
  dt1 <- dt[, lapply(.SD, sum), by = rn, .SDcols = ]
  abundance96 <- unlist(dt1[, -1])
  names(abundance96) <- dt1$rn
  return(abundance96)
}

#' @rdname CollapseCatalog
#' @export
Collapse144CatalogTo78 <- function(catalog) {
  dt144 <- data.table(catalog)
  ref <- substr(rownames(catalog), 1, 2)
  alt <- substr(rownames(catalog), 3, 4)
  dt144$rn <- CanonicalizeDBS(ref, alt)
  dt78 <- dt144[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat78 <- as.matrix(dt78[ , -1])
  rownames(mat78) <- dt78$rn
  mat78 <- mat78[ICAMS::catalog.row.order$DBS78, , drop = FALSE]
  
  abundance <- attributes(catalog)$abundance
  if (!is.null(abundance)) {
    abundance <- Collapse144AbundanceTo78(abundance)
  }
  
  cat78 <-
    as.catalog(
      object = mat78,
      ref.genome = attributes(catalog)$ref.genome,
      region = attributes(catalog)$region,
      catalog.type = attributes(catalog)$catalog.type,
      abundance = abundance)
    
  return(cat78)
}

#' @keywords internal
Collapse144AbundanceTo78 <- function(abundance144) {
  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  dt <- data.table(abundance144)
  rownames(dt) <- names(abundance144)
  dt$rn <- ifelse(rownames(dt) %in% canonical.ref, rownames(dt), 
                  revc(rownames(dt)))
  dt1 <- dt[, lapply(.SD, sum), by = rn, .SDcols = ]
  abundance78 <- unlist(dt1[, -1])
  names(abundance78) <- dt1$rn
  return(abundance78)
}

#' @keywords internal
IsDensity <- function(x) {
  return(x %in% c("density", "density.signature"))
}

#' @keywords internal
IsCounts <- function(x) {
  return(x %in% c("counts", "counts.signature"))
}

#' @keywords internal
IsSignature <- function(x) {
  return(x %in% c("counts.signature", "density.signature"))
}

#' @keywords internal
CheckAndNormalizeTranCatArgs <-
  function(catalog, 
           target.ref.genome, 
           target.region, 
           target.catalog.type, 
           target.abundance) {
    
    StopIfNrowIllegal(catalog)
    
    s <- list(
      ref.genome    = attr(catalog, "ref.genome",   exact = TRUE),
      catalog.type = attr(catalog,  "catalog.type", exact = TRUE),
      abundance    = attr(catalog,  "abundance",    exact = TRUE),
      region       = attr(catalog,  "region",       exact = TRUE))
    
    if (is.null(target.ref.genome))   target.ref.genome   <- s$ref.genome
    if (is.null(target.catalog.type)) target.catalog.type <- s$catalog.type 
    StopIfCatalogTypeIllegal(target.catalog.type)
    if (is.null(target.region))       target.region       <-s$region
    StopIfRegionIllegal(target.region)

    inferred.abundance <- 
        InferAbundance(catalog,
                       target.ref.genome,
                       target.region,
                       target.catalog.type)
                   
    if (!is.null(target.abundance) &&
        !is.null(inferred.abundance)) {
      if (!all.equal(target.abundance,
                     inferred.abundance)) {
        stop("Caller supplied abundance is different from inferred abundance")
      }
      stopifnot(names(inferred.abundance) == names(target.abundance))
    }
    if (is.null(target.abundance)) target.abundance <- inferred.abundance
    stopifnot(names(s[["abundance"]]) == names(target.abundance))

    if (s$catalog.type == target.catalog.type) {
      if (all(s[["abundance"]] == target.abundance)) { 
        warning("Input and target catalog.type and abundance are the same\n",
                "no transformation")
      }
    }
    return(list(s = s, 
                t = 
                  list(ref.genome   = target.ref.genome,
                       catalog.type = target.catalog.type,
                       region       = target.region,
                       abundance    = target.abundance)))
  }


#' @keywords internal
IsTransformationLegal <- function(s, t) {
 
  if (IsSignature(s$catalog.type) &&
      !IsSignature(t$catalog.type)) {
    stop("cannot transform from a signature to a counts or density catalog")
  }
  
  if (("COMPOSITECatalog" %in% class(s)) ||
       ("COMPOSITECatalog" %in% class(t))) {
    stop("Transformation of class COMPOSITECatalog not supported")
  }
  
  if (IsDensity(s$catalog.type))
    return(TCFromDenSigDen(s, t))
  else
    return(TCFromCouSigCou(s, t))
}


#' Source catalog type is counts or counts.signature
#' 
#' counts.signature -> density.signature, counts.signature
#' counts -> anything
#' 
#' @keywords internal
TCFromCouSigCou <- function(s, t) {
  
  stopifnot(IsCounts(s$catalog.type))
  
  if (IsSignature(s$catalog.type)) {
    # counts.signature -> density.signatyre, counts.signature
    return(TCFromCouSig(s, t))
  } else {
    # counts -> anything
    return(TCFromCou(s, t))
  }
}

#' @keywords internal
#' 
#' counts.signature -> counts.signature, density.signature
TCFromCouSig <- function(s, t) {
  
  if (t$catalog.type == "counts.signature") {
    # counts.signature -> counts.signature
    if (is.null(s[["abundance"]])) {
      stop("Cannot transform from counts.signature -> counts.signature ",
           "source abundance is NULL")
    }
    if (is.null(t[["abundance"]])) {
      stop("Cannot transform from counts.signature -> counts.signature ",
           "target abundance is NULL")
    }
    if (all(s[["abundance"]] == t[["abundance"]])) {
      warning("Transformation from counts.signature -> counts.sinature",
              "with equal abundance is a null operation")
      return("null.op")
    }
    return(TRUE)
  } else if (t$catalog.type == "counts") {
    # counts.signature -> counts
    stop("Cannot transform from counts.signature -> counts")
  } else if (t$catalog.type == "density.signature") {
    # counts.signature -> density.signature
     if (is.null(s[["abundance"]])) {
       stop("Cannot transform from counts.signature -> denisity.signature ",
            "if source abundance is NULL")
     } else {
       return(TRUE)
     }

  } else if (t$catalog.type == "density") {
    # counts.signature -> density
    stop("Cannnot tranform from counts.signature -> density")
  }
  stop("programming error, unexpected catalog.type ", t$catalog.type)
}


#' @keywords internal
#' 
#' counts -> <anything>
TCFromCou <- function(s, t) {
  if (t$catalog.type == "counts.signature") {
    # counts -> counts.signature
    if (is.null(s[["abundance"]])) {
      if (!is.null(t[["abundance"]]))
      stop("Cannot transform from counts -> counts.signature ",
           "source abundance is NULL and target abundance is not NULL")
    }
    if (is.null(t[["abundance"]])) {
      if (!is.null(s[["abundance"]]))
      stop("Cannot transform from counts -> counts.signature ",
           "target abundance is NULL and source abundance in not NULL")
    }
    # Source and target catalog both have NULL abundance
    if (is.null(s[["abundance"]]) && is.null(t[["abundance"]])) {
      # If source and target catalog have different ref.genome, 
      # raise an error
      if (!RefGenomeIsSame(s[["ref.genome"]], t[["ref.genome"]])) {
        stop("Cannot transform from counts -> counts.signature ",
             "source and target catalog both have NULL abundance,",
             "but have different ref.genome")
      }
      # If source and target catalog have different region, 
      # raise an error
      if (!RegionIsSame(s[["region"]], t[["region"]])) {
        stop("Cannot transform from counts -> counts.signature ",
             "source and target catalog both have NULL abundance,",
             "but have different region")
      }
    }
    if (AbundanceIsSame(s[["abundance"]], t[["abundance"]])) {
      return("sig.only")
    }
    return(TRUE)
  } else if (t$catalog.type == "counts") {
    # counts -> counts
    warning("counts -> counts is deprecated; ",
            "it simply infers new counts based on ",
            "changes in abundance; ", 
            "we strongly suggest that you work with densities")
    if (is.null(s[["abundance"]]) || is.null(t[["abundance"]])) {
      stop("Cannot transform from counts -> counts if either ",
           "source or target abundance is null")
    }
    if (all(s[["abundance"]] == t[["abundance"]])) {
      warning("Trasformation from counts -> counts with equal abundances ",
              "is a null operation")
      return("null.op")
    }
    return(TRUE)
        
  } else if (t$catalog.type == "density.signature") {
    # counts -> density.signature
    if (is.null(s[["abundance"]])) {
      stop("Cannot transform from counts -> denisity.signature ",
           "if source abundance is NULL")
    } 
    return(TRUE)
  } else if (t$catalog.type == "density") {
    # counts -> density
    if (is.null(s[["abundance"]])) {
      stop("Cannot transform from counts -> denisity ",
           "if source abundance is NULL")
    } 
    return(TRUE)
    
  }
  stop("programming error, unexpected catalog.type ", t$catalog.type)
}

#' density -> <anything>
#' density.signature -> density.signature, counts.signature
#' 
#' @keywords internal
TCFromDenSigDen <- function(s, t) {

  if (IsSignature(s$catalog.type))
    return(TCFromDenSig(s, t))
  else
    return(TCFromDen(s,t))
}

#' @keywords internal
#' 
#' density.signature -> anything
TCFromDenSig <- function (s, t) {
  if (t$catalog.type == "density") {
    # density.signature -> density
    stop("programming error: density.signature -> density")
  } else if (t$catalog.type == "density.signature") {
    # density.signature -> density.signature
    warning("Transform density.signature -> density signature",
            " is a null operation")
    return("null.op")
  } else if (t$catalog.type == "counts") {
    # density.signature -> counts
    stop("programming error: density.signature -> counts")
  } else if (t$catalog.type == "counts.signature") {
    # density.signature -> counts.signature
    if (is.null(t[["abundance"]])) {
      stop("Cannot transform density.signature -> counts.signature",
           " if target abundance is NULL")
    } else {  return(TRUE)   }
  }
  stop("programming error, unexpected catalog.type ", t$catalog.type)
}

#' @keywords internal
#' density -> anything
TCFromDen <- function (s, t) {  
  if (t$catalog.type == "density") {
    # density -> density
    warning("Transformation from density to density ",
            "is a null operation")
    return("null.op")
  } else if (t$catalog.type == "density.signature") { 
    # density -> density.signature
    return("sig.only")
  } else if (t$catalog.type == "counts") {
    # density -> counts
    if (is.null(t[["abundance"]])) {
      stop("Cannot tranform density -> counts; target abundance is NULL")
    } else {
      return(TRUE)
    }
  } else if (t$catalog.type == "counts.signature") {
    # density -> counts.signature
    if (is.null(t[["abundance"]])) {
      stop("Cannot tranform density -> counts.signature; target abundance is NULL")
    } else {
      return(TRUE)
    }    
  }
  stop("programming error, unexpected catalog.type ", t$catalog.type)
}

#' @keywords internal
AbundanceIsSame <- function(a1, a2) {
  if (is.null(a1)  && is.null(a2))  return(TRUE)
  if (is.null(a1)  && !is.null(a2)) return(FALSE)
  if (!is.null(a1) && is.null(a2))  return(FALSE)
  if (all(a1 == a2) && all(names(a1) == names(a2))) return(TRUE)
  return(FALSE)
}

#' @keywords internal
RefGenomeIsSame <- function(a1, a2) {
  if (is.null(a1)  && is.null(a2))  return(TRUE)
  if (is.null(a1)  && !is.null(a2)) return(FALSE)
  if (!is.null(a1) && is.null(a2))  return(FALSE)
  a1 <- NormalizeGenomeArg(a1)
  a2 <- NormalizeGenomeArg(a2)
  if (a1@pkgname == a2@pkgname) return(TRUE)
  return(FALSE)
}

#' @keywords internal
RegionIsSame <- function(a1, a2) {
  if (is.null(a1)  && is.null(a2))  return(TRUE)
  if (is.null(a1)  && !is.null(a2)) return(FALSE)
  if (!is.null(a1) && is.null(a2))  return(FALSE)
  if (a1 == a2) return(TRUE)
  return(FALSE)
}

#' @keywords internal
CheckCatalogAttributes <- function(catalog, target.catalog.type) {
  ref.genome <- attr(catalog, "ref.genome", exact = TRUE)
  catalog.type <- attr(catalog, "catalog.type", exact = TRUE)
  abundance <- attr(catalog, "abundance", exact = TRUE)
  region <- attr(catalog, "region", exact = TRUE)
  check.result <- TRUE
  
  if (is.null(ref.genome)) {
    check.result <<- TRUE
  } else if (inherits(ref.genome, "BSgenome")) {
    check.result <<- TRUE
  } else if (ref.genome == "mixed") {
    stop('Cannot perform transformation from a catalog with "mixed" ref.genome')
  }
  if (is.null(catalog.type)) {
    stop("Cannot perform transformation from a catalog with NULL catalog.type")
  }
  if (is.null(abundance)) {
    if (!(catalog.type == "counts" && target.catalog.type == "counts.signature")) {
      stop("Cannot perform transformation from a catalog with NULL abundance")
    }
  } else if (!inherits(abundance, "integer")) {
    stop("Cannot perform transformation from a catalog with non integer abundance")
  }
  if (is.null(region)) {
    stop("Cannot perform transformation from a catalog with NULL region")
  } else if (region == "mixed") {
    stop('Cannot perform transformation from a catalog with "mixed" region')
  }
  return(check.result)
}

#' Transform between counts and density spectrum catalogs
#' and counts and density signature catalogs
#'
#' @details Only the following transformations are legal:
#'
#' \enumerate{
#'
#' \item \code{counts -> counts} (deprecated, generates a warning;
#' we strongly suggest that you work with densities if comparing
#' spectra or signatures generated from data with
#' different underlying abundances.)
#' 
#' \item \code{counts -> density}
#'
#' \item \code{counts -> (counts.signature, density.signature)}
#'
#' \item \code{density -> counts} (the semantics are to
#' infer the genome-wide or exome-wide counts based on the
#' densities)
#'
#' \item \code{density -> density} (a null operation, generates
#' a warning)
#'
#' \item \code{density -> (counts.signature, density.signature)}
#'
#' \item \code{counts.signature -> counts.signature} (used to transform
#'    between the source abundance and \code{target.abundance})
#'
#' \item \code{counts.signature -> density.signature}
#' 
#' \item \code{counts.signature -> (counts, density)} (generates an error)
#' 
#' \item \code{density.signature -> density.signature} (a null operation,
#' generates a warning)
#'
#' \item \code{density.signature -> counts.signature}
#' 
#' \item \code{density.signature -> (counts, density)} (generates an error)
#'
#' }
#'
#'
#' @param catalog An SBS or DBS catalog as described in \code{\link{ICAMS}};
#'  must \strong{not} be an ID (small insertion and deletion) catalog.
#'
#' @param target.ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}. If \code{NULL}, then defaults to the
#'   \code{ref.genome} attribute of \code{catalog}.
#'
#' @param target.region A \code{region} argument; see \code{\link{as.catalog}}
#' and \code{\link{ICAMS}}. If \code{NULL}, then defaults to the
#'   \code{region} attribute of \code{catalog}.
#'
#' @param target.catalog.type A character string acting as a catalog type
#'   identifier, one of "counts", "density", "counts.signature",
#'   "density.signature"; see \code{\link{as.catalog}}. If \code{NULL}, then defaults to the
#'   \code{catalog.type} attribute of \code{catalog}.
#'   
#' @param target.abundance  
#'   A vector of counts, one for each source K-mer for mutations (e.g. for
#'   strand-agnostic single nucleotide substitutions in trinucleotide -- i.e.
#'   3-mer -- context, one count each for ACA, ACC, ACG, ... TTT). See
#'   \code{\link{all.abundance}}. If \code{NULL}, the function tries to infer
#'   \code{target.abundace} from the class of \code{catalog} and the value of
#'   the \code{target.ref.genome}, \code{target.region}, and
#'   \code{target.catalog.type}. If the \code{target.abundance} can be inferred
#'   and is different from a supplied non-\code{NULL} value of
#'   \code{target.abundance}, raise an error.
#'   
#' @return A catalog as defined in \code{\link{ICAMS}}.
#'
#' @export
#' 
#' @examples 
#' file <- system.file("extdata",
#'                     "strelka.regress.cat.sbs.96.csv",
#'                     package = "ICAMS")
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catSBS96.counts <- ReadCatalog(file, ref.genome = "hg19", 
#'                                  region = "genome",
#'                                  catalog.type = "counts")
#'   catSBS96.density <- TransformCatalog(catSBS96.counts,
#'                                        target.ref.genome = "hg19",
#'                                        target.region = "genome",
#'                                        target.catalog.type = "density")}
TransformCatalog <-
  function(catalog, 
           target.ref.genome   = NULL,
           target.region       = NULL, 
           target.catalog.type = NULL,
           target.abundance    = NULL) {
    # Check the attributes of the catalog
    stopifnot(CheckCatalogAttributes(catalog, target.catalog.type))
    
    # Check and normalize the arguments
    args <-
      CheckAndNormalizeTranCatArgs(
        catalog             = catalog,
        target.ref.genome   = target.ref.genome,
        target.catalog.type = target.catalog.type,
        target.region       = target.region,
        target.abundance    = target.abundance)
    rm(target.ref.genome, target.catalog.type, target.region, target.abundance)
    s <- args$s
    t <- args$t
    
    # Check if the transformation is legal
    test <- IsTransformationLegal(s, t)
    if (test == "null.op") return(catalog)
    if (test == "sig.only") {
      out <- apply(catalog, MARGIN = 2, function (x) x / sum(x))
      return(as.catalog(out, t$ref.genome,
                        t[["region"]], t[["catalog.type"]], t[["abundance"]])) 
    }

    factor <- t[["abundance"]] / s[["abundance"]]
    if (any(is.infinite(factor)) || any(is.na(factor))) {
      factor[is.infinite(factor)] <- 0
      factor[is.na(factor)] <- 0
    }

    names(factor) <- names(t[["abundance"]])
    out.catalog <- catalog

    for(source.n.mer in names(s[["abundance"]])) {
      # CAUTION: this loop depends on how mutations are encoded in
      # the row names of catalogs!
      
      # For 96 and 192 SBS, source.n.mer is e.g. "ACT" (for the encoding of ACT >
      # AGT as "ACTG"); for SBS1536 the n-mer for AACAG > AATAG is AACAG, in the
      # encoding AACAGT. For DBS78 and DBS144 TGGA represents TG >GA, and the
      # source n-mer is TG. For DBS136, TTGA represents TTGA > TNNA, and the
      # source n-mer is TTGA.
      
      # First, get the rows with the given source.n.mer
      rows <- grep(paste("^", source.n.mer, sep=''), rownames(out.catalog))
      # Then update those rows using the factor for that source.n.mer
      
      out.catalog[rows, ] <- out.catalog[rows, ] * factor[source.n.mer]
    }

    if (t[["catalog.type"]] %in% c("counts.signature", "density.signature")) {
      out2 <- apply(out.catalog, MARGIN = 2, function (x) x / sum(x))
    } else {
      out2 <- out.catalog
    }
    return(as.catalog(out2, t[["ref.genome"]], t[["region"]],
                      t[["catalog.type"]], t[["abundance"]]))

  }

#' Standardize the chromosome name annotations for a data frame.
#'
#' @param df A data frame whose first column contains the Chromosome name
#' 
#' @param file The name/path of the VCF file, or a complete URL.
#'
#' @return A \strong{list} whose first element is a data frame whose chromosome names are
#'   only in the form of 1:22, "X" and "Y". The second element
#'   \code{discarded.variants} \strong{only} appears if there are variants
#'   whose chromosome names belong to the following groups:
#'   1. Chromosome names that contain "GL".
#'   2. Chromosome names that contain "KI".
#'   3. Chromosome names that contain "random".
#'   4. Chromosome names that contain "Hs".
#'   5. Chromosome names that contain "M".
#' @md
#'   
#' @keywords internal
StandardChromName <- function(df, file) {
  # Create an empty data frame for discarded variants
  discarded.variants <- df[0, ]
  
  # Is there any row in df whose Chromosome names have "GL"?
  if (sum(grepl("GL", df[[1]])) > 0) {
    warning("In ", file, " ", sum(grepl("GL", df[[1]])), " row out of ",
            nrow(df), " had chromosome names that contain 'GL' and ", 
            "were removed. ",
            "See the discarded variants in the return value for more details")
    df1 <- df[-grep("GL", df[[1]]), ]
    discarded.variants <- 
      dplyr::bind_rows(discarded.variants, df[grep("GL", df[[1]]), ])
  } else {
    df1 <- df
  }
  
  # Is there any row in df whose Chromosome names have "KI"?
  if (sum(grepl("KI", df1[[1]])) > 0) {
    warning("In ", file, " ", sum(grepl("KI", df1[[1]])), " row out of ",
            nrow(df), " had chromosome names that contain 'KI' and ", 
            "were removed. ",
            "See the discarded variants in the return value for more details")
    df2 <- df1[-grep("KI", df1[[1]]), ]
    discarded.variants <- 
      dplyr::bind_rows(discarded.variants, df1[grep("KI", df1[[1]]), ])
  } else {
    df2 <- df1
  }
  
  # Is there any row in df whose Chromosome names have "random"?
  if (sum(grepl("random", df2[[1]])) > 0) {
    warning("In ", file, " ", sum(grepl("random", df2[[1]])), " row out of ",
            nrow(df), " had chromosome names that contain 'random' and ", 
            "were removed. ",
            "See the discarded variants in the return value for more details")
    df3 <- df2[-grep("random", df2[[1]]), ]
    discarded.variants <- 
      dplyr::bind_rows(discarded.variants, df2[grep("random", df2[[1]]), ])
  } else {
    df3 <- df2
  }
  
  # Is there any row in df whose Chromosome names are "Hs37D5"?
  if (sum(grepl("^Hs", df3[[1]])) > 0) {
    warning("In ", file, " ", sum(grepl("^Hs", df3[[1]])), " row out of ",
            nrow(df), " had chromosome names that contain 'Hs' and ", 
            "were removed. ",
            "See the discarded variants in the return value for more details")
    df4 <- df3[-grep("^Hs", df3[[1]]), ]
    discarded.variants <- 
      dplyr::bind_rows(discarded.variants, df3[-grep("^Hs", df3[[1]]), ])
  } else {
    df4 <- df3
  }
  
  # Is there any row in df whose Chromosome names contain "M"?
  if (sum(grepl("M", df4[[1]])) > 0) {
    warning("In ", file, " ", sum(grepl("M", df4[[1]])), " row out of ",
            nrow(df), " had chromosome names that contain 'M' and ", 
            "were removed. ",
            "See the discarded variants in the return value for more details")
    df5 <- df4[-grep("M", df4[[1]]), ]
    discarded.variants <- 
      dplyr::bind_rows(discarded.variants, df4[grep("M", df4[[1]]), ])
  } else {
    df5 <- df4
  }
  
  # Remove the "chr" character in the Chromosome's name
  df5[, 1] <- sub(pattern = "chr", replacement = "", df5[[1]])
  
  if (nrow(discarded.variants) == 0) {
    return(list(df = df5))
  } else {
    return(list(df = df5, discarded.variants = discarded.variants))
  }
}

#' Check and, if possible, correct the chromosome names in a VCF \code{data.frame}.
#' 
#' @param vcf.df A VCF as a \code{data.frame}. Check the names in column
#' \code{CHROM}.
#' 
#' @param name.of.VCF Name of the VCF file.
#' 
#' @param ref.genome The reference genome with the chromosome names to check
#' \code{vcf.df$CHROM} against; must be a Bioconductor 
#' \code{\link[BSgenome]{BSgenome}}, e.g.
#' \code{\link[BSgenome.Hsapiens.UCSC.hg38]{BSgenome.Hsapiens.UCSC.hg38}}.
#' 
#' @return If the \code{vcf.df$CHROM} values are correct or
#' can be corrected, then a vector of chromosome names
#' that can be used as a replacement for \code{vcf.df$CHROM}.
#' If the names in \code{vcf.df$CHROM} cannot be made to
#' be consistent with the chromosome names in \code{ref.genome},
#' then \code{stop}.
#' 
#' @keywords internal
CheckAndFixChrNames <- function(vcf.df, ref.genome, name.of.VCF = NULL) {
  names.to.check <- unique(vcf.df$CHROM)
  # Check whether the naming of chromosomes in vcf.df is consistent
  if(!sum(grepl("^chr", names.to.check)) %in% c(0, length(names.to.check))) {
    stop("\nNaming of chromosomes in VCF ", dQuote(name.of.VCF), 
         " is not consistent: ",
         paste(names.to.check, collapse = " "))
  }
  
  ref.genome.names <- seqnames(ref.genome)
  
  not.matched <- setdiff(names.to.check, ref.genome.names)
  
  # The names match -- we leave well-enough alone
  if (length(not.matched) == 0) return(vcf.df$CHROM)
  
  vcf.has.chr.prefix <- any(grepl(pattern = "^chr", names.to.check))
  ref.has.chr.prefix <- any(grepl(pattern = "^chr", ref.genome.names))
  
  new.chr.names <- vcf.df$CHROM
  if (ref.has.chr.prefix && !vcf.has.chr.prefix) {
    names.to.check <- paste0("chr", names.to.check)
    new.chr.names <- paste0("chr", new.chr.names)
    not.matched1 <- setdiff(names.to.check, ref.genome.names)
    if (length(not.matched1) == 0) return(new.chr.names)
  }
  
  if (!ref.has.chr.prefix && vcf.has.chr.prefix) {
    names.to.check <- gsub("chr", "", names.to.check)
    new.chr.names <- gsub("chr", "", names.to.check)
    not.matched2 <- setdiff(names.to.check, ref.genome.names)
    if (length(not.matched2) == 0) return(new.chr.names)
  }
  
  organism <- BSgenome::organism(ref.genome)
  
  CheckForPossibleMatchedChrName <- function(chr1, chr2) {
    if (chr1 %in% names.to.check) {
      # If chr2 is already in names.to.check, then stop
      if (chr2 %in% names.to.check) {
        stopmessage <- function() {
          stop("\n", chr1, " and ", chr2, " both are chromosome names in VCF ",
               dQuote(name.of.VCF),
               ", which should not be the case for ", organism, ". Please check ",
               "your data or specify the correct ref.genome argument")
        }
        if (vcf.has.chr.prefix) {
          if (grepl(pattern = "^chr", chr1)) {
            stopmessage()
          } else {
            chr1 <- paste0("chr", chr1)
            chr2 <- paste0("chr", chr2)
            stopmessage()
          }
        } else {
          if (!grepl(pattern = "^chr", chr1)) {
            stopmessage()
          } else {
            chr1 <- gsub("chr", "", chr1)
            chr2 <- gsub("chr", "", chr2)
            stopmessage()
          }
        }
      }
      
      new.chr.names[new.chr.names == chr1] <<- chr2
      names.to.check <- setdiff(names.to.check, chr1)
      names.to.check <<- unique(c(names.to.check, chr2))
    }
    }
      
  if (organism == "Homo sapiens") {
    
    # Maybe the problem is that X and Y are encoded as chr23 and chr24
    CheckForPossibleMatchedChrName("chr23", "chrX")
    CheckForPossibleMatchedChrName("chr24", "chrY")
    
    # Maybe the problem is that X and Y are encoded as 23 and 24
    CheckForPossibleMatchedChrName("23", "X")
    CheckForPossibleMatchedChrName("24", "Y")
  }
  
  if (organism == "Mus musculus") {
    
    # Maybe the problem is that X and Y are encoded as chr20 and chr21
    CheckForPossibleMatchedChrName("chr20", "chrX")
    CheckForPossibleMatchedChrName("chr21", "chrY")
    
    # Maybe the problem is that X and Y are encoded as 20 and 21
    CheckForPossibleMatchedChrName("20", "X")
    CheckForPossibleMatchedChrName("21", "Y")
  }
  
  not.matched3 <- setdiff(names.to.check, ref.genome.names)
  if (length(not.matched3) == 0) return(new.chr.names)
  
  stop("\nChromosome names in VCF ", dQuote(name.of.VCF), 
       " not in ref.genome for ", organism, ": ", 
       # We report the _original_ list of not matched names
       paste(not.matched, collapse = " "))
}

#' Create a transcript range file from the raw GFF3 File
#'
#' @param file The name/path of the raw GFF3 File, or a complete URL.
#'
#' @importFrom  stringi stri_split_fixed
#'
#' @return A data table which contains chromosome name, start, end position,
#'   strand information and gene name. It is keyed by chrom, start, and
#'   end. Only genes that are associated with a CCDS ID are kept for
#'   transcriptional strand bias analysis.
#'
#' @keywords internal
CreateTransRanges <- function(file) {
  df <- read.csv(file, header = FALSE, fill = TRUE, nrows = 20)
  # Count the number of comment lines
  n <- sum(grepl("#", df[, 1]))

  # Read in the raw GFF3 File while skipping the comment lines
  dt <- data.table::fread(file, header = FALSE, sep = "\t", fill = TRUE, skip = n)

  # Extract the gene ID associated with CCDS ID
  dt1 <- dt[grep("CCDS", dt$V9), ]
  info <- stri_split_fixed(dt1$V9, ";")
  gene.id.CCDS.idx <- lapply(info, grep, pattern = "gene_id")
  ExtractInfo <- function(idx, list1, list2) {
    return(list1[[idx]][list2[[idx]]])
  }
  gene.id.CCDS <-
    unique(sapply(1:length(info), ExtractInfo,
                  list1 = info, list2 = gene.id.CCDS.idx))
  gene.id.CCDS <- sapply(stri_split_fixed(gene.id.CCDS, "="), "[", 2)

  # Extract the gene ID and gene name from entries in dt which belong to the
  # type "gene"
  dt2 <- dt[dt$V3 == "gene", ]
  attributes.info <- stri_split_fixed(dt2$V9, ";")
  gene.id.idx <- lapply(attributes.info, grep, pattern = "gene_id")
  gene.name.idx <- lapply(attributes.info, grep, pattern = "gene_name")
  len <- length(attributes.info)
  gene.id <-
    sapply(1:len, ExtractInfo, list1 = attributes.info, list2 = gene.id.idx)
  gene.id <- sapply(stri_split_fixed(gene.id, "="), "[", 2)
  gene.name <-
    sapply(1:len, ExtractInfo, list1 = attributes.info, list2 = gene.name.idx)
  gene.name <- sapply(stri_split_fixed(gene.name, "="), "[", 2)
  dt3 <- dt2[, c("gene.id", "gene.name") := .(gene.id, gene.name)]

  # Only keep the entries in dt3 with gene ID that is associated with CCDS ID.
  # Select the necessary columns and standardize the chromosome names.
  dt4 <- dt3[gene.id %in% gene.id.CCDS, c(1, 4, 5, 7, 10, 11)]
  colnames(dt4) <- c("chrom", "start", "end", "strand", 
                     "Ensembl.gene.ID", "gene.name")
  dt5 <- StandardChromName(dt4)

  # Reorder dt5 according to chromosome name, start and end position
  chrOrder <- c((1:22), "X", "Y")
  dt5$chrom <- factor(dt5$chrom, chrOrder, ordered = TRUE)
  setkeyv(dt5, c("chrom", "start", "end"))

  # Check whether one gene.name have multiple entries in dt5
  dt6 <- dt5[, count := .N, by = .(chrom, gene.name)]
  dt7 <- dt6[count == 1, ]
  dt8 <- dt6[count != 1, ]
  if (nrow(dt8) == 0) {
    return(dt7[, c(1:6)])
  } else {
    # Check and remove entries in dt8 if it has readthrough transcripts
    dt9 <- dt[grep("readthrough", dt$V9), ]
    info1 <- stri_split_fixed(dt9$V9, ";")
    gene.id.readthrough.idx <- lapply(info1, grep, pattern = "gene_id")
    gene.id.readthrough <-
      unique(sapply(1:length(info1), ExtractInfo,
                    list1 = info1, list2 = gene.id.readthrough.idx))
    gene.id.readthrough1 <- 
      sapply(stri_split_fixed(gene.id.readthrough, "="), "[", 2)
    dt8$readthrough <- FALSE
    dt9 <- dt8[Ensembl.gene.ID %in% gene.id.readthrough1, readthrough := TRUE]
    dt10 <- dt9[readthrough == FALSE, ]
    dt11 <- dplyr::distinct(dt10, chrom, gene.name, .keep_all = TRUE)
    
    # Rbind dt7 and dt11
    dt12 <- rbind(dt7[, c(1:6)], dt11[, c(1:6)], use.names = FALSE)
    
    # Remove the string after dot in the Ensembl gene ID
    Ensembl.gene.ID.withoutdot <- gsub("\\..*", "", dt12$Ensembl.gene.ID)
    dt12$Ensembl.gene.ID <- Ensembl.gene.ID.withoutdot
    
    return(setkeyv(dt12, c("chrom", "start", "end")))
  }
}

#' Create exome transcriptionally stranded regions
#'
#' @param file Path to a SureSelect BED file which contains unstranded exome
#'   ranges.
#'
#' @param trans.ranges A data.table which contains transcript range and strand
#'   information. Please refer to \code{\link{TranscriptRanges}} for more
#'   details.
#'
#' @import data.table
#'
#' @return A data table which contains chromosome name, start, end position,
#'   strand information. It is keyed by chrom, start, and
#'   end.
#'
#' @keywords internal
CreateExomeStrandedRanges <- function(file, trans.ranges) {
  exome.ranges <- ReadBedRanges(file)
  colnames(exome.ranges) <- c("chrom", "exome.start", "exome.end")

  # Remove ranges which fall on transcripts on both strands and get
  # transcriptionally stranded regions
  trans.ranges <- RemoveRangesOnBothStrand(trans.ranges)

  # Find range overlaps between the exome.ranges and trans.ranges
  dt <- foverlaps(exome.ranges, trans.ranges, type = "any", mult = "all")
  dt1 <- dt[!is.na(strand)]

  # Find out and remove exome ranges that fall on transcripts on both strands
  dt2 <- dt1[, bothstrand := "+" %in% strand && "-" %in% strand,
             by = .(chrom, exome.start, exome.end)]
  dt3 <- dt2[bothstrand == FALSE]

  # One unstranded exome range can be represented by more than 1 row in dt3 if it
  # overlaps with multiple transcriptionally stranded regions. When creating the
  # exome transcriptionally stranded regions, we only need to count those ranges once.
  dt4 <- dt3[, count := .N, by = .(chrom, exome.start, exome.end)]
  dt4 <- dt3[, .(strand = strand[1]), by = .(chrom, exome.start, exome.end)]

  # Intersect dt4 ranges with the transcriptionally stranded regions
  gr1 <- GenomicRanges::makeGRangesFromDataFrame(dt4)
  gr2 <- GenomicRanges::makeGRangesFromDataFrame(trans.ranges)
  gr3 <- GenomicRanges::intersect(gr1, gr2)

  dt5 <- as.data.table(gr3)[, c(1:3, 5)]
  colnames(dt5)[1] <- "chrom"
  chrOrder <- c((1:22), "X", "Y")
  dt5$chrom <- factor(dt5$chrom, chrOrder, ordered = TRUE)
  return(setkeyv(dt5, c("chrom", "start", "end")))
}

#' @keywords internal
PyrTri <- function(mutstring) {
  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  output <-
    ifelse(substr(mutstring, 2, 2) %in% c("A", "G"),
           paste0(revc(substr(mutstring, 1, 3)),
                  revc(substr(mutstring, 4, 4))),
           mutstring)
  return(output)
}

#' @keywords internal
PyrPenta <- function(mutstring) {
  stopifnot(nchar(mutstring) == rep(6, length(mutstring)))
  output <-
    ifelse(substr(mutstring, 3, 3) %in% c("A", "G"),
           paste0(revc(substr(mutstring, 1, 5)),
                  revc(substr(mutstring, 6, 6))),
           mutstring)
  return(output)
}

#' Reverse complement every string in \code{string.vec}
#' 
#' Based on \code{\link{reverseComplement}}.
#' Handles IUPAC ambiguity codes but not "u" (uracil). \cr
#' (see <https://en.wikipedia.org/wiki/Nucleic_acid_notation>).
#'
#' @param string.vec A character vector.
#'
#' @importFrom Biostrings reverseComplement DNAStringSet
#'
#' @return A character vector with the reverse complement of every
#'   string in \code{string.vec}.
#'
#' @export
#' 
#' @examples 
#' revc("aTgc") # GCAT
#' 
#' # A vector and strings with ambiguity codes
#' revc(c("ATGC", "aTGc", "wnTCb")) # GCAT GCAT VGANW
#' 
#' \dontrun{
#' revc("ACGU") # An error}
revc <- function(string.vec) {
  return(
    as.character(reverseComplement(DNAStringSet(string.vec)))
  )
}

#' @title Reverse complement strings that represent stranded SBSs
#'
#' @param mutstring A vector of 4-character strings representing
#' stranded SBSs in trinucleotide context,
#' for example "AATC" represents AAT > ACT mutations.
#'
#' @return Return the vector of
#' reverse complements of the first 3 characters
#' concatenated with the reverse complement of the
#' last character, e.g. "AATC" returns "ATTG".
#'
#' @keywords internal
RevcSBS96 <- function(mutstring) {
  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  context <- revc(substr(mutstring, 1, 3))
  target  <- revc(substr(mutstring, 4, 4))
  return(paste0(context, target))
}

#' @title Reverse complement strings that represent stranded DBSs
#'
#' @param mutstring A vector of 4-character strings representing
#' stranded DBSs, for example "AATC" represents AA > TC mutations.
#'
#' @return Return the vector of
#' reverse complements of the first 2 characters
#' concatenated with the reverse complement of the second
#' 2 characters, e.g. "AATC" returns "TTGA".
#'
#' @keywords internal
RevcDBS144 <- function(mutstring) {
  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  context <- revc(substr(mutstring, 1, 2))
  target  <- revc(substr(mutstring, 3, 4))
  return(paste0(context, target))
}

#' Read transcript ranges and strand information from a gff3 format file.
#' Use this one for the new, cut down gff3 file (2018 11 24)
#'
#' @param file Path to the file with the transcript information with 1-based
#'   start end positions of genomic ranges.
#'
#' @return A data.table keyed by chrom, start, and end.
#'
#' @keywords internal
ReadTranscriptRanges <- function(file) {
  dt <- data.table::fread(file)
  colnames(dt) <- c("chrom", "start", "end", "strand", 
                    "Ensembl.gene.ID", "gene.symbol")
  
  # Check whether the transcript ranges come from human or mouse genomes
  if (length(unique(dt$chrom)) == 24) {
    chrOrder <- c((1:22), "X", "Y")
  } else if (length(unique(dt$chrom)) == 21) {
    chrOrder <- c((1:19), "X", "Y")
  }
  
  dt$chrom <- factor(dt$chrom, chrOrder, ordered = TRUE)
  data.table::setkeyv(dt, c("chrom", "start", "end"))
  return(dt)
}

#' Read chromosome and position information from a bed format file.
#'
#' @param file Path to the file in bed format.
#'
#' @return A data.table keyed by chrom, start, and end. It uses one-based
#'   coordinates.
#'   
#' @keywords internal
ReadBedRanges <- function(file) {
  dt <- data.table::fread(file)
  dt1 <- StandardChromName(dt[, 1:3])
  colnames(dt1) <- c("chrom", "start", "end")

  # Delete duplicate entries in the BED file
  dt2 <- dplyr::distinct(dt1, chrom, start, end, .keep_all = TRUE)

  # Bed file are 0 based start and 1 based end (an oversimplification).
  # We need to add 1L and not 1, otherwise the column turns to a double
  # we get a warning from data.table.
  dt2$start <- dt2$start + 1L

  chrOrder <- c((1:22), "X", "Y")
  dt2$chrom <- factor(dt2$chrom, chrOrder, ordered = TRUE)
  return(data.table::setkeyv(dt2, c("chrom", "start", "end")))
}

#' Create trinucleotide abundance
#'
#' @param file Path to the file with the nucleotide abundance information with 3
#'   base pairs.
#'
#' @return A numeric vector whose names indicate 32 different types of 3 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreateTrinucAbundance <- function(file) {
  dt <- fread(file)
  colnames(dt) <- c("3bp", "occurrences")
  dt$type <-
    ifelse(substr(dt[[1]], 2, 2) %in% c("A", "G"), revc(dt[[1]]), dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create stranded trinucleotide abundance
#'
#' @param file Path to the file with the nucleotide abundance information with 3
#'   base pairs.
#'
#' @return A numeric vector whose names indicate 64 different types of 3 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreateStrandedTrinucAbundance <- function(file) {
  dt <- fread(file)
  colnames(dt) <- c("3bp", "occurrences")
  dt$type <- dt[[1]]
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create dinucleotide abundance
#'
#' @param file Path to the file with the nucleotide abundance information with 2
#'   base pairs.
#'
#' @import data.table
#'
#' @return A numeric vector whose names indicate 10 different types of 2 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreateDinucAbundance <- function(file) {
  dt <- fread(file)
  colnames(dt) <- c("2bp", "occurrences")
  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  dt$type <-
    ifelse((dt[[1]]) %in% canonical.ref, dt[[1]], revc(dt[[1]]))
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create stranded dinucleotide abundance
#'
#' @param file Path to the file with the nucleotide abundance information with 2
#'   base pairs.
#'
#' @import data.table
#'
#' @return A numeric vector whose names indicate 16 different types of 2 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreateStrandedDinucAbundance <- function(file) {
  dt <- fread(file)
  colnames(dt) <- c("2bp", "occurrences")
  dt$type <- dt[[1]]
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create tetranucleotide abundance
#'
#' @param file Path to the file with the nucleotide abundance information with 4
#'   base pairs.
#'
#' @import data.table
#'
#' @return A numeric vector whose names indicate 136 different types of 4 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreateTetranucAbundance <- function(file) {
  dt <- fread(file)
  colnames(dt) <- c("4bp", "occurrences")
  dt$type <- CanonicalizeQUAD(dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create pentanucleotide abundance
#'
#' @param file Path to the file with the nucleotide abundance information
#'   with 5 base pairs.
#'
#' @import data.table
#'
#' @return A numeric vector whose names indicate 512 different types of 5 base
#'   pairs combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreatePentanucAbundance <- function(file) {
  dt <- fread(file)
  colnames(dt) <- c("5bp", "occurrences")
  dt$type <-
    ifelse(substr(dt[[1]], 3, 3) %in% c("A", "G"), revc(dt[[1]]), dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Take strings representing a genome and return the \code{\link{BSgenome}} object.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @return If \code{ref.genome} is 
#' a \code{\link{BSgenome}} object, return it.
#' Otherwise return the \code{\link{BSgenome}} object identified by the
#' string \code{ref.genome}.
#'
#' @keywords internal
NormalizeGenomeArg <- function(ref.genome) {
  
  if (is.null(ref.genome)) stop("Need a non-NULL ref.genome")
  if (class(ref.genome) == "BSgenome") return(ref.genome)
  
  stopifnot(class(ref.genome) == "character")
  
  if (ref.genome %in%
      c("GRCh38", "hg38", "BSgenome.Hsapiens.UCSC.hg38")) {
    if ("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38")) {
      stop("\nPlease install BSgenome.Hsapiens.UCSC.hg38:\n",
           "BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")")
    }
    stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE))
    ref.genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  } else if (ref.genome %in%
             c("GRCh37", "hg19", "BSgenome.Hsapiens.1000genomes.hs37d5")) {
    if ("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5")) {
      stop("\nPlease install BSgenome.Hsapiens.1000genomes.hs37d5:\n",
           "BiocManager::install(\"BSgenome.Hsapiens.1000genomes.hs37d5\")")
    }
    stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE))
    ref.genome <- 
      BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
  } else if (ref.genome %in%
             c("GRCm38", "mm10", "BSgenome.Mmusculus.UCSC.mm10")) {
    if ("" == system.file(package = "BSgenome.Mmusculus.UCSC.mm10")) {
      stop("\nPlease install BSgenome.Mmusculus.UCSC.mm10:\n",
           "BiocManager::install(\"BSgenome.Mmusculus.UCSC.mm10\")")
    }
    stopifnot(requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE))
    ref.genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  } else {
    stop("Unrecoginzed ref.genome:\n", ref.genome,
         "\nNeed a BSgenome reference genome\n",
         "or one of the character strings GRCh38, hg38, GRCh37, hg19, ",
         "GRCm38, mm10")
  }
  
  return(ref.genome)
}

#' Stop if \code{region} is illegal for an in-transcript catalogs
#'
#' @param region The region to test (a character string)
#' 
#' @keywords internal

StopIfTranscribedRegionIllegal <- function(region) {
  if (!region %in% c("transcript", "exome", "unknown")) {
    stop("Require region to be one of transcript, exome, or unknown for ",
         "SBS192Catalog and DBS144Catalog\n",
         "Got ", region)
  }
}

#' Stop if the number of rows in \code{object} is illegal
#' 
#' @param object A \code{catalog}, numeric \code{matrix}, or numeric \code{data.fram}
#' 
#' @keywords internal
StopIfNrowIllegal <- function(object) {
  # Will call stop() if the number or rows is illegal.
  InferCatalogClassString(object)

}

#' Stop if \code{region} is illegal.
#'
#' @param region Character string to check.
#' 
#' @keywords internal
StopIfRegionIllegal <- function(region) {
  if (!region %in% c("genome", "exome", "transcript", "unknown")) {
    stop("Unrecoginzed region identifier: ", region,
         "\nNeed one of genome, exome, transcript, unknown")
  }
 return(NULL)
}

#' Stop if \code{catalog.type} is illegal.
#'
#' @param catalog.type Character string to check.
#'
#' @keywords internal
StopIfCatalogTypeIllegal <- function(catalog.type) {
  if (is.null(catalog.type)) return (NULL)
  if (!catalog.type %in% c("counts", "density",
                           "counts.signature", "density.signature")) {
    stop("Unrecoginzed catalog type identifier: ", catalog.type,
         "\nNeed one of counts, density, counts.signature, density.signature")
  }
  return(NULL)
}

#' Infer the class of catalog in a file.
#'
#' @param file Path to a catalog on disk in the standardized format.
#'
#' @return A character string with class appropriate for the catalog on disk.
#'
#' @keywords internal
InferClassOfCatalogForRead <- function(file) {
  cos <- data.table::fread(file)
  StopIfNrowIllegal(cos)
  class.string <- InferCatalogClassString(cos) 
  return(structure("ClassofCatalog", class = class.string))
}

#' @keywords internal
InferCatalogClassPrefix <- function(object) {
  
  nrow <- nrow(object)
  
  if (nrow == 96)   return("SBS96")
  if (nrow == 192)  return("SBS192")
  if (nrow == 1536) return("SBS1536")
  if (nrow == 78)   return("DBS78")
  if (nrow == 144)  return("DBS144")
  if (nrow == 136)  return("DBS136")
  if (nrow == 83)   return("ID")
  if (nrow == 1697) return("COMPOSITE")
  
  stop("\nThe number of rows in the input object must be one of\n",
       "96 192 1536 78 144 136 83 1697\ngot ", nrow)
  
}


#' @keywords internal
InferCatalogClassString <- function(object) {
  
  prefix <- InferCatalogClassPrefix(object)
  if (prefix == "ID") prefix <- "Indel"
  return(paste0(prefix, "Catalog"))
}


#' Infer the correct rownames for a matrix based on its number of rows
#' 
#' @keywords internal
InferRownames <- function(object) {
  prefix <- InferCatalogClassPrefix(object)
  return(ICAMS::catalog.row.order[[prefix]])
}

#' Check whether the rownames of \code{object} are correct, if yes then put the
#' rows in the correct order.
#'
#' @keywords internal
CheckAndReorderRownames <- function(object) {
  prefix <- InferCatalogClassPrefix(object)
  correct.rownames <- ICAMS::catalog.row.order[[prefix]]
  difference <- setdiff(correct.rownames, rownames(object))
  if(identical(difference, character(0))) {
    return(object[correct.rownames, , drop = FALSE])
  } else {
    stop("\nThe input object does not have correct rownames to denote the\n", 
         "mutation types. Mutation types that are missing are\n", 
         paste(difference, collapse = " "))
  }
}

#' Test if object is \code{BSgenome.Hsapiens.1000genome.hs37d5}.
#'
#' @param x Object to test.
#' 
#' @return TRUE if \code{x} is \code{BSgenome.Hsapiens.1000genome.hs37d5}.
#' 
#' @keywords internal
IsGRCh37 <- function(x) {
  if (is.null(x)) return(FALSE)
  return(NormalizeGenomeArg(x)@pkgname == 
           "BSgenome.Hsapiens.1000genomes.hs37d5")
}

#' Test if object is \code{BSgenome.Hsapiens.UCSC.hg38}.
#'
#' @param x Object to test.
#' 
#' @return TRUE if \code{x} is \code{BSgenome.Hsapiens.UCSC.hg38}.
#' 
#' @keywords internal
IsGRCh38 <- function(x) {
  if (is.null(x)) return(FALSE)
  return(NormalizeGenomeArg(x)@pkgname == 
           "BSgenome.Hsapiens.UCSC.hg38")
}

#' Test if object is \code{BSgenome.Mmusculus.UCSC.mm10}.
#'
#' @param x Object to test.
#' 
#' @return TRUE if \code{x} is \code{BSgenome.Mmusculus.UCSC.mm10}.
#' 
#' @keywords internal
IsGRCm38 <- function(x) {
  if (is.null(x)) return(FALSE)
  return(NormalizeGenomeArg(x)@pkgname == 
           "BSgenome.Mmusculus.UCSC.mm10")
}


#' Infer \code{abundance} given a matrix-like \code{object} and additional information.
#'
#' @param object A numeric matrix, numeric data frame, or \code{catalog}.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @param catalog.type A character string for \code{catalog.type}
#' as described in \code{\link{ICAMS}}.
#'
#' @return A value that can be set as the abundance attribute of
#' a \code{catalog} (which may be \code{NULL} if no abundance
#' can be inferred).
#'
#' @keywords internal
InferAbundance <- function(object, ref.genome, region, catalog.type) {
    
    StopIfNrowIllegal(object)
    StopIfRegionIllegal(region)
    StopIfCatalogTypeIllegal(catalog.type)
    if (nrow(object) == 1697) {
      return(NULL)
      # There are no meaningful abundances for COMPOSITE catalogs.
    }

    if (IsDensity(catalog.type)) {
      ab <- flat.abundance[[as.character(nrow(object))]]
      stopifnot(!is.null(ab))
      return(ab)
    }
    
    if (is.null(ref.genome)) return(NULL)
    ref.genome <- NormalizeGenomeArg(ref.genome)
    ab <- ICAMS::all.abundance[[ref.genome@pkgname]]
    if (is.null(ab)) return(NULL)

    ab2 <- ab[[region]]
    if (is.null(ab2)) return(NULL)

    ab3 <- ab2[[as.character(nrow(object))]]

    return(ab3)

  }

#' Create a catalog from a \code{matrix}, \code{data.frame}, or \code{vector}
#'
#' @param object A numeric \code{matrix}, numeric \code{data.frame},
#' or \code{vector}.
#' If a \code{vector}, converted to a 1-column \code{matrix}
#' with rownames taken from the element names of the \code{vector}
#' and with column name \code{"Unknown"}.
#' If argument \code{infer.rownames}
#'  is \code{FALSE} than this argument must have
#'   rownames to denote the mutation types. See \code{\link{CatalogRowOrder}}
#'   for more details.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region A character string designating a region, one of
#' \code{genome}, \code{transcript}, \code{exome}, \code{unknown};
#' see \code{\link{ICAMS}}.
#'
#' @param catalog.type One of "counts", "density", "counts.signature",
#'   "density.signature".
#'
#' @param abundance If \code{NULL}, then
#'  inferred if \code{ref.genome}
#'  is one of
#'  the reference genomes known to ICAMS and \code{region}
#'  is not \code{unknown}. See \code{\link{ICAMS}}.
#'  The argument \code{abundance} should
#'  contain the counts of different source sequences for mutations
#'  in the same format as the numeric vectors in \code{\link{all.abundance}}.
#'  
#' @param infer.rownames If \code{TRUE}, and \code{object} has no
#' rownames, then assume the rows of \code{object} are in the
#' correct order and add the rownames implied by the number of rows
#' in \code{object} (e.g. rownames for SBS 192 if there are 192 rows).
#' If \code{TRUE}, \strong{be sure the order of rows is correct.}
#'
#' @return A catalog as described in \code{\link{ICAMS}}.
#'
#' @export
#' 
#' @examples
#' # Create an SBS96 catalog with all mutation counts equal to 1.  
#' object <- matrix(1, nrow = 96, ncol = 1, 
#'                  dimnames = list(catalog.row.order$SBS96))
#' catSBS96 <- as.catalog(object)
             
as.catalog <- function(object, 
                       ref.genome = NULL, 
                       region = "unknown", 
                       catalog.type = "counts", 
                       abundance = NULL,
                       infer.rownames = FALSE) {
  if (!is.matrix(object)) {
    if (is.data.frame(object)) {
      object <- as.matrix(object)
    } else if (is.vector(object)) {
      obj2 <- matrix(object, ncol = 1)
      rownames(obj2) <- names(object)
      colnames(obj2) <- "Unknown"
      object <- obj2
    } else {
      stop("object must be numeric matrix, vector, or data frame")
    }
  }
  stopifnot(mode(object) == "numeric")

  if (is.null(rownames(object))) {
    if (!infer.rownames) {
      stop("Require correct rownames on object unless infer.rownames == TRUE")
    }
    rownames(object) <- InferRownames(object)
  } else {
    object <- CheckAndReorderRownames(object)
  }

  StopIfRegionIllegal(region)
  
  # Will call stop() if nrow(object) is illegal
  class.string  <- InferCatalogClassString(object)
  
  StopIfCatalogTypeIllegal(catalog.type)
  
  if (!is.null(ref.genome)) {
    ref.genome <- NormalizeGenomeArg(ref.genome)
  }
  attr(object, "ref.genome") <- ref.genome

  attr(object, "catalog.type") <- catalog.type
  
  if (is.null(abundance)) {
    abundance <- InferAbundance(object, ref.genome, region, catalog.type)
  } 
  attr(object, "abundance") <- abundance

  class(object) <- append(class(object), class.string, after = 0)
  class(object) <- unique(attributes(object)$class)
  
  # TODO(Steve): check that his has been moved 
  # to the VCF -> catalog functions
  if (attributes(object)$class[1] %in% c("SBS192Catalog", "DBS144Catalog")) {
    StopIfTranscribedRegionIllegal(region)
  }
  attr(object, "region") <- region
  
  return(object)
}

#' Generate all possible k-mers of length k.
#'
#' @param k Length of k-mers (k>=2)
#'
#' @return Character vector containing all possible k-mers.
#'
#' @keywords internal
GenerateKmer <- function(k) {
  base <- c("A", "C", "G", "T")
  list.of.base <- list()
  for (i in 1:k){
    list.of.base[[i]] <- base
  }
  permutation <- expand.grid(list.of.base, stringsAsFactors = FALSE)
  all.kmer.list <- character(nrow(permutation))
  for (i in 1:k) {
    all.kmer.list <- stringi::stri_c(all.kmer.list, permutation[[i]])
  }
  all.kmer.list <- stringi::stri_sort(all.kmer.list)
}

#' Generate an empty matrix of k-mer abundance
#'
#' @param k Length of k-mers (k>=2)
#'
#' @return An empty matrix of k-mer abundance
#'
#' @keywords internal
GenerateEmptyKmerCounts <- function(k) {
  all.kmer.list <- GenerateKmer(k)
  kmer.counts <- matrix(0, nrow = length(all.kmer.list))
  rownames(kmer.counts) <- all.kmer.list
  colnames(kmer.counts) <- "Freq"
  return(kmer.counts)
}

#' Generate k-mer abundance from given nucleotide sequences
#'
#' @param sequences A vector of nucleotide sequences
#'
#' @param k Length of k-mers (k>=2)
#'
#' @return Matrix of the counts of each k-mer inside \code{sequences}
#'
#' @keywords  internal
GetSequenceKmerCounts <- function(sequences, k) {
  kmer.counts <- GenerateEmptyKmerCounts(k)

  for(start_idx in 1:k){
    temp.seqs <- substring(sequences, start_idx, nchar(sequences))
    temp.kmers <-
      stringi::stri_extract_all_regex(
        temp.seqs, pattern = paste(rep(".", each = k), collapse = ""))
    temp.kmers <- unlist(temp.kmers)
    temp.kmer.counts <- data.frame(table(temp.kmers))
    row.names(temp.kmer.counts) <- temp.kmer.counts[, 1]
    if (any(grepl("N", temp.kmer.counts[, 1]))) {
      temp.kmer.counts <- temp.kmer.counts[-grep("N", temp.kmer.counts[, 1]), ]
    }

    kmer.counts[row.names(temp.kmer.counts), ] <-
      kmer.counts[row.names(temp.kmer.counts), ]  +
      temp.kmer.counts$Freq
  }
  return(kmer.counts)
}

#' Generate k-mer abundance from a given genome
#'
#' @param k Length of k-mers (k>=2)
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param filter.path If given, homopolymers will be masked from
#'   genome(sequence). Only simple repeat masking is accepted now.
#'   
#' @param verbose If \code{TRUE}, generate progress messages.
#'
#' @importFrom GenomicRanges GRanges
#'
#' @importFrom IRanges IRanges
#'
#' @return Matrix of the counts of each k-mer across the \code{ref.genome}
#'
#' @keywords internal
GetGenomeKmerCounts <- function(k, ref.genome, filter.path, verbose = FALSE) {
  kmer.counts <- GenerateEmptyKmerCounts(k)
  genome <- NormalizeGenomeArg(ref.genome)

  # Remove decoyed chromosomes and mitochondrial DNA
  chr.list <- seqnames(genome)[which(nchar(seqnames(genome)) <= 5)]
  if (any(grepl("M", chr.list))) {
    chr.list <- chr.list[-grep("M", chr.list)]
  }

  if (!missing(filter.path)) {
    filter.df <- fread(filter.path, header = FALSE, stringsAsFactors = FALSE)
    filter.df <- filter.df[filter.df$V6 <= 6]
    filter.df <- StandardChromName(filter.df[, 2:ncol(filter.df)])
    # Check whether chromosome names in filter.df are the same as in ref.genome
    if (!(filter.df$V2[1] %in% seqnames(genome))){
      filter.df$V2 <- paste0("chr", filter.df$V2)
    }
  }

  if (verbose) message("Start counting by chromosomes")

  for (idx in 1:length(chr.list)) {
    if (verbose) message(chr.list[idx])

    if (!missing(filter.path)) {
      chr.filter.df <- filter.df[which(filter.df$V2 == chr.list[idx]), ]
      filter.chr <- with(chr.filter.df, GRanges(V2, IRanges(V3 + 1, V4)))
      genome.ranges.chr <-
        GRanges(chr.list[idx],
                IRanges(1, as.numeric(GenomeInfoDb::seqlengths(genome)[idx])))
      filtered.genome.ranges.chr <- GenomicRanges::setdiff(genome.ranges.chr, filter.chr)
      genome.seq <- BSgenome::getSeq(genome, filtered.genome.ranges.chr,
                                     as.character = TRUE)
      #Filter shorter homopolymer and microsatellites by regex
      genome.seq <- gsub(homopolymer.ms.regex.pattern, "N", genome.seq)

    } else {
      genome.seq <- BSgenome::getSeq(genome, chr.list[idx], as.character = TRUE)
    }

    kmer.counts <- kmer.counts + GetSequenceKmerCounts(genome.seq, k)
  }
  return(kmer.counts)
}

#' Remove ranges that fall on both strands
#'
#' @param stranded.ranges A keyed data table which has stranded ranges information.
#' It has four columns: chrom, start, end and strand.
#'
#' @return A data table which has removed ranges that fall on both strands from
#'   the input \code{stranded.ranges}.
#'
#' @keywords internal
RemoveRangesOnBothStrand <- function(stranded.ranges) {
  dt <- as.data.table(stranded.ranges)
  dt.plus <- dt[strand == "+", ]
  dt.minus <- dt[strand == "-", ]
  gr.plus <-
    GenomicRanges::makeGRangesFromDataFrame(dt.plus, keep.extra.columns = TRUE)
  gr.minus <-
    GenomicRanges::makeGRangesFromDataFrame(dt.minus, keep.extra.columns = TRUE)
  gr.plus.reduced <- GenomicRanges::reduce(gr.plus, with.revmap = TRUE)
  gr.minus.reduced <- GenomicRanges::reduce(gr.minus, with.revmap = TRUE)

  # Find ranges that have transcripts on both strands and remove these
  # ranges from each strand
  gr <-
    GenomicRanges::intersect(gr.plus.reduced, gr.minus.reduced,
                             ignore.strand = TRUE)

  if (length(gr) != 0) {
    gr1 <- gr
    gr1@strand@values <- "+"
    gr2 <- gr
    gr2@strand@values <- "-"
    gr3 <- GenomicRanges::setdiff(gr.plus.reduced, gr1)
    gr4 <- GenomicRanges::setdiff(gr.minus.reduced, gr2)
  } else {
    gr3 <- gr.plus.reduced
    gr4 <- gr.minus.reduced
  }

  # Get a new data table which does not have ranges on both strands
  dt1 <- as.data.table(gr3)
  dt2 <- as.data.table(gr4)
  dt3 <- rbind(dt1, dt2)
  dt4 <- dt3[, c(1:3, 5)]
  colnames(dt4) <- c("chrom", "start", "end", "strand")
  return(setkeyv(dt4, c("chrom", "start", "end")))
}

#' Generate stranded k-mer abundance from a given genome and gene annotation file
#'
#' @param k Length of k-mers (k>=2)
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param filter.path If given, homopolymers will be masked from
#'   genome(sequence). Only simple repeat masking is accepted now.
#'
#' @param stranded.ranges A keyed data table which has stranded ranges
#'   information. It has four columns: chrom, start, end and strand. It should
#'   use one-based coordinate system.
#'   
#' @param verbose If \code{TRUE} generate progress messages.
#'
#' @importFrom stats start end
#'
#' @importFrom GenomicRanges GRanges
#'
#' @importFrom IRanges IRanges
#'
#' @return Matrix of the counts of each stranded k-mer across the \code{ref.genome}
#'
#' @keywords internal
GetStrandedKmerCounts <- 
  function(k, ref.genome, stranded.ranges, filter.path, verbose = FALSE) {
  stranded.ranges <- RemoveRangesOnBothStrand(stranded.ranges)
  stranded.ranges <- StandardChromName(stranded.ranges)
  genome <- NormalizeGenomeArg(ref.genome)
  kmer.counts <- GenerateEmptyKmerCounts(k)

  # Check whether chromosome names in stranded.ranges are the same as in ref.genome
  if (!(stranded.ranges$chrom[1] %in% seqnames(genome))) {
    stranded.ranges$chrom <- paste0("chr", stranded.ranges$chrom)
  }

  if (!missing(filter.path)) {
    filter.df <- fread(filter.path, header = FALSE, stringsAsFactors = FALSE)
    filter.df <- filter.df[filter.df$V6 <= 6]
    filter.df <- StandardChromName(filter.df[, 2:ncol(filter.df)])
    # Check whether chromosome names in filter.df are the same as in ref.genome
    if (!(filter.df$V2[1] %in% seqnames(genome))){
      filter.df$V2 <- paste0("chr", filter.df$V2)
    }
  }

  if (verbose) message("Start counting by chromosomes")

  for (chr in unique(stranded.ranges$chrom)) {
    if (verbose) message(chr)
    temp.stranded.ranges <- stranded.ranges[stranded.ranges$chrom == chr, ]
    stranded.ranges.chr <-
      with(temp.stranded.ranges,
           GRanges(chrom, IRanges(start, end), strand = strand))

    # Remove the overlapping ranges in stranded.ranges.chr
    stranded.ranges.chr <- IRanges::reduce(stranded.ranges.chr)

    if (!missing(filter.path)) {
      chr.filter.df <- filter.df[which(filter.df$V2 == chr), ]

      # Add strand information to chr.filter.df because
      # the setdiff function is strand specific
      filter.chr <-
        c(with(chr.filter.df, GRanges(V2, IRanges(V3 + 1, V4), strand = "+")),
          with(chr.filter.df, GRanges(V2, IRanges(V3 + 1, V4), strand = "-")))
      filtered.stranded.ranges.chr <-
        GenomicRanges::setdiff(stranded.ranges.chr, filter.chr)

      stranded.seq <- BSgenome::getSeq(genome, filtered.stranded.ranges.chr,
                                       as.character = TRUE)
      # Filter shorter homopolymer and microsatellites by regex
      stranded.seq <- gsub(homopolymer.ms.regex.pattern, "N", stranded.seq)
    } else {
      stranded.seq <- BSgenome::getSeq(genome, stranded.ranges.chr,
                                       as.character = TRUE)
    }
    kmer.counts <- kmer.counts + GetSequenceKmerCounts(stranded.seq, k)
  }
  return(kmer.counts)
}

#' Generate custom k-mer abundance from a given reference genome
#'
#' @param k Length of k-mers (k>=2)
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param custom.range A keyed data table which has custom ranges information. It
#'   has three columns: chrom, start and end. It should use one-based coordinate
#'   system. You can use the internal function in this package
#'   \code{ICAMS:::ReadBedRanges} to read a BED file in 0-based coordinates and
#'   convert it to 1-based coordinates.
#'   
#' @param filter.path If given, homopolymers will be masked from
#'   genome(sequence). Only simple repeat masking is accepted now.
#'   
#' @param verbose If \code{TRUE} generate progress messages.
#'
#' @importFrom stats start end
#'
#' @importFrom GenomicRanges GRanges
#'
#' @importFrom IRanges IRanges
#'
#' @return Matrix of the counts of custom k-mer across the \code{ref.genome}
#'
#' @keywords internal
GetCustomKmerCounts <- function(k, ref.genome, custom.ranges, filter.path, 
                                verbose = FALSE) {
  custom.ranges <- StandardChromName(custom.ranges)
  genome <- NormalizeGenomeArg(ref.genome)
  kmer.counts <- GenerateEmptyKmerCounts(k)
  
  # Check whether chromosome names in custom.ranges are the same as in ref.genome
  if (!(custom.ranges$chrom[1] %in% seqnames(genome))) {
    custom.ranges$chrom <- paste0("chr", custom.ranges$chrom)
  }
  
  if (!missing(filter.path)) {
    filter.df <- fread(filter.path, header = FALSE, stringsAsFactors = FALSE)
    filter.df <- filter.df[filter.df$V6 <= 6]
    filter.df <- StandardChromName(filter.df[, 2:ncol(filter.df)])
    # Check whether chromosome names in filter.df are the same as in ref.genome
    if (!(filter.df$V2[1] %in% seqnames(genome))){
      filter.df$V2 <- paste0("chr", filter.df$V2)
    }
    
  }
  if (verbose) message("Start counting by chromosomes")
  
  for (chr in unique(custom.ranges$chrom)) {
    if (verbose) message(chr)
    temp.custom.ranges <- custom.ranges[custom.ranges$chrom == chr, ]
    custom.range.chr <-
      with(temp.custom.ranges, GRanges(chrom, IRanges(start, end)))
    
    # Remove the overlapping ranges in custom.range.chr
    custom.range.chr <- IRanges::reduce(custom.range.chr)
    
    if (!missing(filter.path)) {
      chr.filter.df <- filter.df[which(filter.df$V2 == chr), ]
      filter.chr <- with(chr.filter.df, GRanges(V2, IRanges(V3 + 1, V4)))
      filtered.custom.range.chr <-
        GenomicRanges::setdiff(custom.range.chr, filter.chr)
      custom.seq <- BSgenome::getSeq(genome, filtered.custom.range.chr,
                                     as.character = TRUE)
      #Filter shorter homopolymer and microsatellites by regex
      custom.seq <- gsub(homopolymer.ms.regex.pattern, "N", custom.seq)
      
    } else {
      custom.seq <- BSgenome::getSeq(genome, custom.range.chr,
                                     as.character = TRUE)
    }
    kmer.counts <- kmer.counts + GetSequenceKmerCounts(custom.seq, k)
  }
  return(kmer.counts)
}

# Redefine the [ methods for catalogs
#' @export
`[.SBS96Catalog` <- function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==
                                1) {
  y <- NextMethod("[")
  if (inherits(y, c("integer", "numeric"))) {
    return(y)
  } else {
    class(y) <- class(x)
    for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
      attr(y, at) <- attr(x, at, exact = TRUE)
    }
    return(y)
  }
}

#' @export
`[.SBS192Catalog` <- function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==
                                1) {
  y <- NextMethod("[")
  if (inherits(y, c("integer", "numeric"))) {
    return(y)
  } else {
    class(y) <- class(x)
    for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
      attr(y, at) <- attr(x, at, exact = TRUE)
    }
    return(y)
  }
}

#' @export
`[.SBS1536Catalog` <- function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==
                                 1) {
  y <- NextMethod("[")
  if (inherits(y, c("integer", "numeric"))) {
    return(y)
  } else {
    class(y) <- class(x)
    for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
      attr(y, at) <- attr(x, at, exact = TRUE)
    }
    return(y)
  }
}

#' @export
`[.DBS78Catalog` <- function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==
                                  1) {
  y <- NextMethod("[")
  if (inherits(y, c("integer", "numeric"))) {
    return(y)
  } else {
    class(y) <- class(x)
    for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
      attr(y, at) <- attr(x, at, exact = TRUE)
    }
    return(y)
  }
}

#' @export
`[.DBS144Catalog` <- function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==
                                1) {
  y <- NextMethod("[")
  if (inherits(y, c("integer", "numeric"))) {
    return(y)
  } else {
    class(y) <- class(x)
    for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
      attr(y, at) <- attr(x, at, exact = TRUE)
    }
    return(y)
  }
}

#' @export
`[.DBS136Catalog` <- function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==
                                 1) {
  y <- NextMethod("[")
  if (inherits(y, c("integer", "numeric"))) {
    return(y)
  } else {
    class(y) <- class(x)
    for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
      attr(y, at) <- attr(x, at, exact = TRUE)
    }
    return(y)
  }
}

#' @export
`[.IndelCatalog` <- function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==
                                 1) {
  y <- NextMethod("[")
  if (inherits(y, c("integer", "numeric"))) {
    return(y)
  } else {
    class(y) <- class(x)
    for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
      attr(y, at) <- attr(x, at, exact = TRUE)
    }
    return(y)
  }
}

#' @keywords internal
CheckAndAssignAttributes <- function(x, list0) {
  for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
    
    GetAttribute2 <- function(item, at) {
      attr(item, at, exact = TRUE)
    }
    attributes.list <- lapply(list0, FUN = GetAttribute2, at = at)
    
    CheckNullAttribute <- function(x, attributes.list) {
      is.null.result <- sapply(attributes.list, FUN = is.null)
      # If all the attributes are NULL, then x will also have NULL attribute
      if (all(is.null.result)) {
        attr(x, at) <- NULL
        return(list(is.null = TRUE, x = x))
      } else if (any(is.null.result)) {
        # If some attributes are NULL, then x will have a "mixed" attribute
        attr(x, at) <- "mixed"
        return(list(is.null = TRUE, x = x))
      }
    }
    
    if (at == "ref.genome") {
      result.list <- CheckNullAttribute(x, attributes.list)
      if (isTRUE(result.list$is.null)) {
        x <- result.list$x
      } else {
        GetRefGenomeName <- function(object) {
          return(object@pkgname)
        }
        ref.genome.check.result <- 
          sapply(attributes.list, FUN = GetRefGenomeName)
        if (length(unique(ref.genome.check.result)) == 1) {
          attr(x, at) <- attributes.list[[1]]
        } else {
          attr(x, at) <- "mixed"
        }
      }
    }
    
    if (at == "catalog.type") {
      is.null.result <- sapply(attributes.list, FUN = is.null)
      if (any(is.null.result)) {
        # If any attribute is NULL, then stop
        stop("Cannot perform cbind operation to catalogs with NULL",
             " catalog.type")
      } else {
        catalog.type.check.result <- do.call("c", attributes.list)
        if (length(unique(catalog.type.check.result)) == 1) {
          attr(x, at) <- attributes.list[[1]]
        } else {
          stop("Cannot perform cbind operation to catalogs with different",
               " catalog.type")
        }
      }
    }
    
    if (at == "abundance") {
      result.list <- CheckNullAttribute(x, attributes.list)
      if (isTRUE(result.list$is.null)) {
        x <- result.list$x
      } else {
        abundance.length.result <- sapply(attributes.list, FUN = length)
        if (length(unique(abundance.length.result)) == 1) {
          GetAbundanceNames <- function(abundance) {
            return(sort(names(abundance)))
          }
          
          abundance.names.list <- 
            lapply(attributes.list, FUN = GetAbundanceNames)
          CheckNamesOfAbundance <- function(abundance.names.list) {
            num <- length(abundance.names.list)
            ref.names <- abundance.names.list[[1]]
            for (i in 2:num) {
              if (!all(abundance.names.list[[i]] == ref.names)) {
                return(FALSE)
              } 
            }
            return(TRUE)
          }
          
          CheckValuesOfAbundance <- function(abundance.list) {
            num <- length(abundance.list)
            ref.values <- abundance.list[[1]][abundance.names.list[[1]]]
            for (i in 2:num) {
              if (!all(abundance.list[[i]][abundance.names.list[[i]]] == 
                       ref.values)) {
                return(FALSE)
              }
            }
            return(TRUE)
          }
          
          if (CheckNamesOfAbundance(abundance.names.list) && 
              CheckValuesOfAbundance(attributes.list)) {
            attr(x, at) <- attributes.list[[1]]
          } else {
            attr(x, at) <- "mixed"
          }
        } else {
          attr(x, at) <- "mixed"
        }
      }
    }
    
    if (at == "region") {
      result.list <- CheckNullAttribute(x, attributes.list)
      if (isTRUE(result.list$is.null)) {
        x <- result.list$x
      } else {
        region.check.result <- do.call("c", attributes.list)
        if (length(unique(region.check.result)) == 1) {
          attr(x, at) <- attributes.list[[1]]
        } else {
          attr(x, at) <- "mixed"
        }
      }
    }
  }
  return(x)
}

# Redefine the cbind methods for catalogs
#' @export
`cbind.SBS96Catalog` <- function (..., deparse.level = 1) {
  x <- base::cbind.data.frame(..., deparse.level = deparse.level)
  x <- data.matrix(x)
  class(x) <- class(..1)
  list0 <- list(...)
  x <- CheckAndAssignAttributes(x, list0)
  return(x)
}

#' @export
`cbind.SBS192Catalog` <- function (..., deparse.level = 1) {
  x <- base::cbind.data.frame(..., deparse.level = deparse.level)
  x <- data.matrix(x)
  class(x) <- class(..1)
  list0 <- list(...)
  x <- CheckAndAssignAttributes(x, list0)
  return(x)
}

#' @export
`cbind.SBS1536Catalog` <- function (..., deparse.level = 1) {
  x <- base::cbind.data.frame(..., deparse.level = deparse.level)
  x <- data.matrix(x)
  class(x) <- class(..1)
  list0 <- list(...)
  x <- CheckAndAssignAttributes(x, list0)
  return(x)
}

#' @export
`cbind.DBS78Catalog` <- function (..., deparse.level = 1) {
  x <- base::cbind.data.frame(..., deparse.level = deparse.level)
  x <- data.matrix(x)
  class(x) <- class(..1)
  list0 <- list(...)
  x <- CheckAndAssignAttributes(x, list0)
  return(x)
}

#' @export
`cbind.DBS144Catalog` <- function (..., deparse.level = 1) {
  x <- base::cbind.data.frame(..., deparse.level = deparse.level)
  x <- data.matrix(x)
  class(x) <- class(..1)
  list0 <- list(...)
  x <- CheckAndAssignAttributes(x, list0)
  return(x)
}

#' @export
`cbind.DBS136Catalog` <- function (..., deparse.level = 1) {
  x <- base::cbind.data.frame(..., deparse.level = deparse.level)
  x <- data.matrix(x)
  class(x) <- class(..1)
  list0 <- list(...)
  x <- CheckAndAssignAttributes(x, list0)
  return(x)
}

#' @export
`cbind.IndelCatalog` <- function (..., deparse.level = 1) {
  x <- base::cbind.data.frame(..., deparse.level = deparse.level)
  x <- data.matrix(x)
  class(x) <- class(..1)
  list0 <- list(...)
  x <- CheckAndAssignAttributes(x, list0)
  return(x)
}

#' Calculate base counts from three mer abundance
#' @keywords internal
CalBaseCountsFrom3MerAbundance <- function(three.mer.abundance) {
  base.counts <- integer(4)
  names(base.counts) <- c("A", "C", "G", "T")
  
  tmp <- three.mer.abundance
  names(tmp) <- substr(names(tmp), 2, 2)
  for (base in c("A", "C", "G", "T")) {
    base.counts[base] <- sum(tmp[names(tmp) == base])
  }
  return(base.counts)
}

#' @keywords internal
IsBinomialTestApplicable <- function(catalog) {
  catalog.type <- attributes(catalog)$catalog.type
  abundance <- attributes(catalog)$abundance
  
  if(catalog.type == "counts" && !is.null(abundance)) {
    if (length(abundance) == 64) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(FALSE)
  }
}

#' @keywords internal
AssignNumberOfAsterisks <- function(value) {
  label <- NULL
  if (value < 0.001) {
    label <- "***"
  } else if (value < 0.01) {
    label <- "**"
  } else if (value < 0.05) {
    label <- "*"
  }
  return(label)
}