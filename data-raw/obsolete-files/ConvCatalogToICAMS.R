#' Methods to convert to ICAMS-formatted catalog
#'
#' See the generic function \code{\link{ConvCatalogToICAMS}}
#'
#' @inheritParams ConvCatalogToICAMS
#' 
#' @title ConvCatalogToICAMS.Methods
#' @name ConvCatalogToICAMS.Methods
#' 
NULL


#' Convert a COSMIC or SigProfiler-formatted catalog to ICAMS-format
#'
#' Read a catalog in standardized format from path.
#'
#' @param file Path to a catalog on disk in the standardized format.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @param catalog.type One of "counts", "density", "counts.signature",
#'   "density.signature".
#'
#' @param strict If TRUE, do additional checks on the input, and stop if the
#'   checks fail.
#'   
#' @param stop.on.error If TRUE, call \code{stop} on error; otherwise
#'   return a 1-column matrix of NA's with the attribute "error"
#'   containing error information. The number of rows may not
#'   be the correct number for the expected catalog type.
#'
#' @return A catalog as an S3 object; see \code{\link{as.catalog}}.
#'
#' @note In ID (small insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#' 
#' @examples 
#' NULL
#' @export
ConvCatalogToICAMS <- function(
  file, ref.genome = NULL, region = "unknown", 
  catalog.type = "counts", strict = TRUE,
  stop.on.error = TRUE) {
  
    tryCatch(
    return(ConvCatalogToICAMSInternal(
      file = file,
      ref.genome = ref.genome,
      region = region,
      catalog.type = catalog.type,
      strict = strict,
      stop.on.error = stop.on.error)),
        error = function(err) {
      return(
        ReadCatalogErrReturn(
          err.info      = err,
          nrow          = 1, # We do not know what type of catalog
          stop.on.error = stop.on.error))})
}

#' Internal catalog conversion function to be wrapped in a tryCatch
#' @inheritParams ConvCatalogToICAMS
#' @keywords internal
ConvCatalogToICAMSInternal <- function(
  file, ref.genome = NULL, region = "unknown", 
  catalog.type = "counts", strict = TRUE,
  stop.on.error = TRUE) {
  #StopIfRegionIllegal(region)
  #StopIfCatalogTypeIllegal(catalog.type)
  class.of.catalog <- InferClassOfCatalogForReadNonICAMS(file)
  UseMethod(generic = "ConvCatalogToICAMS", 
            object = class.of.catalog)
}





#' @rdname ConvCatalogToICAMS.Methods
#' @export
ConvCatalogToICAMS.COSMICSBS96Catalog <- function(
  file, ref.genome = NULL, region = "unknown", 
  catalog.type = "counts", strict = TRUE,
  stop.on.error = TRUE) {
  cos <- data.table::fread(file)
  stopifnot(nrow(cos) == 96)
  if (strict) {
    stopifnot(names(cos)[1] == "Type")
    stopifnot(names(cos)[2] == "Subtype")
  }
  ref.gt.var       <- unlist(cos[, 1])
  before.ref.after <- unlist(cos[, 2])
  var <- substring(ref.gt.var, 3, 3)
  out <- cos[, -(1 : 2), drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- paste0(before.ref.after, var)
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[ICAMS::catalog.row.order$SBS96, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}


#' @rdname ConvCatalogToICAMS.Methods
#' @export
ConvCatalogToICAMS.SigProfilerSBS96Catalog <- function(
  file, ref.genome = NULL, region = "unknown", 
  catalog.type = "counts", strict = TRUE,
  stop.on.error = TRUE) {
  cos <- data.table::fread(file)
  stopifnot(nrow(cos) == 96)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("MutationType", "MutationsType"))
  }
  before <- substring(unlist(cos[, 1]),1,1)
  ref <- substring(unlist(cos[, 1]),3,3)
  after <- substring(unlist(cos[, 1]),7,7)
  before.ref.after <- paste0(before,ref,after)
  var <- substring(unlist(cos[, 1]), 5, 5)
  out <- cos[, -1, drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- paste0(before.ref.after, var)
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[ICAMS::catalog.row.order$SBS96, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}



#' @rdname ConvCatalogToICAMS.Methods
#' @export
ConvCatalogToICAMS.COSMICSBS192Catalog <-
  function(file, ref.genome = NULL, region = "unknown", 
           catalog.type = "counts", strict = TRUE,
           stop.on.error = TRUE) {
  
  StopIfTranscribedRegionIllegal(region)
  
  cos <- data.table::fread(file)
  # cos.copy <- cos # For debugging, testing
  stopifnot(nrow(cos) == 192)
  if (strict) {
    stopifnot(names(cos)[2] == "Type")
    stopifnot(names(cos)[3] == "Subtype")
    stopifnot(names(cos)[1] == "Strand")
  }
  ref.gt.var       <- unlist(cos[, 2])
  before.ref.after <- unlist(cos[, 3])

  ## Find the rows labeled with "T", indicating the
  ## SBS is on the transcribed (which is the *antisense*) strand.
  transcribed.strand.pos <- which(cos[, 1] == 'T')

  before.ref.after[transcribed.strand.pos] <-
    revc(before.ref.after[transcribed.strand.pos])

  var <- substring(ref.gt.var, 3, 3)
  var[transcribed.strand.pos] <- revc(var[transcribed.strand.pos])

  tmp <- paste0(before.ref.after, var)
  out <- cos[, -(1 : 3), drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- tmp
  out <- out[ICAMS::catalog.row.order$SBS192, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}

#' @rdname ConvCatalogToICAMS.Methods
#' @export
ConvCatalogToICAMS.SigProfilerSBS1536Catalog <- function(
  file, ref.genome = NULL, region = "unknown", 
  catalog.type = "counts", strict = TRUE,
  stop.on.error = TRUE) {
  cos <- data.table::fread(file)
  stopifnot(nrow(cos) == 1536)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("MutationType", "MutationsType"))
  }
  
  ## Convert rownames to ICAMS format
  ref.gt.var <- substring(unlist(cos[,1]), 4, 6)
  ref <- substring(ref.gt.var, 1, 1)
  var <- substring(ref.gt.var, 3, 3)
  before <- substring(unlist(cos[,1]), 1, 2)
  after <- substring(unlist(cos[,1]), 8, 9)
  out <- as.matrix(cos[ , -1], drop = FALSE)
  rownames(out) <- paste0(before, ref, after, var)

  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[ICAMS::catalog.row.order$SBS1536, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}


#' @rdname ConvCatalogToICAMS.Methods
#' @export
ConvCatalogToICAMS.COSMICDBS78Catalog <- 
  function(file, ref.genome = NULL, region = "unknown", 
           catalog.type = "counts", strict = TRUE,
           stop.on.error = TRUE) {
  cos <- data.table::fread(file)
  stopifnot(nrow(cos) == 78)
  if (strict) {
    stopifnot(names(cos)[1] == "Type")
  }
  Ref <- substring(unlist(cos[,1]), 1, 2)
  Var <- substring(unlist(cos[,1]), 4, 5)

  out <- cos[, -1, drop = FALSE]
  out <- as.matrix(out)
  if(ncol(out) == 1) colnames(out) <- colnames(cos)[2]
  
  rn <- paste0(cos$Ref, cos$Var)
  diff1 <- sort(setdiff(rn, ICAMS::catalog.row.order$DBS78))
  if ( (length(diff1) > 0)
       &&
       (diff1 == c("CGAA", "CGAC", "CGGA", "TAAC", "TAAG", "TACC"))
       &&
       (sort(setdiff(ICAMS::catalog.row.order$DBS78, rn) ==
             c("CGGT", "CGTC", "CGTT", "TACT", "TAGG", "TAGT")))
  ) {
    warning("using temporary hack to handle old DBS canonicalization")
    # CGAA -> CGTT
    rn[rn == "CGAA"] <- "CGTT"

    # CGAC -> CGGT
    rn[rn == "CGAC"] <- "CGGT"

    # CGGA -> CGTC
    rn[rn == "CGGA"] <- "CGTC"

    # TAAC -> TAGT
    rn[rn == "TAAC"] <- "TAGT"

    # TAAG -> TACT
    rn[rn == "TAAG"] <- "TACT"

    # TACC -> TAGG
    rn[rn == "TACC"] <- "TAGG"
  }
  rownames(out) <- rn
  out <- out[ICAMS::catalog.row.order$DBS78, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}

#' @rdname ConvCatalogToICAMS.Methods
#' @export
ConvCatalogToICAMS.SigProfilerDBS78Catalog <- 
  function(file, ref.genome = NULL, region = "unknown", 
           catalog.type = "counts", strict = TRUE,
           stop.on.error = TRUE) {
    cos <- data.table::fread(file)
    stopifnot(nrow(cos) == 78)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("MutationType", "MutationsType"))
  }

  Ref <- substring(unlist(cos[,1]), 1, 2)
  Var <- substring(unlist(cos[,1]), 4, 5)

  out <- cos[, -1, drop = FALSE]
  out <- as.matrix(out)
  if(ncol(out) == 1) colnames(out) <- colnames(cos)[2]


  rn <- paste0(cos$Ref, cos$Var)
  diff1 <- sort(setdiff(rn, ICAMS::catalog.row.order$DBS78))
  if ( (length(diff1) > 0)
       &&
       (diff1 == c("CGAA", "CGAC", "CGGA", "TAAC", "TAAG", "TACC"))
       &&
       (sort(setdiff(ICAMS::catalog.row.order$DBS78, rn) ==
             c("CGGT", "CGTC", "CGTT", "TACT", "TAGG", "TAGT")))
  ) {
    warning("using temporary hack to handle old DBS canonicalization")
    # CGAA -> CGTT
    rn[rn == "CGAA"] <- "CGTT"

    # CGAC -> CGGT
    rn[rn == "CGAC"] <- "CGGT"

    # CGGA -> CGTC
    rn[rn == "CGGA"] <- "CGTC"

    # TAAC -> TAGT
    rn[rn == "TAAC"] <- "TAGT"

    # TAAG -> TACT
    rn[rn == "TAAG"] <- "TACT"

    # TACC -> TAGG
    rn[rn == "TACC"] <- "TAGG"
  }
  rownames(out) <- rn
  out <- out[ICAMS::catalog.row.order$DBS78, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}

#' @rdname ConvCatalogToICAMS.Methods
#' @export
ConvCatalogToICAMS.COSMICSigProfilerDBS78Catalog <- 
  function(file, ref.genome = NULL, region = "unknown", 
           catalog.type = "counts", strict = TRUE,
           stop.on.error = TRUE) {
    cos <- data.table::fread(file)
    stopifnot(nrow(cos) == 78)
    if (strict) {
      stopifnot(names(cos)[1] %in% c("Type", "MutationType", "MutationsType"))
    }
    
    Ref <- substring(unlist(cos[,1]), 1, 2)
    Var <- substring(unlist(cos[,1]), 4, 5)
    
    out <- cos[, -1, drop = FALSE]
    out <- as.matrix(out)
    if(ncol(out) == 1) colnames(out) <- colnames(cos)[2]
    
    
    rn <- paste0(Ref, Var)
    diff1 <- sort(setdiff(rn, ICAMS::catalog.row.order$DBS78))
    if ( (length(diff1) > 0)
         &&
         (diff1 == c("CGAA", "CGAC", "CGGA", "TAAC", "TAAG", "TACC"))
         &&
         (sort(setdiff(ICAMS::catalog.row.order$DBS78, rn) ==
               c("CGGT", "CGTC", "CGTT", "TACT", "TAGG", "TAGT")))
    ) {
      warning("using temporary hack to handle old DBS canonicalization")
      # CGAA -> CGTT
      rn[rn == "CGAA"] <- "CGTT"
      
      # CGAC -> CGGT
      rn[rn == "CGAC"] <- "CGGT"
      
      # CGGA -> CGTC
      rn[rn == "CGGA"] <- "CGTC"
      
      # TAAC -> TAGT
      rn[rn == "TAAC"] <- "TAGT"
      
      # TAAG -> TACT
      rn[rn == "TAAG"] <- "TACT"
      
      # TACC -> TAGG
      rn[rn == "TACC"] <- "TAGG"
    }
    rownames(out) <- rn
    out <- out[ICAMS::catalog.row.order$DBS78, , drop = FALSE]
    return(as.catalog(out, ref.genome, region, catalog.type))
  }

#' @rdname ConvCatalogToICAMS.Methods
#' @export
ConvCatalogToICAMS.COSMICSigProfilerID83Catalog <- function(file, ref.genome = NULL, region = "unknown", 
                                     catalog.type = "counts", strict = TRUE,
                                     stop.on.error = TRUE) {
  # null.out <- matrix(NA, ncol = 1, nrow = length(ICAMS::catalog.row.order$ID))
  cos <- data.table::fread(file)
  ConvCatalogToICAMS.COSMICSigProfilerID83CatalogDT(catalog.fread.dt = cos, 
                                                    ref.genome = ref.genome, 
                                                    region = region, 
                                                    catalog.type = catalog.type, 
                                                    strict = strict,
                                                    stop.on.error = stop.on.error)
}


#' @keywords internal
ConvCatalogToICAMS.COSMICSigProfilerID83CatalogDT <- function(
    catalog.fread.dt, ref.genome = NULL, region = "unknown", 
    catalog.type = "counts", strict = TRUE,
    stop.on.error = TRUE) {
    tryCatch({
    cos <- catalog.fread.dt
    if (nrow(cos) != 83) {
      stop("Expected 83 rows in catalog.fread.dt, got ", nrow(cos))
    }
    
    if (any(grepl("Del:M:1", cos[ , 1]))) {
      message("Interpreting catalog.fread.dt as a COSMIC/SigProfiler insertion/deletion catalog")
      rn <- TransRownames.ID.SigPro.PCAWG(unlist(cos[ , 1]))
      out <- as.matrix(cos[ , -1, drop = FALSE])
    } else {
      stop("catalog.fread.dt is not a COSMIC/SigProfiler ID83 catalog")
    }
    
    # null.out <- matrix(NA, ncol = ncol(out), nrow = nrow(out))
    if (!setequal(rn, ICAMS::catalog.row.order$ID)) {
      msg <- 
        paste("The row names are not correct:\n",
              "got", paste(rn, collapse = ", "),
              "\nexpected", paste(ICAMS::catalog.row.order$ID,
                                  collapse = ", "))
      stop(msg)
    }
    rownames(out) <- rn
    #   if (ncol(out) == 1) colnames(out) <- colnames(cos)[3] 
    out <- out[ICAMS::catalog.row.order$ID, , drop = FALSE]
    return(as.catalog(out, ref.genome, region, catalog.type))
  },
  error = function(e) { 
    return(
      ReadCatalogErrReturn(e, 
                           nrow = length(ICAMS::catalog.row.order$ID),
                           stop.on.error = stop.on.error))}
  ) # tryCatch


  
}

									 
#' @rdname ConvCatalogToICAMS.Methods
#' @export
ConvCatalogToICAMS.COSMICSigProfilerID96Catalog <- function(file, ref.genome = NULL, region = "unknown", 
                                     catalog.type = "counts", strict = TRUE,
                                     stop.on.error = TRUE) {
  
  # null.out <- matrix(NA, ncol = 1, nrow = length(ICAMS::catalog.row.order$ID))
  cos <- data.table::fread(file)
  
  if(nrow(cos) != 96) stop("Expected 96 rows in catalog file, got ", nrow(cos))

  if (cos[84,1] == "2:Ins:M:1" & nrow(cos) == 96) {
    out <- cos[-(84:96), , drop = FALSE]
	ConvCatalogToICAMS.COSMICSigProfilerID83CatalogDT(catalog.fread.dt = out, 
	                                                  ref.genome = ref.genome, 
	                                                  region = region, 
	                                                  catalog.type = catalog.type, 
	                                                  strict = strict,
	                                                  stop.on.error = stop.on.error)
  } else {
    stop("Not an COSMIC/SigProfiler-formatted ID96 catalog.")
  }

}

