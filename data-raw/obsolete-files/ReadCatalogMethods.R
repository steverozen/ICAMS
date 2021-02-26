#' Read catalog methods
#'
#' See the generic function \code{\link{ReadCatalog}}
#'
#' @inheritParams ReadCatalog
#' 
#' @title ReadCatalog.Methods
#' @name ReadCatalog.Methods
#' 
NULL

#' @rdname ReadCatalog.Methods
#' @keywords internal
ReadCatalog.SBS96Catalog <- function(file, ref.genome = NULL, region = "unknown", 
                                     catalog.type = "counts", strict = TRUE,
                                     stop.on.error = TRUE) {
  cos <- data.table::fread(file)
  stopifnot(nrow(cos) == 96)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("Mutation type", "Mutation Type",
                                   "Mutation.type", "Mutation.Type"))
    stopifnot(names(cos)[2] == "Trinucleotide")
  }
  ref.gt.var       <- unlist(cos[, 1])
  before.ref.after <- unlist(cos[, 2])
  var <- substring(ref.gt.var, 3, 3)
  out <- cos[, -(1 : 2), drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- paste0(before.ref.after, var)
  if (strict) {
    stopifnot(rownames(out) == ICAMS::catalog.row.order$SBS96)
  }
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[ICAMS::catalog.row.order$SBS96, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}

#' @keywords internal
ReadSBS96CatalogFromTsv <- function(file, ref.genome = NULL, region = "unknown", 
                                    catalog.type = "counts", strict = TRUE) {
  cos <- data.table::fread(file)
  stopifnot(nrow(cos) == 96)
  if (strict) {
    stopifnot(names(cos)[1] == "Bef")
    stopifnot(names(cos)[2] == "Ref")
    stopifnot(names(cos)[3] == "After")
    stopifnot(names(cos)[4] == "Var")
  }
  before.ref.after <- 
    paste0(unlist(cos[, 1]), unlist(cos[, 2]), unlist(cos[, 3]))
  var <- unlist(cos[, 4])
  out <- cos[, -(1 : 4), drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- paste0(before.ref.after, var)
  if (strict) {
    stopifnot(rownames(out) == ICAMS::catalog.row.order$SBS96)
  }
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[5]
  out <- out[ICAMS::catalog.row.order$SBS96, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}

#' @rdname ReadCatalog.Methods
#' @keywords internal
ReadCatalog.SBS192Catalog <- function(file, ref.genome = NULL, region = "unknown", 
                                      catalog.type = "counts", strict = TRUE,
                                      stop.on.error = TRUE) {
  
  StopIfTranscribedRegionIllegal(region)
  
  cos <- data.table::fread(file)
  # cos.copy <- cos # For debugging, testing
  stopifnot(nrow(cos) == 192)
  if (strict) {
    stopifnot(names(cos)[2] %in% c("Mutation type", "Mutation.type"))
    stopifnot(names(cos)[3] == "Trinucleotide")
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
  if (strict) {
    stopifnot(tmp == ICAMS::catalog.row.order$SBS192)
  }
  out <- cos[, -(1 : 3), drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- tmp
  out <- out[ICAMS::catalog.row.order$SBS192, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}

#' @rdname ReadCatalog.Methods
#' @keywords internal
ReadCatalog.SBS1536Catalog <- function(file, ref.genome = NULL, region = "unknown", 
                                       catalog.type = "counts", strict = TRUE,
                                       stop.on.error = TRUE) {
  cos <- data.table::fread(file)
  stopifnot(nrow(cos) == 1536)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("Mutation type", "Mutation Type",
                                   "Mutation.type", "Mutation.Type"))
    stopifnot(names(cos)[2] == "Pentanucleotide")
  }
  names(cos)[1:2] <- c("Mutation type", "Pentanucleotide")
  ref.gt.var       <- cos[["Mutation type"]]
  before.ref.after <- cos[["Pentanucleotide"]]
  var <- substring(ref.gt.var, 3, 3)
  out <- as.matrix(cos[ , -(1 : 2)], drop = FALSE)
  rownames(out) <- paste0(before.ref.after, var)
  if (strict) {
    stopifnot(rownames(out) == ICAMS::catalog.row.order$SBS1536)
  }
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[ICAMS::catalog.row.order$SBS1536, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}


#' @rdname ReadCatalog.Methods
#' @keywords internal
ReadCatalog.DBS78Catalog <- 
  function(file, ref.genome = NULL, region = "unknown", 
           catalog.type = "counts", strict = TRUE,
           stop.on.error = TRUE) {
    cos <- data.table::fread(file)
    stopifnot(nrow(cos) == 78)
    if (strict) {
      stopifnot(names(cos)[1 : 2] == c("Ref", "Var"))
    }
    names(cos)[1 : 2] <- c("Ref", "Var")
    out <- cos[, -(1 : 2), drop = FALSE]
    out <- as.matrix(out)
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
    if (strict) {
      stopifnot(rownames(out) == ICAMS::catalog.row.order$DBS78)
    }
    out <- out[ICAMS::catalog.row.order$DBS78, , drop = FALSE]
    return(as.catalog(out, ref.genome, region, catalog.type))
  }


#' @rdname ReadCatalog.Methods
#' @keywords internal
ReadCatalog.DBS144Catalog <- function(file, ref.genome = NULL, region = "unknown", 
                                      catalog.type = "counts", strict = TRUE,
                                      stop.on.error = TRUE) {
  
  StopIfTranscribedRegionIllegal(region)
  
  cos <- data.table::fread(file)
  stopifnot(nrow(cos) == 144)
  if (strict) {
    stopifnot(names(cos)[1 : 2] == c("Ref", "Var"))
  }
  names(cos)[1 : 2] <- c("Ref", "Var")
  out <- cos[, -(1 : 2), drop = FALSE]
  out <- as.matrix(out)
  rn <- paste0(cos$Ref, cos$Var)
  rownames(out) <- rn
  if (strict) {
    stopifnot(rownames(out) == ICAMS::catalog.row.order$DBS144)
  }
  out <- out[ICAMS::catalog.row.order$DBS144, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}


#' @rdname ReadCatalog.Methods
#' @keywords internal
ReadCatalog.DBS136Catalog <- function(file, ref.genome = NULL, region = "unknown", 
                                      catalog.type = "counts", strict = TRUE,
                                      stop.on.error = TRUE) {
  cos <- data.table::fread(file)
  stopifnot(nrow(cos) == 136)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("Quad", "quad", "QUAD"))
  }
  names(cos)[1] <- "Quad"
  out <- cos[, -1, drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- cos$Quad
  if (strict) {
    stopifnot(rownames(out) == ICAMS::catalog.row.order$DBS136)
  }
  out <- out[ICAMS::catalog.row.order$DBS136, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}


#' @rdname ReadCatalog.Methods
#' @keywords internal
ReadCatalog.IndelCatalog <- function(file, ref.genome = NULL, region = "unknown", 
                                     catalog.type = "counts", strict = TRUE,
                                     stop.on.error = TRUE) {
  
  tryCatch({
    # null.out <- matrix(NA, ncol = 1, nrow = length(ICAMS::catalog.row.order$ID))
    cos <- data.table::fread(file)
    
    
    if (nrow(cos) != 83) {
      stop("Expected 83 rows in catalog file, got ", nrow(cos))
    }
    
    if (any(grepl("Del:M:1", cos[ , 1]))) {
      if (strict) {
        stop("Cannot interpret ", file, 
             " as a SigProfiler ID catalog when strict = TRUE")
      } 
      warning("Interpreting ", file, 
              " as a SigProfiler insertion/deletion catalog")
      rn <- TransRownames.ID.SigPro.PCAWG(unlist(cos[ , 1]))
      out <- as.matrix(cos[ , -1, drop = FALSE])
    } else {
      cn <- names(cos)
      ex.cn <- c("Type", "Subtype", "Indel_size", "Repeat_MH_size")
      # Repeat_MH_size is the size of repeat OR microhomology (MH)
      # if (strict) { for (i in 1:4) { stopifnot(cn[i] == ex.cn[i]) } }
      if (strict) stopifnot(cn[1:4] == ex.cn)
      names(cos)[1:4] <- ex.cn
      rn <- apply(cos[ , 1:4], MARGIN = 1, paste, collapse = ":")
      out <- as.matrix(cos[ , -(1:4), drop = FALSE])
    }
    
    # null.out <- matrix(NA, ncol = ncol(out), nrow = nrow(out))
    if ((length(setdiff(rn, ICAMS::catalog.row.order$ID)) > 0) ||
        (length(setdiff(ICAMS::catalog.row.order$ID, rn)) > 0)) {
      msg <- 
        paste("The row names are not correct:\n",
              "got", paste(rn, collapse = ", "),
              "\nexpected", paste(ICAMS::catalog.row.order$ID,
                                  collapse = ", "))
      stop(msg)
    }
    if (strict) {
      stopifnot(rn == ICAMS::catalog.row.order$ID)
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

#' @rdname ReadCatalog.Methods
ReadCatalog.COMPOSITECatalog <-
  function(file, ref.genome = NULL, region = "unknown", 
           catalog.type = "counts", strict = TRUE,
           stop.on.error = TRUE)
    
  {
    dt <- data.table::fread(file)
    names <- dt[[1]]
    dt1 <- dt[, -1]
    mat <- as.matrix(dt1)
    rownames(mat) <- names
    return(as.catalog(mat,
                      ref.genome = ref.genome,
                      region = region,
                      catalog.type = catalog.type))
  }