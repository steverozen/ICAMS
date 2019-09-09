#' Read a 192-channel spectra (or signature) catalog in Duke-NUS format.
#' 
#' WARNING: will not work with \code{region = "genome"}. For this 
#' you must first read with \code{region = "unknown"}, then
#' convert the \code{cat96} return to \code{"genome"} and
#' ignore the \code{cat192} return, which is nonsensical.
#' 
#' The file needs to have the column names Before	Ref	After	Var
#' in the first 4 columns 
#' @keywords internal
#' @return A list with two elements
#' @importFrom utils read.table

ReadDukeNUSCat192 <- function(file,
                              ref.genome   = NULL,
                              region       = "unknown",
                              catalog.type = "counts",
                              abundance    = NULL) {
  df <- read.table(file,
                   stringsAsFactors = FALSE,
                   as.is            = TRUE,
                   header           = TRUE, 
                   check.names      = FALSE)
  # Careful, df has 192 rows
  stopifnot(nrow(df)==192)
  
  rownames(df) <- paste(df$Before, df$Ref, df$After, df$Var, sep='')
  df <- df[ , -(1:4)]
  df <- df[ICAMS::catalog.row.order$SBS192,]
  cat192 <- as.catalog(df, ref.genome = ref.genome, region = region,
                       catalog.type = catalog.type, abundance = abundance)
  cat96  <- Collapse192CatalogTo96(cat192)
  
  list(cat96 = cat96, cat192 = cat192)
}

#' Convert 96-channel mutations-type identifiers like this \code{"A[C>A]T" -> "ACTA"}.
#' 
#' @param c1 A vector of character strings with the mutation indicated by
#' e.g. \code{[C>A]} in the middle.
#' 
#' @keywords internal
Unstaple96 <- function(c1) {
  retval <-
    paste0(substr(c1, 1, 1),
           substr(c1, 3, 3),
           substr(c1, 7, 7),
           substr(c1, 5, 5))
  return(retval)
}

#' Convert 96-channel mutation-type identifiers like this \code{"ACTA" -> "A[C>A]T"}.
#' 
#' This is an internal function needed for generating
#' "non-canonical" row number formats for catalogs.
#' 
#' @param c1 A vector of character strings with the first 3 characters
#' being the source trinucleotide and the last character being the
#' mutated (center) nucleotide. E.g. \code{ACTA} means a mutation from
#' \code{ACT > AAT}.
#' 
#' @keywords internal

Restaple96 <- function(c1) {
  retval <-
    paste0(substr(c1, 1, 1),
           "[",
           substr(c1, 2, 2),
           ">",
           substr(c1, 4, 4),
           "]",
           substr(c1, 3, 3))
  return(retval)
}


#' Read a 96-channel spectra (or signature) catalog where rownames are e.g. "A[C>A]T".
#' 
#' The file needs to have the rownames in the first column.
#' @keywords internal
#' @importFrom utils read.table

ReadStapleGT96SBS <- function(file,
                              ref.genome   = NULL,
                              region       = "unknown",
                              catalog.type = "counts",
                              abundance    = NULL,
                              sep          = "\t") {
  df <- read.table(file,
                   stringsAsFactors = FALSE,
                   as.is            = TRUE,
                   header           = TRUE, 
                   check.names      = FALSE,
                   sep              = sep)
  # Careful, df has 96 rows
  stopifnot(nrow(df)==96)
  c1 <- df[ , 1]

  # E.g. "A[C>A]T" -> "ACTA"
  rn <- Unstaple96(c1)
  rownames(df) <- rn
  df <- df[ , -1]
  df.cat <- as.catalog(df, ref.genome = ref.genome, region = region,
                       catalog.type = catalog.type, abundance = abundance)
  return(df.cat)
}



# 1:Del:C:0 -> DEL:C:1:0
# 
# SigProID2ICAMS("DEL:C:1:0")

TransRownames.ID.SigPro.ICAMS <- function(vector.of.rownames) {
  retval <- 
    matrix(unlist(strsplit(vector.of.rownames, ":")), nrow = 4)
  rownames(retval) <- c("Indel_size", "Type", "Subtype", "Repeat_MH_size")
  # NOTE: MH_size is the size of the repeat unit OR microhomology (MH)
  ty <- retval["Type", , drop = FALSE]
  ty <- toupper(ty)
  ty <- sub("R", "repeats", ty)
  ty <- sub("M", "MH", ty)
  # Not done
  return(retval)
}

TransRownames.ID.ICAMS.SigPro <- function(vector.of.rownames) {
  retval <- 
    matrix(unlist(strsplit(vector.of.rownames, ":")), nrow = 4)
  # Not done
  return(retval)
  
}
