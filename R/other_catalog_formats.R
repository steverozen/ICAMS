#' Read a 192-channel spectra (or signature) catalog in Duke-NUS format
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

#' Convert SBS96-channel mutations-type identifiers like this \code{"A[C>A]T" -> "ACTA"}
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

#' Convert SBS1536-channel mutations-type identifiers like this \code{"AC[C>A]GT" -> "ACCGTA"}
#'
#' @inheritParams Unstaple96
#'
#' @keywords internal
Unstaple1536 <- function(c1) {
  retval <-
    paste0(substr(c1, 1, 1),
           substr(c1, 2, 2),
           substr(c1, 4, 4),
           substr(c1, 8, 8),
           substr(c1, 9, 9),
           substr(c1, 6, 6))
  return(retval)
}


#' Convert DBS78-channel mutations-type identifiers like this \code{"AC>GA" -> "ACGA"}
#'
#' @param c1 A vector of character strings with a \code{>} sign separating
#' reference and variant context e.g. \code{AC>GA}.
#'
#' @keywords internal
Unstaple78 <- function(c1) {
  retval <-
    paste0(substr(c1, 1, 1),
           substr(c1, 2, 2),
           substr(c1, 4, 4),
           substr(c1, 5, 5))
  return(retval)
}


#' Convert 96-channel mutation-type identifiers like this \code{"ACTA" -> "A[C>A]T"}
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

#' Convert 1536-channel mutation-type identifiers like this \code{"ACCGTA" -> "AC[C>A]GT"}
#'
#' This is an internal function needed for generating
#' "non-canonical" row number formats for catalogs.
#'
#' @param c1 A vector of character strings with the first 5 characters
#' being the source trinucleotide and the last character being the
#' mutated (center) nucleotide. E.g. \code{ACCGTA} means a mutation from
#' \code{ACCGT > ACAGT}.
#'
#' @keywords internal

Restaple1536 <- function(c1) {
  retval <-
    paste0(substr(c1, 1, 1),
           substr(c1, 2, 2),
           "[",
           substr(c1, 3, 3),
           ">",
           substr(c1, 6, 6),
           "]",
           substr(c1, 4, 4),
           substr(c1, 5, 5))
  return(retval)
}



#' Read a 96-channel spectra (or signature) catalog where rownames are e.g. "A[C>A]T"
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
  stopifnot(nrow(df) == 96)
  
  if (ncol(df) == 1) {
    stop("Only one column after reading ", file, 
         "; do you need to set sep to something other than \"", sep, "\"?")
    
  }
  
  c1 <- df[ , 1]

  # E.g. "A[C>A]T" -> "ACTA"
  rn <- Unstaple96(c1)
  rownames(df) <- rn
  df <- df[ , -1]
  df.cat <- as.catalog(df, ref.genome = ref.genome, region = region,
                       catalog.type = catalog.type, abundance = abundance)
  return(df.cat)
}


#' For indels, convert SigProfiler rownames into ICAMS/PCAWG7 rownames
#'
#' @examples 
#' ICAMS:::TransRownames.ID.SigPro.PCAWG("1:Del:C:0") # DEL:C:1:0;
#' ICAMS:::TransRownames.ID.SigPro.PCAWG("2:Ins:R:5") # INS:repeat:2:5+
#'
#' @keywords  internal
#' 
TransRownames.ID.SigPro.PCAWG <- function(vector.of.rownames) {
  retval <-
    matrix(unlist(strsplit(vector.of.rownames, ":")), nrow = 4)
  rownames(retval) <- c("Indel_size", "Type", "Subtype", "Repeat_MH_size")
  # NOTE: MH_size is the size of the repeat unit OR microhomology (MH)
  
  ty <- retval["Type", ]
  ty <- toupper(ty)
  retval["Type", ] <- ty

  subty <- retval["Subtype", ]
  subty <- sub("R", "repeats", subty)
  subty <- sub("M", "MH", subty)
  retval["Subtype", ] <- subty

  indsi <- retval["Indel_size", ]
  indsi <- sub("5", "5\\+", indsi)
  retval["Indel_size", ] <- indsi

  repmh <- retval["Repeat_MH_size", ]
  repmh <- sub("5", "5\\+", repmh)
  retval["Repeat_MH_size", ] <- repmh
  
  retval <- paste(ty,subty,indsi,repmh,sep = ":")
  #  Each name in retval has the order of an ICAMS/PCAWG7 name.
  
  return(retval)
}


#' For indels, convert ICAMS/PCAWG7 rownames into SigProfiler rownames
#'
#' @examples
#' ICAMS:::TransRownames.ID.PCAWG.SigPro("DEL:C:1:0") # 1:Del:C:0;
#' ICAMS:::TransRownames.ID.PCAWG.SigPro("INS:repeat:2:5+") # 2:Ins:R:5
#'
#' @keywords  internal
TransRownames.ID.PCAWG.SigPro <- function(vector.of.rownames) {
  retval <-
    matrix(unlist(strsplit(vector.of.rownames, ":")), nrow = 4)
  
  rownames(retval) <- c("Type", "Subtype", "Indel_size", "Repeat_MH_size")
  
  ## Not done
  ty <- retval["Type", ]
  ty <- sub("INS", "Ins", ty)
  ty <- sub("DEL", "Del", ty)
  retval["Type", ] <- ty
  
  subty <- retval["Subtype", ]
  subty <- sub("repeats", "R", subty)
  subty <- sub("MH", "M", subty)
  retval["Subtype", ] <- subty
  
  indsi <- retval["Indel_size", ]
  indsi <- sub("5\\+", "5", indsi)
  retval["Indel_size", ] <- indsi
  
  repmh <- retval["Repeat_MH_size", ]
  repmh <- sub("5\\+", "5", repmh)
  retval["Repeat_MH_size", ] <- repmh
  
  retval <- paste(indsi,ty,subty,repmh,sep = ":")
  # Each name in retval has the order of a SigProfiler name.
  
  return(retval)
}


#' Write Indel Catalogs in SigProExtractor format
#' 
#' Write Indel Catalogs in SigProExtractor format to a file.
#' 
#' @param catalog A catalog as defined in \code{\link{ICAMS}};
#' see also \code{\link{as.catalog}}.
#'
#' @param file The path to the file to be created.
#'
#' @param strict If TRUE, do additional checks on the input, and stop if the
#' checks fail.
#' 
#' @param sep Separator to use in the output file. In older version 
#'  SigProfiler read comma-separated files; as of May 2020 it
#'  reads tab-separated files.
#'
#' @note In ID (small insertion and deletion) catalogs in SigProExtractor format, 
#' deletion repeat sizes range from 0 to 5, rather than 0 to 5+.
#'
#' @keywords internal
WriteCatalogIndelSigPro <- function(catalog, file, strict = TRUE, sep = "\t"){
  mut.categories <- rownames(catalog)
  stopifnot(nrow(catalog) == 83)
  if (strict) {
    stopifnot(mut.categories == ICAMS::catalog.row.order$ID)
  }
  ## Change row headers to SigProExtractor-format
  mut.categories <- TransRownames.ID.PCAWG.SigPro(mut.categories)
  rownames(catalog) <- mut.categories
  ## Change row order to SigProExtractor-order  
  catalog <- catalog[ICAMS.to.SigPro.ID[,1], , drop = FALSE]
  DT <- data.table("Mutation Types" = rownames(catalog),catalog)
  ## Write the SigProExtractor-formatted catalog into a csv file
  fwrite(DT, file = file, sep = sep)
}

## Linker from PCAWG(ICAMS)-formatted to SigProExtractor-formatted indel names
##
## This data is designed for converting ICAMS-formatted indel names to
## SigProExtractor-formatted indel names.
##
## @format A 83*1 matrix. Its contents (first column) contain SigProExtractor
## formatted indel names in SigProExtractor order. Its rownames refer to the 
## corresponding PCAWG(ICAMS)-formatted indel names.
##
## @name ICAMS.to.SigPro.ID
## 
## @examples 
## ICAMS.to.SigPro.ID
## #              SigPro.ID.names      
## # DEL:C:1:0    "1:Del:C:0"    
## # DEL:C:1:1    "1:Del:C:1"
## # DEL:C:1:2    "1:Del:C:2"  
## # DEL:C:1:3    "1:Del:C:3" 
## # DEL:C:1:4    "1:Del:C:4"  
## #       ...           ...      
"ICAMS.to.SigPro.ID"


## Linker from SigProExtractor-formatted to PCAWG(ICAMS)-formatted indel names
##
## This data is designed for converting SigProExtractor-formatted indel names
## to ICAMS-formatted indel names.
##
## @format A 83*1 matrix. Its contents (first column) contain PCAWG(ICAMS)-
## formatted indel names in PCAWG(ICAMS) order. Its rownames refer to the 
## corresponding SigProExtractor indel names.
##
## @name SigPro.to.ICAMS.ID
## 
## @examples 
## SigPro.to.ICAMS.ID
## #              ICAMS.ID.names      
## # 1:Del:C:0    "DEL:C:1:0"    
## # 1:Del:C:1    "DEL:C:1:1"
## # 1:Del:C:2    "DEL:C:1:2"  
## # 1:Del:C:3    "DEL:C:1:3" 
## # 1:Del:C:4    "DEL:C:1:4"  
## #       ...           ...  
"SigPro.to.ICAMS.ID"

#' Covert an ICAMS SBS96 Catalog to SigProfiler format
#' 
#' @param input.catalog Either a character string, in which case this is the
#'   path to a file containing a catalog in \code{\link[ICAMS]{ICAMS}}
#'   format, or an in-memory \code{\link[ICAMS]{ICAMS}} catalog.
#'   
#' @param file The path of the file to be written.
#' 
#' @param sep Separator to use in the output file. 
#'
#' @importFrom data.table as.data.table
#' 
#' @importFrom utils write.table 
#' 
#' @keywords internal
ConvertICAMSCatalogToSigProSBS96 <- function(input.catalog, file, sep = "\t") {
  if (inherits(input.catalog, "character")) {
    input.catalog <- ICAMS::ReadCatalog(input.catalog)
  } 
  new.list <- lapply(row.names(input.catalog), function(x){
    new <- paste(substring(x, 1, 1), "[",
                 substring(x, 2, 2), ">",
                 substring(x, 4, 4), "]",
                 substring(x, 3, 3), sep = "")
    
  })
  
  row.names(input.catalog) <- unlist(new.list)
  input.catalog <- 
    input.catalog[catalog.row.headers.sp$SBS96, , drop = FALSE]
  DT <- as.data.table(input.catalog)
  input.sigpro <- cbind("MutationType" = catalog.row.headers.sp$SBS96, DT)
  write.table(input.sigpro, file, sep = sep, col.names = TRUE,
              row.names = FALSE, quote = FALSE)
}

#' Covert an ICAMS Catalog to SigProfiler format
#' 
#' @param input.catalog Either a character string, in which case this is the
#'   path to a file containing a catalog in \code{\link[ICAMS]{ICAMS}}
#'   format, or an in-memory \code{\link[ICAMS]{ICAMS}} catalog.
#'   
#' @param file The path of the file to be written.
#' 
#' @param sep Separator to use in the output file. 
#'
#' @importFrom data.table as.data.table
#' 
#' @importFrom utils write.table 
#' 
#' @note This function can only transform SBS96, DBS78 and ID ICAMS catalog
#' to SigProfiler format.
#' 
#' @keywords internal
ConvertCatalogToSigProfilerFormat <- function(input.catalog, file, sep = "\t") {
  if (inherits(input.catalog, "character")) {
    input.catalog <- ICAMS::ReadCatalog(input.catalog)
  } 
  if (nrow(input.catalog) == 96) {
    mutation.type <- "SBS96"
  } else if (nrow(input.catalog) == 78) {
    mutation.type <- "DBS78"
  } else if (nrow(input.catalog) == 83) {
    mutation.type <- "ID83"
  } else {
    stop("Can only convert SBS96, DBS78 and ID ICAMS catalog to SigProfiler format")
  }
  
  if (mutation.type == "SBS96") {
    new.list <- lapply(row.names(input.catalog), function(x){
      new <- paste(substring(x, 1, 1), "[",
                   substring(x, 2, 2), ">",
                   substring(x, 4, 4), "]",
                   substring(x, 3, 3), sep = "")
    })
    row.names(input.catalog) <- unlist(new.list)
  } else if (mutation.type == "DBS78") {
    new.list <- lapply(row.names(input.catalog), function(x){
      new <- paste(substring(x, 1, 2), ">",
                   substring(x, 3, 4), sep = "")
    })
    row.names(input.catalog) <- unlist(new.list)
  } else if (mutation.type == "ID83") {
    row.names(input.catalog) <- 
      TransRownames.ID.PCAWG.SigPro(row.names(input.catalog))
  }
  
  input.catalog <- 
    input.catalog[catalog.row.headers.sp[[mutation.type]], , drop = FALSE]
  DT <- as.data.table(input.catalog)
  input.sigpro <- 
    cbind("MutationType" = catalog.row.headers.sp[[mutation.type]], DT)
  write.table(input.sigpro, file, sep = sep, col.names = TRUE,
              row.names = FALSE, quote = FALSE)
}
