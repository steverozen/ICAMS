#' This function converts an data.table imported
#' from external catalog text file into ICAMS
#' internal catalog object of appropriate type.
#'
#' 
#' @keywords internal
InferCatalogInfo <- function(object) {
  
  nrow <- nrow(object)

  if (nrow == 96)  {
    return(Make96Catalog(object))
  } 
  if (nrow == 192)  {
    return(MakeSBS192Catalog(object))
  }
  if (nrow == 1536) {
    return(MakeSBS1536Catalog(object))
  }
  if (nrow == 78)   {
    return(MakeDBS78Catalog(object))
  }
  if (nrow == 144){
    return(MakeDBS144Catalog(object))
  }
  if (nrow == 136){
    return(MakeDBS136Catalog(object))
  }
  if (nrow == 83)  {
    return(MakeID83Catalog(object))
  }
  if (nrow == 166)  {
    return(MakeID166Catalog(object))
  }
  if (nrow == 1697){ # TODO(Wuyang)
    return(MakeCOMPOSITECatalog(object))
  }  
  stop("\nThe number of rows in the input object must be one of\n",
       "96 192 1536 78 144 136 83 166 1697\ngot ", nrow)
  
}

#' These two functions is applicable only for 
#' internal ICAMS-formatted catalog object.
#'
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
  if (nrow == 166)   return("ID166")
  if (nrow == 1697) return("COMPOSITE")
  
  stop("\nThe number of rows in the input object must be one of\n",
       "96 192 1536 78 144 136 83 166 1697\ngot ", nrow)
  
}


InferCatalogClassString <- function(object) {
  
  prefix <- InferCatalogClassPrefix(object)
  if (prefix == "ID") prefix <- "Indel"
  return(paste0(prefix, "Catalog"))
}


#' Check whether an R object contains one of the ICAMS catalog classes 
#'
#' @param object An R object.
#'
#' @return A logical value.
#' 
#' @export
#'
#' @examples
#' # Create a matrix with all values being 1
#' object <- matrix(1, nrow = 96, ncol = 1, 
#'                  dimnames = list(catalog.row.order$SBS96))
#' IsICAMSCatalog(object) # FALSE
#' 
#' # Use as.catalog to add class attribute to object
#' catalog <- as.catalog(object)
#' IsICAMSCatalog(catalog) # TRUE      
IsICAMSCatalog <- function(object) {
  supported.classes <-
    paste0(
    c(paste0("SBS", c(96, 192, 1536)),
      paste0("DBS", c(78, 144, 136)),
      "Indel", "ID166", "COMPOSITE"),
    "Catalog")
  return(inherits(x = object, what = supported.classes))
}


## Convert external catalog files with 96 rows into ICAMS internal catalog format.
#' @keywords internal
Make96Catalog <- function(object) {
  # ICAMS and COSMIC SBS96 csv format
  if("C>A" %in% unlist(object[,1]) && "TTT" %in% unlist(object[,2])) {
    # Convert ICAMS / COSMIC csv into ICAMS internal format
    return(MakeSBS96CatalogFromICAMSExt(object))
  }
  ## SigPro SBS96 txt format
  if("A[C>A]A" %in% unlist(object[,1])) {
    return(MakeSBS96CatalogFromSigPro(object))
  }
  ## COSMIC / SigPro ID96 text format
  if(all(c("1:Del:C:0","1:Ins:C:0","1:Del:T:0","1:Ins:T:0",
           "2:Del:R:0","5:Ins:R:0","5:Del:M:5","5:Ins:M:5",
           "complex","non_matching") %in% unlist(object[,1]))) {
    return(MakeID83ICAMSFromSigProID96(object))
  }
  stop("96 mutation types, but not an SBS96 catalog in ICAMS/COSMIC/SigProfiler",
       " or ID96 catalog in ICAMS/COSMIC format")
}
## Both COSMIC and ICAMS use ICAMS-external SBS96 format: 
## C/T>D/V NYN
MakeSBS96CatalogFromICAMSExt <- function(cos) {
  ref.gt.var       <- unlist(cos[, 1])
  before.ref.after <- unlist(cos[, 2])
  var <- substring(ref.gt.var, 3, 3)
  out <- cos[, -(1 : 2), drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- paste0(before.ref.after, var)
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  
  if(!setequal(rownames(out),ICAMS::catalog.row.order$SBS96)){
    stop("The mutation types in this SBS96 catalog is not in correct ICAMS format",
      " check ICAMS::catalog.row.order$SBS96 for more details.")
  }
  out <- out[ICAMS::catalog.row.order$SBS96, , drop = FALSE]
  out <- as.matrix(out)
  class(out) <- c("SBS96Catalog", class(out))
  return(out)
}
## N[C/T>D/V]N
MakeSBS96CatalogFromSigPro <- function(cos) {
  rownames <- Unstaple96(unlist(cos[ , 1]))
  cos <- cos[ , -1]
  out <- as.matrix(cos)
  rownames(out) <- rownames
  if(!setequal(rownames(out),ICAMS::catalog.row.order$SBS96)){
    stop("The mutation types in this SBS96 catalog is not in correct ICAMS format",
      " check ICAMS::catalog.row.order$SBS96 for more details.")
  }
  out <- out[ICAMS::catalog.row.order$SBS96, , drop = FALSE]
  class(out) <- c("SBS96Catalog", class(out))
  return(out)
}
## See ICAMS:::TransRownames.ID.SigPro.PCAWG() for more details
## on ID96/ID83 format.
MakeID83ICAMSFromSigProID96 <- function(cos) {

  ## Row names before conversion.
  headerPrev <- unlist(cos[,1])
  names(headerPrev) <- NULL

  ## Row headers not in ID83
  headerIDsToBeDeleted <- which(headerPrev %in% c("complex","non_matching","non-matching", "non.matching"))
  headerIDsToBeDeleted <- c(headerIDsToBeDeleted, grep("Ins:M",headerPrev))
  cos <- cos[setdiff(1:96,headerIDsToBeDeleted),]
  headers <- unlist(cos[,1])
  headers <- TransRownames.ID.SigPro.PCAWG(headers)

  out <- cos[, -1] # No need to add drop = FALSE
  out <- as.matrix(out)
  rownames(out) <- headers
  
  if(!setequal(rownames(out),ICAMS::catalog.row.order$ID)){
    stop("The mutation types in this ID83 catalog is not in correct ICAMS format",
      " check ICAMS::catalog.row.order$ID for more details.")
  }
  
  # Need drop = FALSE, otherwise out will not be a matrix and has attributes
  # other than names. Then in ReadCatalogInternal(), as.catalog() will throw an
  # error as is.vector(out) is FALSE.
  out <- out[ICAMS::catalog.row.order$ID, , drop = FALSE]
  ## TODO(Wuyang): check if the matrix is unified.
  class(out) <- c("IndelCatalog", class(out))
  return(out)
}

## Convert external catalog files with 192 rows (should always be SBS192)
## into ICAMS SBS192 internal catalog object.
MakeSBS192Catalog <- function(object) {
  # ICAMS / COSMIC SBS192 csv format
  if(setequal(c("T","U"), unique(unlist(object[,1]))) && "C>A" %in% unlist(object[,2]) && "TTT" %in% unlist(object[,3])) {
    # Convert ICAMS / COSMIC csv into ICAMS internal format
    return(MakeSBS192CatalogFromICAMSExt(object))
  }
  stop("192 mutation types, but not an SBS192 catalog in",
       " ICAMS or COSMIC format")
}
## COSMIC, ICAMS and SigPro all use ICAMS-external SBS192 format:
## T/U C/T>D/V NYN
MakeSBS192CatalogFromICAMSExt <- function(cos) {
  
  stopifnot(nrow(cos) == 192)
  ref.gt.var       <- unlist(cos[, 2])
  before.ref.after <- unlist(cos[, 3])

  ## Find the rows labeled with "T", indicating the
  ## SBS is on the transcribed (which is the *antisense*) strand.
  transcribed.strand.pos <- which(cos[, 1] == 'T')

  before.ref.after[transcribed.strand.pos] <-
    revc(before.ref.after[transcribed.strand.pos])

  var <- substring(ref.gt.var, 3, 3)
  var[transcribed.strand.pos] <- revc(var[transcribed.strand.pos])

  internalHeaders <- paste0(before.ref.after, var)
  stopifnot(setequal(internalHeaders, ICAMS::catalog.row.order$SBS192))

  out <- cos[, -(1 : 3), drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- internalHeaders
  if(!setequal(rownames(out),ICAMS::catalog.row.order$SBS192)){
    stop("The mutation types in this SBS192 catalog is not in correct ICAMS format",
      " check ICAMS::catalog.row.order$SBS192 for more details.")
  }
  out <- out[ICAMS::catalog.row.order$SBS192, , drop = FALSE]
  class(out) <- c("SBS192Catalog", class(out))
  return(out)
}


## Convert external catalog files with 1536 rows (should always be SBS1536)
## into ICAMS SBS1536 internal catalog object.
MakeSBS1536Catalog <- function(object) {
  # ICAMS SBS1536 csv format
  if("C>A" %in% unlist(object[,1]) && "TTTGT" %in% unlist(object[,2])) {
    # Convert ICAMS / COSMIC csv into ICAMS internal format
    return(MakeSBS1536CatalogFromICAMSExt(object))
  }
  stop("1536 mutation types, but not an SBS1536 catalog in",
       " ICAMS format")
}
## ICAMS uses ICAMS-external SBS1536 format:
## C/T>D/V NNYNN
MakeSBS1536CatalogFromICAMSExt <- function(cos) {
  
  stopifnot(nrow(cos) == 1536)
  names(cos)[1:2] <- c("Mutation type", "Pentanucleotide")
  ref.gt.var       <- cos[["Mutation type"]]
  before.ref.after <- cos[["Pentanucleotide"]]
  var <- substring(ref.gt.var, 3, 3)
  out <- as.matrix(cos[ , -(1 : 2)], drop = FALSE)
  rownames(out) <- paste0(before.ref.after, var)
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  # Check whether the row names are setequal to standard internal row names.
  if(!setequal(rownames(out),ICAMS::catalog.row.order$SBS1536)){
    stop("The mutation types in this SBS1536 catalog is not in correct ICAMS format",
      " check ICAMS::catalog.row.order$SBS1536 for more details.")
  }
  out <- out[ICAMS::catalog.row.order$SBS1536, , drop = FALSE]
  return(out)
}


## Convert external catalog files with 78 rows into ICAMS DBS78 internal catalog format.
MakeDBS78Catalog <- function(object) {
  # ICAMS DBS78 csv format
  if("AC" %in% unlist(object[,1]) && "GA" %in% unlist(object[,2])) {
    # Convert ICAMS-external csv into ICAMS internal format
    return(MakeDBS78CatalogFromICAMSExt(object))
  }
  ## SigPro DBS78 txt format
  if("AC>GA" %in% unlist(object[,1])) {
    return(MakeDBS78CatalogFromSigPro(object))
  }
  stop("78 mutation types, but not a DBS78 catalog in",
       " ICAMS, COSMIC, or SigProfiler format")
}
## NN(Ref) NN(Var)
MakeDBS78CatalogFromICAMSExt <- function(cos) {
  ref <- unlist(cos[, 1])
  var <- unlist(cos[, 2])
  out <- cos[, -(1 : 2), drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- paste0(ref, var)
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  
  if(!setequal(rownames(out),ICAMS::catalog.row.order$DBS78)){
    stop("The mutation types in this DBS78 catalog is not in correct ICAMS format",
      " check ICAMS::catalog.row.order$DBS78 for more details.")
  }
  out <- out[ICAMS::catalog.row.order$DBS78, , drop = FALSE]
  out <- as.matrix(out)
  class(out) <- c("DBS78Catalog", class(out))
  return(out)
}
## Both COSMIC and SigPro use SigPro DBS78 format:
## NN>NN
MakeDBS78CatalogFromSigPro <- function(cos) {
  rownames <- Unstaple78(unlist(cos[ , 1]))
  cos <- cos[ , -1]
  out <- as.matrix(cos)
  rownames(out) <- rownames
  
  if(!setequal(rownames(out),ICAMS::catalog.row.order$DBS78)){
    stop("The mutation types in this DBS78 catalog is not in correct ICAMS format",
      " check ICAMS::catalog.row.order$DBS78 for more details.")
  }
  out <- out[ICAMS::catalog.row.order$DBS78, , drop = FALSE]
  class(out) <- c("DBS78Catalog", class(out))
  return(out)
}

## Convert external catalog files with 136 rows into ICAMS DBS136 internal catalog format.
MakeDBS136Catalog <- function(object) {
  # ICAMS DBS136 csv format
  if("AATA" %in% unlist(object[,1])) {
    # Convert ICAMS-external csv into ICAMS internal format
    return(MakeDBS136CatalogFromICAMSExt(object))
  }
  stop("136 mutation types, but not a DBS136 catalog in",
       " ICAMS format")
}
## Quad
## NNNN
MakeDBS136CatalogFromICAMSExt <- function(cos) {

  names(cos)[1] <- "Quad"
  out <- cos[, -1, drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- cos$Quad

  if(!setequal(rownames(out), ICAMS::catalog.row.order$DBS136)){
    stop("The mutation types in this DBS136 catalog is not in correct ICAMS format",
      " check ICAMS::catalog.row.order$DBS136 for more details.")
  }
  out <- out[ICAMS::catalog.row.order$DBS136, , drop = FALSE]
  return(out)
}


## Convert external catalog files with 144 rows into ICAMS DBS144 internal catalog format.
MakeDBS144Catalog <- function(object) {
  # ICAMS DBS144 csv format
  if("AC" %in% unlist(object[,1]) && "GA" %in% unlist(object[,2])) {
    # Convert ICAMS-external csv into ICAMS internal format
    return(MakeDBS144CatalogFromICAMSExt(object))
  }
  stop("144 mutation types, but not a DBS144 catalog in",
       " ICAMS format")
}
## Ref Var
## NN > NN 
MakeDBS144CatalogFromICAMSExt <- function(cos) {
  
  stopifnot(nrow(cos) == 144)
  names(cos)[1 : 2] <- c("Ref", "Var")
  out <- cos[, -(1 : 2), drop = FALSE]
  out <- as.matrix(out)
  rn <- paste0(cos$Ref, cos$Var)
  rownames(out) <- rn
  if(!setequal(rownames(out), ICAMS::catalog.row.order$DBS144)){
    stop("The mutation types in this DBS144 catalog is not in correct ICAMS format",
      " check ICAMS::catalog.row.order$DBS144 for more details.")
  }
  out <- out[ICAMS::catalog.row.order$DBS144, , drop = FALSE]
  return(out)
}


## Convert external catalog files with 83 rows into ICAMS ID83 internal catalog format.
MakeID83Catalog <- function(object) {
  # ICAMS ID83 csv format
  if(setequal(c("DEL","INS"),unique(unlist(object[,1]))) &&
     "repeats" %in% unlist(object[,2]) && 
     "5+" %in% unlist(object[,3]) && 
     "5+" %in% unlist(object[,4])) {
    # Convert ICAMS-external csv into ICAMS internal format
    return(MakeID83CatalogFromICAMSExt(object))
  }
  ## COSMIC and SigPro ID83 txt format
  if("3:Del:R:5" %in% unlist(object[,1])) {
    return(MakeID83CatalogFromSigPro(object))
  }
  stop("83 mutation types, but not a ID83 catalog in",
       " ICAMS, COSMIC, or SigProfiler format")
}

# Convert external catalog files with 166 rows into ICAMS ID166 internal catalog format
# ID166 is genic-intergenic indel catalog, there is an additional first column "Region" added
# to the canonical ID83 catalog CSV file
# "G" stands for "Genic", indicating the indel happens on genic region
# "I" stands for "Intergenic", indicating the indel happens on Intergenic region
MakeID166Catalog <- function(cos) {
  cn <- names(cos)
  ex.cn <- c("Region", "Type", "Subtype", "Indel_size", "Repeat_MH_size")
  names(cos)[1:5] <- ex.cn
  rn <- apply(cos[ , 1:5], MARGIN = 1, paste, collapse = ":")
  out <- as.matrix(cos[ , -(1:5), drop = FALSE])
  
  if (!setequal(rn, ICAMS::catalog.row.order$ID166)) {
    msg <- 
      paste("The row names are not correct:\n",
            "got", paste(rn, collapse = ", "),
            "\nexpected", paste(ICAMS::catalog.row.order$ID166,
                                collapse = ", "))
    stop(msg)
  }
  
  rownames(out) <- rn
  out <- out[ICAMS::catalog.row.order$ID166, , drop = FALSE]
  return(out)
}

## Type Subtype Indel_size Repeat_MH_size
## INS/DEL C/T/repeats/MH 1~5+ 0~5+
MakeID83CatalogFromICAMSExt <- function(cos) {
     
    cn <- names(cos)
    ex.cn <- c("Type", "Subtype", "Indel_size", "Repeat_MH_size")
    # Repeat_MH_size is the size of repeat OR microhomology (MH)
    # if (strict) { for (i in 1:4) { stopifnot(cn[i] == ex.cn[i]) } }
    names(cos)[1:4] <- ex.cn
    rn <- apply(cos[ , 1:4], MARGIN = 1, paste, collapse = ":")
    out <- as.matrix(cos[ , -(1:4), drop = FALSE])
    
    
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
    return(out)
}
## <Indel_size:Type:Subtype:Repeat_MH_size>
## 1~5:Ins/Del:C/T/R/M:0~5
MakeID83CatalogFromSigPro <- function(cos) {
  rn <- unlist(cos[,1])
  if(!setequal(rn, catalog.row.headers.sp$ID83)) {
      msg <- 
        paste("The row names are not correct:\n",
              "got", paste(rn, collapse = ", "),
              "\nexpected", paste(catalog.row.headers.sp$ID83,
                                  collapse = ", "))
      stop(msg)
  }
  internalHeader <- TransRownames.ID.SigPro.PCAWG(rn)
  cos <- as.matrix(cos[,-1])
  rownames(cos) <- internalHeader
  return(cos)
}

## Convert external catalog files with 1697 rows into ICAMS COMPOSITE internal catalog format.
MakeCOMPOSITECatalog <- function(object) {
  # ICAMS COMPOSITE csv format
  if("TCTAAA" %in% unlist(object[,1])) {
    # Convert ICAMS-external csv into ICAMS internal format
    return(MakeCOMPOSITECatalogFromICAMSExt(object))
  }
  stop("1697 mutation types, but not a COMPOSITE catalog in",
       " ICAMS format")
}
## NN(Ref) NN(Var)
MakeCOMPOSITECatalogFromICAMSExt <- function(cos) {
    
  names <- cos[[1]]
  cos <- cos[, -1]
  out <- as.matrix(cos)
  rownames(out) <- names
  
  if(!setequal(rownames(out),ICAMS::catalog.row.order$COMPOSITE)){
    stop("The mutation types in this COMPOSITE catalog is not in correct ICAMS format",
      " check ICAMS::catalog.row.order$COMPOSITE for more details.")
  }
  out <- out[ICAMS::catalog.row.order$COMPOSITE, , drop = FALSE]
  out <- as.matrix(out)
  class(out) <- c("COMPOSITECatalog", class(out))
  return(out)
}




#' Infer the correct rownames for a matrix based on its number of rows
#' 
#' @keywords internal
InferRownames <- function(object) {
  prefix <- InferCatalogClassPrefix(object)
  return(ICAMS::catalog.row.order[[prefix]])
}
