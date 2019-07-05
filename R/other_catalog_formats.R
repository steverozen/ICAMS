#' Read a 192-channel spectra (or signature) catalog in Duke-NUS format.
#' 
#' The file needs to have the column names Before	Ref	After	Var
#' in the first 4 columns 
#' @keywords internal
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
  stopifnot(nrow(df)==96)
  c1 <- df[ , 1]
  
  # E.g. "A[C>A]T" -> "ACTA"
  rn <-
    paste0(substr(c1, 1, 1),
           substr(c1, 3, 3),
           substr(c1, 7, 7),
           substr(c1, 5, 5))
  
  rownames(df) <- rn
  df <- df[ , -1]
  df.cat <- as.catalog(df, ref.genome = ref.genome, region = region,
                       catalog.type = catalog.type, abundance = abundance)
  return(df.cat)
}

