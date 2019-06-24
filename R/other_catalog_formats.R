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
                   stringsAsFactors = F,
                   as.is=T,
                   header=T, 
                   check.names = F)
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

