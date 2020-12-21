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
