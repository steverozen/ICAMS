#' Infer the class of catalog in a file.
#'
#' @param file Path to a catalog on disk in the standardized format.
#'
#' @return A character string with class appropriate for the catalog on disk.
#'
#' @keywords internal
InferClassOfCatalogForRead <- function(file) {
  cos <- data.table::fread(file) # a one-column data.table never drops to a vector
  # StopIfNrowIllegal(cos)
  catalog.info <- InferCatalogInfo(cos) 
  return(catalog.info)
}