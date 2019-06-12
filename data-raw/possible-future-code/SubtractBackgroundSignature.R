#' @keywords internal
SubtractBackground1 <- function(spectrum, background) {
     
  
  
}



#' Subtract a background signature of unknown intensity from mutational specta
#' 
#' @param spectra An \code{ICAMS}
#'  \code{counts} or \code{density} \code{catalog}.
#'
#' @param background An \code{ICAMS} 
#'  \code{counts.signature} or \code{density.signature} \code{catalog}. 
#'  
#'  \itemize{
#'  
#'  \item If
#'  the \code{catalog.type} of
#'  \code{signature} is \code{density.signature} then the \code{catalog.type}
#'  of \code{spectra} must be \code{density}.
#'  
#'  \item If the \code{catalog.type} of
#'  \code{signature} is \code{counts.signature} then the \code{catalog.type}
#'  of \code{spectra} must be \code{counts} \strong{and} the
#'  \code{abundance} attributes must match.
#'  
#'  } 
#'  
#' @export
SubtractBackgroundSignature <- function(spectra, signature) {
    out.spectra <- 
      apply(X = colnames(spectra),
            MARGIN = 2,
            function(spectrum, signature) {
              SubtractBackground1(spectrum, signature)
            }
      )
    return(out.spectra) # Not sure how to make this a catalog again
}

