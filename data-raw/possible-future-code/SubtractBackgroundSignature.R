#' @keywords internal
LikelihoodBackgroundOnly <- function(spectra, background.sig) {
     # Find max likelihood and count to generate each
     # spectrum in spectra using only background (the model could contain 
     # an expected distribution of intensity of background mutagenesis)
  
     # will compare this with a model with another signature
     # (the same across all samples) with variable count
     # per sample
  
     # To find the background signature, find the shape and
     # and intensity that gives the observed backgrounds.
  
}



#' Subtract background if required by likelihood ratio test for presenece of non-backgound signature
#' 
#' @param spectra An \code{ICAMS}
#'  \code{counts} or \code{density} \code{catalog}.
#'
#' @param background.sig An \code{ICAMS} 
#'  \code{counts.signature} or \code{density.signature} \code{catalog}. 
#'  
#'  \itemize{
#'  
#'  \item If
#'  the class and \code{background.sig} and \code{catalog}
#'  must be the same, and
#'  
#'  \item if \code{catalog.type} of
#'  \code{background.sig} is \code{density.signature} then the \code{catalog.type}
#'  of \code{spectra} must be \code{density}, and
#'  
#'  \item If the \code{catalog.type} of
#'  \code{background.sig} is \code{counts.signature} then the \code{catalog.type}
#'  of \code{spectra} must be \code{counts} \strong{and} the
#'  \code{abundance} attributes must match.
#'  
#'  } 
#'  
#' @export
LSubtractBackground <- function(spectra, background.sig) {
  out.spectra <- 
    apply(X = colnames(spectra),
          MARGIN = 2,
          function(spectrum) {
            SubtractBackground1(spectrum, background.sig, background.count)
          }
    )
  return(out.spectra) # Not sure how to make this a catalog again
}




#' Subtract background if required by a mutation count test methods
#' 
#' @inheritParams LSubtractBackground
#'
#' @param background.count Expected total number of mutatations due to 
#' \code{background.sig}.
#'  
#' @export
CSubtractBackground <- function(spectra, background.sig, background.count) {
    out.spectra <- 
      apply(X = colnames(spectra),
            MARGIN = 2,
            function(spectrum) {
              SubtractBackground1(spectrum, background.sig, background.count)
            }
      )
    return(out.spectra) # Not sure how to make this a catalog again
}

