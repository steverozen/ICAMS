% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/other_catalog_formats.R
\name{TransRownames.ID.SigPro.PCAWG}
\alias{TransRownames.ID.SigPro.PCAWG}
\title{For indels, convert SigProfiler rownames into ICAMS/PCAWG7 rownames}
\usage{
TransRownames.ID.SigPro.PCAWG(vector.of.rownames)
}
\description{
For indels, convert SigProfiler rownames into ICAMS/PCAWG7 rownames
}
\examples{
ICAMS:::TransRownames.ID.SigPro.PCAWG("1:Del:C:0") # DEL:C:1:0;
ICAMS:::TransRownames.ID.SigPro.PCAWG("2:Ins:R:5") # INS:repeat:2:5+

}
\keyword{internal}
