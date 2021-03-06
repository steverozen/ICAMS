% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VCF_to_catalog_functions.R
\name{CheckSBSClassInVCF}
\alias{CheckSBSClassInVCF}
\title{Check SBS mutation class in VCF with the corresponding SBS mutation matrix}
\usage{
CheckSBSClassInVCF(vcf, mat, sample.id)
}
\arguments{
\item{vcf}{An annotated SBS VCF with columns of SBS mutation
classes added by \code{AddSBSClass}.}

\item{mat}{The SBS mutation count matrix.}

\item{sample.id}{Usually the sample id, but defaults to "count".}
}
\description{
Check SBS mutation class in VCF with the corresponding SBS mutation matrix
}
\keyword{internal}
