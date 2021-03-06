% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VCF_to_catalog_functions.R
\name{AddTranscript}
\alias{AddTranscript}
\title{Add transcript information to a data frame with mutation records}
\usage{
AddTranscript(df, trans.ranges = NULL, ref.genome, name.of.VCF = NULL)
}
\arguments{
\item{df}{A data frame storing mutation records of a VCF file.}

\item{trans.ranges}{A \code{\link[data.table]{data.table}} which contains
transcript range and strand information. Please refer to
\code{\link{TranscriptRanges}} for more details.}

\item{ref.genome}{A \code{ref.genome} argument as described in
\code{\link{ICAMS}}.}

\item{name.of.VCF}{Name of the VCF file.}
}
\value{
A data frame with new columns added to the input data frame,
    which contain the mutated gene's name, range and strand information.
}
\description{
Add transcript information to a data frame with mutation records
}
\keyword{internal}
