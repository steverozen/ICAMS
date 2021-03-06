% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ID_functions.R
\name{FindMaxRepeatIns}
\alias{FindMaxRepeatIns}
\title{Return the number of repeat units in which an insertion
is embedded.}
\usage{
FindMaxRepeatIns(context, rep.unit.seq, pos)
}
\arguments{
\item{context}{A string into which \code{rep.unit.seq} was
inserted at position \code{pos}.}

\item{rep.unit.seq}{The inserted sequence and candidate repeat unit
sequence.}

\item{pos}{\code{rep.unit.seq} is understood to be inserted between
positions \code{pos} and \code{pos + 1}.}
}
\value{
If same sequence as \code{rep.unit.seq} occurs ending at
  \code{pos} or starting at \code{pos + 1} then the number of
  repeat units before the insertion, otherwise 0.
}
\description{
Return the number of repeat units in which an insertion
is embedded.
}
\details{
For example

\preformatted{

rep.unit.seq = ac
pos = 2
context = xyaczt
return 1

rep.unit.seq = ac
pos = 4
context = xyaczt
return 1

rep.unit.seq = cgct
pos = 2
rep.unit.seq = at
return 0

context = gacacacacg
rep.unit.seq = ac
pos = any of 1, 3, 5, 7, 9
return 4
}

If 
\code{substr(context, pos, pos + nchar(rep.unit.seq) - 1) != rep.unit.seq},
then stop.
}
\keyword{internal}
