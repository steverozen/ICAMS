% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny_related_functions.R
\name{ReadAndSplitVCFs}
\alias{ReadAndSplitVCFs}
\title{Read and split VCF files}
\usage{
ReadAndSplitVCFs(
  files,
  variant.caller = "unknown",
  num.of.cores = 1,
  names.of.VCFs = NULL,
  tumor.col.names = NA,
  filter.status = NULL,
  get.vaf.function = NULL,
  ...,
  max.vaf.diff = 0.02,
  suppress.discarded.variants.warnings = TRUE,
  always.merge.SBS = FALSE
)
}
\arguments{
\item{files}{Character vector of file paths to the VCF files.}

\item{variant.caller}{Name of the variant caller that produces the VCF, can
be either \code{"strelka"}, \code{"mutect"}, \code{"freebayes"} or
\code{"unknown"}. This information is needed to calculate the VAFs (variant
allele frequencies). If variant caller is \code{"unknown"}(default) and
\code{get.vaf.function} is NULL, then VAF and read depth will be NAs. If
variant caller is \code{"mutect"}, do \strong{not} merge SBSs into DBS.}

\item{num.of.cores}{The number of cores to use. Not available on Windows
unless \code{num.of.cores = 1}.}

\item{names.of.VCFs}{Character vector of names of the VCF files. The order
of names in \code{names.of.VCFs} should match the order of VCF file paths
in \code{files}. If \code{NULL}(default), this function will remove all of
the path up to and including the last path separator (if any) and file
paths without extensions (and the leading dot) will be used as the names of
the VCF files.}

\item{tumor.col.names}{Optional. Only applicable to \strong{Mutect} VCFs.
Vector of column names or column indices in \strong{Mutect} VCFs which
contain the tumor sample information. The order of elements in
\code{tumor.col.names} should match the order of \strong{Mutect} VCFs
specified in \code{files}. If \code{tumor.col.names} is equal to
\code{NA}(default), this function will use the 10th column in all the
\strong{Mutect} VCFs to calculate VAFs. See \code{\link{GetMutectVAF}} for
more details.}

\item{filter.status}{The status indicating a variant has passed all filters.
An example would be \code{"PASS"}. Variants which don't have the specified
\code{filter.status} in the \code{FILTER} column in VCF will be removed. If
\code{NULL}(default), no variants will be removed from the original VCF.}

\item{get.vaf.function}{Optional. Only applicable when \code{variant.caller} is
\strong{"unknown"}. Function to calculate VAF(variant allele frequency) and read
depth information from original VCF. See \code{\link{GetMutectVAF}} as an example.
If \code{NULL}(default) and \code{variant.caller} is "unknown", then VAF
and read depth will be NAs.}

\item{...}{Optional arguments to \code{get.vaf.function}.}

\item{max.vaf.diff}{\strong{Not} applicable if \code{variant.caller =
  "mutect"}. The maximum difference of VAF, default value is 0.02. If the
absolute difference of VAFs for adjacent SBSs is bigger than \code{max.vaf.diff},
then these adjacent SBSs are likely to be "merely" asynchronous single base
mutations, opposed to a simultaneous doublet mutation or variants involving
more than two consecutive bases.}

\item{suppress.discarded.variants.warnings}{Logical. Whether to suppress
warning messages showing information about the discarded variants. Default
is TRUE.}

\item{always.merge.SBS}{If \code{TRUE} merge adjacent SBSs as DBSs
regardless of VAFs and regardless of the value of \code{max.vaf.diff}
and regardless of the value of \code{get.vaf.function}. It is an
error to set this to \code{TRUE} when \code{variant.caller = "mutect"}.}
}
\description{
Read and split VCF files
}
\section{Value}{
 A list containing the following objects:
\itemize{
\item \code{SBS}: List of VCFs with only single base substitutions.
\item \code{DBS}: List of VCFs with only doublet base substitutions.
\item \code{ID}: List of VCFs with only small insertions and deletions.
\item \code{discarded.variants}: \strong{Non-NULL only if} there are variants
that were excluded from the analysis. See the added extra column
\code{discarded.reason} for more details.
}
}

\examples{
file <- c(system.file("extdata/Mutect-vcf",
                      "Mutect.GRCh37.s1.vcf",
                      package = "ICAMS"))
list.of.vcfs <- ReadAndSplitVCFs(file, variant.caller = "mutect")
}
\seealso{
\code{\link{VCFsToCatalogs}}
}
