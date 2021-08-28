#' Get all the sequence contexts of the indels in a given 1 base-pair indel class.
#'
#' @param annotated.vcf An in-memory \code{data.frame} or similar table
#'  containing "VCF" (variant call format) data as created by
#'  \code{\link{VCFsToIDCatalogs}}
#'   with argument \code{return.annotated.vcfs = TRUE}.
#'
#' @param indel.class A single character string that denotes a 1 base pair
#' insertion or deletion, as taken from \code{ICAMS::catalog.row.order$ID}.
#' Insertions or deletions into / from 5+ base-pair homopolymers are
#' not supported.
#'
#' @return A list of all sequence contexts for the specified
#' \code{indel.class}.
#'
#' @export
ExtendSeqContextForOnebpINDEL <- function(annotated.vcf, indel.class){

  if(!indel.class %in% ICAMS::catalog.row.order$ID[c(1:5,7:11,13:17,19:23)]){
    stop("Argument indel.class value ", indel.class, " not supported")
  }

  if (!"ID.class" %in% colnames(annotated.vcf)) {
    stop("Argument annotated.vcf does not have columnn ID.class; ",
         "use ICAMS::VCFsToIDCatalogs with argument return.annotated.vcfs = TRUE")
  }

  annotated.vcf.this.class <-
    annotated.vcf[annotated.vcf$ID.class%in% indel.class,]

  ##extend the ref seq context from 13 to 21.
  annotated.vcf.this.class <-
    ICAMS::AnnotateIDVCF(ID.vcf = annotated.vcf.this.class,
                         ref.genome="hg19",
                         seq.context.width = 21)

  annotated.vcf.this.class <- annotated.vcf.this.class$annotated.vcf

  extended_sequence_context <-
    apply(annotated.vcf.this.class,1,
          function(x){
            Get1BPIndelFlanks(sequence = x["seq.context"],
                              ref = x["REF"],
                              alt = x["ALT"],
                              indel.class = indel.class)})

  return(extended_sequence_context)
}


#' Get all the sequence contexts of the indels in a given 1 base-pair indel class.
#'
#' @param sequence A string from \code{seq.context} column from in-memory \code{data.frame} or similar table
#'  containing "VCF" (variant call format) data as created by
#'  \code{\link{AnnotateIDVCF}}.
#'
#' @param ref A string from \code{REF} column from in-memory \code{data.frame} or similar table
#'  containing "VCF" (variant call format) data as created by
#'  \code{\link{AnnotateIDVCF}}.

#' @param alt A string from \code{ALT} column from in-memory \code{data.frame} or similar table
#'  containing "VCF" (variant call format) data as created by
#'  \code{\link{AnnotateIDVCF}}.
#'
#' @param indel.class A single character string that denotes a 1 base pair
#' insertion or deletion, as taken from \code{ICAMS::catalog.row.order$ID}.
#' Insertions or deletions into / from 5+ base-pair homopolymers are
#' not supported.
#'
#' @param flank.length The length of flanking bases from the region targeted by
#' insertion or deletion
#'
#' @return A string for the specified \code{sequence} and \code{indel.class}.
#'

Get1BPIndelFlanks <- function(sequence, ref,alt,indel.class, flank.length = 5){

  #cat(
  #  paste0('Get1BPIndelFlanks("', sequence, '", "', ref, "\",  \"", alt, "\",  \"",indel.class, "\")\n")
  #)

  # indel.class is string such as "DEL:T:1:3"
  split_indel.class <- unlist(strsplit(indel.class,":"))

  indel.base <- split_indel.class[2] # The base inserted or deleted (T or C)

  homopolymer.length <-  as.numeric(split_indel.class[4]) # a digit, 0:4

  #stopifnot(nchar(sequence) %% 2 == 1) # Must be an odd number. we don't need this one

  mid.base <- (nchar(sequence)+1)/2

  ins.or.del <- split_indel.class[1]


  if (ins.or.del == "INS" & homopolymer.length == 0) {
    if (nchar(alt) == 2) {
      alt <- substr(alt, 2, 2)
    }
    stopifnot(nchar(alt) == 1)
    seq.context <- substring(sequence, mid.base-flank.length + 1, mid.base+flank.length)

    # alt could be A, C, G, T, need to "normalize" to C or T
    if(alt != indel.base) {
      seq.context <-  ICAMS::revc(seq.context)
    }
  } else {

    if(ins.or.del == "DEL"){
      homopolymer.length <- homopolymer.length + 1
    }

    ##except for de nove insertion, we need to check if the ref is at the center


    if(nchar(alt) == 2 & substring(sequence,mid.base,mid.base)!= ref){ #means we are looking at an insertion
      stop("REF not at the center")
    }
    if(nchar(alt) == 1 & substring(sequence,mid.base,mid.base+1)!= ref){ #means we are looking at a deletion
      stop("REF not at the center")
    }



    homopolymer.starts <- mid.base+1

    homopolymer.ends <- mid.base+homopolymer.length

    ## normalize the insertion context to the middle

    seq.context <-  substring(sequence,homopolymer.starts-flank.length,homopolymer.ends+flank.length)

    if(substring(sequence,homopolymer.starts,homopolymer.starts)!= indel.base){
      seq.context <-  ICAMS::revc(substring(sequence,homopolymer.starts-flank.length,homopolymer.ends+flank.length))

    }
  }
  return(seq.context)

}





