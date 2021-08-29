#' Get all the sequence contexts of the indels in a given 1 base-pair indel class from a VCF.
#'
#' @param annotated.vcf An in-memory \code{data.frame} or similar table
#'  containing "VCF" (variant call format) data as created by
#'  \code{\link{VCFsToIDCatalogs}}
#'   with argument \code{return.annotated.vcfs = TRUE}.
#'
#' @param indel.class A single character string that denotes a 1 base pair
#' insertion or deletion, as taken from \code{ICAMS::catalog.row.order$ID}.
#' Insertions or deletions into or from 5+ base-pair homopolymers are
#' not supported.
#'
#' @param flank.length The length of flanking bases around the position
#' or homopolymer targeted by the indel.
#'
#' @return A list of all sequence contexts for the specified
#' \code{indel.class}.
#'
#' @export
SymmetricalContextsFor1BPIndel <- function(annotated.vcf, indel.class, flank.length = 5){

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
            Get1BPIndelFlanks(sequence     = x["seq.context"],
                              ref          = x["REF"],
                              alt          = x["ALT"],
                              indel.class  = indel.class,
                              flank.length = flank.length)})

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
#' @param flank.length The length of flanking bases around the position
#' or homopolymer targeted by the indel.
#'
#' @return A string for the specified \code{sequence} and \code{indel.class}.
#'

Get1BPIndelFlanks <- function(sequence, ref, alt, indel.class, flank.length = 5){

  # Un-comment the following to generate test function calls.
  # cat(
  #  paste0('Get1BPIndelFlanks("', sequence, '", "', ref, "\",  \"", alt, "\",  \"",indel.class, "\")\n")
  # )

  # sanity check; assume the input VCF provides one base of context to the left of the indel,
  # e.g. insertion A -> AT, deletion CT -> C
  stopifnot(nchar(ref) + nchar(alt) == 3)

  # indel.class is string such as "DEL:T:1:3"
  split_indel.class <- unlist(strsplit(indel.class,":"))

  indel.base <- split_indel.class[2] # The base inserted or deleted (T or C)

  homopolymer.length <-  as.numeric(split_indel.class[4])
  stopifnot(homopolymer.length %in% 0:4)

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

    ##except for de novo insertion, we need to check if the ref is at the center

    # I don't think this is the only check we need
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

    homopolymer.seq <- paste(rep(indel.base, homopolymer.length), collapse = "")
    re <- paste0("[ACGT]{", flank.length - 1, "}[^", indel.base, "]", homopolymer.seq,
                 "[^", indel.base, "][ACGT]{", flank.length - 1, "}")
    if (!grepl(re, seq.context, perl = TRUE)) {
      stop("Extracted sequence ", seq.context, " does not have the expected form ",
           "(does not match the RE '", re, "')\n",
           "Possibly the variant caller is not standardizing the position of the indel in the homopolymer")
    }

  }
  return(seq.context)

}

#' Generate PFMmatrix from a given list of sequences.
#'
#' @param sequences A list of strings returned from \code{\link{SymmetricalContextsFor1BPIndel}}.
#'
#' @param indel.class A single character string that denotes a 1 base pair
#' insertion or deletion, as taken from \code{ICAMS::catalog.row.order$ID}.
#' Insertions or deletions into or from 5+ base-pair homopolymers are
#' not supported.
#'
#' @param flank.length The length of flanking bases around the position
#' or homopolymer targeted by the indel.+
#'
#' @param plot.dir If provided, make a dot-line plot for PFMmatrix.
#'
#' @param plot.title The title of the dot-line plot
#'
#' @return A matrix recording the frequency of each base (A, C, G, T) on each position of the sequence.
#' @export
GeneratePlotPFMmatrix <- function(sequences,flank.length = 5,indel.class,plot.dir=NULL,plot.title=NULL){

  if(length(unique(nchar(sequences)))>1){
    stop("All sequences must have the same length")
  }

  indel.base <- unlist(strsplit(indel.class,":"))[2]

  indel.context <- as.numeric(unlist(strsplit(indel.class,":"))[4])

  target.seq <- NULL

  if(indel.context>0){target.seq <- paste0(indel.base,1:indel.context)}

  positions <- c(paste0("-",(flank.length:1)),
                 paste0("+",1:flank.length))

  if(!is.null(target.seq)){
    positions <- c(paste0("-",(flank.length:1)),
                   target.seq,
                   paste0("+",1:flank.length))
  }

  classes<-c("A","C","G","T")
  PFMmatrix<-matrix(data=NA,nrow=length(positions),ncol=length(classes),
                    dimnames = list(positions,classes))

  for(row in 1:nrow(PFMmatrix)){
    tmp<-substr(sequences,row,row)
    PFMmatrix[row,"A"]<-sum(tmp == "A")
    PFMmatrix[row,"C"]<-sum(tmp == "C")
    PFMmatrix[row,"G"]<-sum(tmp == "G")
    PFMmatrix[row,"T"]<-sum(tmp == "T")
  }

  if(!is.null(plot.dir)){

    plot.title <- "Dot-line plot for PFMmatrix"

    grDevices::pdf(plot.dir)

    if(!is.null(plot.title)){
      dot.line.plot <- PlotPFMmatrix(PFMmatrix = PFMmatrix,
                                     title = plot.title)
    }

    grDevices::dev.off()
  }
  return(PFMmatrix)



}
#' Generate dot-line plot for sequence contest of 1bp indel.
#'
#' @param PFMmatrix An object return from \code{\link{GeneratePlotPFMmatrix}}.
#'
#' @param title A string provides the title of the plot
#'
#' @return An \strong{invisible} list.
#'
PlotPFMmatrix<-function(PFMmatrix,title){

  number.of.rows <- nrow(PFMmatrix)

  #set x-axis positions for the plot
  x <- c((1-number.of.rows/2):(number.of.rows/2))

  if((number.of.rows %% 2) != 0){
    x <- c(-((number.of.rows-1)/2):((number.of.rows-1)/2))

  }

  plot(x,PFMmatrix[,"A"]/sum(PFMmatrix[1,]),
       main=title,
       xlab="",
       ylab="frequency",xaxt="n",
       col="darkgreen",ylim=c(0,1),type="b",pch=20,lwd=2,cex.main=1)
  axis(1, at=x,labels=rownames(PFMmatrix),
       tick=F,outer=F,las=2,font=1,par(cex.axis=1))

  lines(x,PFMmatrix[,"C"]/sum(PFMmatrix[1,]),col="blue",type="b",pch=20,lwd=2)
  lines(x,PFMmatrix[,"G"]/sum(PFMmatrix[1,]),col="black",type="b",pch=20,lwd=2)
  lines(x,PFMmatrix[,"T"]/sum(PFMmatrix[1,]),col="red",type="b",pch=20,lwd=2)
  abline(h=0.25,col="grey50")
  legend("topright",pch=16,cex=1,ncol=4,bty="n",
         legend=c("A","C","G","T"),col=c("darkgreen","blue","black","red"))


}




