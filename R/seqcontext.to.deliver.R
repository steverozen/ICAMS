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
            Get1BPIndelFlanks(x["seq.context"],
                                x["ALT"],
                                indel.class = indel.class)})
  
  return(extended_sequence_context)
}


##input_line a line of ICAMS vcf. must contain "ALT" and "seq.context"
##indel.class an ICAMS catalog row order

Get1BPIndelFlanks <- function(sequence, alt, indel.class, flank.length = 5){
  
  cat(
    paste0('Get1BPIndelFlanks("', sequence, '", "', alt, "\",  \"", indel.class, "\")\n")
  )
  
  # indel.class is string such as "DEL:T:1:3"
  split_indel.class <- unlist(strsplit(indel.class,":"))

  indel_base <- split_indel.class[2] # The base inserted or deleted (T or C)

  homopolymer_length <-  as.numeric(split_indel.class[4]) # a digit, 0:4

  stopifnot(nchar(sequence) %% 2 == 1) # Must be an odd number
  mid.base <- (nchar(sequence)+1)/2
  
  ins.or.del <- split_indel.class[1]
  
  if (ins.or.del == "INS" & homopolymer_length == 0) {
    if (nchar(alt) == 2) {
      alt <- substr(alt, 2, 2)
    }
    stopifnot(nchar(alt) == 1)
    seq.context <- substring(sequence, mid.base-flank.length + 1, mid.base+flank.length)
    
    # alt could be A, C, G, T, need to "normalize" to C or T
    if(alt != indel_base) { 
      seq.context <-  ICAMS::revc(seq.context)
    }
  } else {

    homopolymer.starts <- mid.base+1

    homopolymer.ends <- mid.base+homopolymer_length

    ## normalize the insertion context to the middle
    if(substring(sequence,homopolymer.starts,homopolymer.starts)== ICAMS::revc(indel_base)){
      seq.context <-  ICAMS::revc(substring(sequence,homopolymer.starts-flank.length,homopolymer.ends+flank.length))

    }else{
      seq.context <-  substring(sequence,homopolymer.starts-flank.length,homopolymer.ends+flank.length)
    }

  }
  return(seq.context)

}







plotPWM<-function(PWMmatrix,main){

  number.of.rows <- nrow(PWMmatrix)

  if((number.of.rows %% 2) != 0){
    x <- c(-((number.of.rows-1)/2):((number.of.rows-1)/2))

  }else{
    x <- c((1-number.of.rows/2):(number.of.rows/2))
  }



  plot(x,PWMmatrix[,"A"]/sum(PWMmatrix[1,]),
       main=main,
       xlab="",
       ylab="frequency",xaxt="n",
       col="darkgreen",ylim=c(0,1),type="b",pch=20,lwd=2,cex.main=1)
  axis(1, at=x,labels=rownames(PWMmatrix),
       tick=F,outer=F,las=2,font=1,par(cex.axis=1))

  lines(x,PWMmatrix[,"C"]/sum(PWMmatrix[1,]),col="blue",type="b",pch=20,lwd=2)
  lines(x,PWMmatrix[,"G"]/sum(PWMmatrix[1,]),col="black",type="b",pch=20,lwd=2)
  lines(x,PWMmatrix[,"T"]/sum(PWMmatrix[1,]),col="red",type="b",pch=20,lwd=2)
  abline(h=0.25,col="grey50")
  legend("topright",pch=16,cex=1,ncol=4,bty="n",
         legend=c("A","C","G","T"),col=c("darkgreen","blue","black","red"))

}




