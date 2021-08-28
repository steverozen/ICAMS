Get1BPIndelFlanks("TCAGTGGGAAGAAGTGTTGGGGCAGGGTTGCAAAGACTCAGAA" , "GT", "INS:T:1:0")
Get1BPIndelFlanks("ACCCCCGCAGGACAGTCAGAGGACCCCGTCCAGACCTGAAGCC", "GA", "INS:T:1:1")
Get1BPIndelFlanks("ATTTAATTCTAAATTTCATATGAAAATATGAACGACAGAGAAT", "GA",  "INS:T:1:4")
Get1BPIndelFlanks("CTTGTTGTTTTAAGCCACTAGGTTTTAGGGTCTTTTGTTATGC", "GT",  "INS:T:1:4")
Get1BPIndelFlanks("AATTTAATTCTAAATTTCATATGAAAATATGAACGACAGAGAA", "AA",  "INS:T:1:4") # wrong
Get1BPIndelFlanks("GTTATTTGCTGAAACTCCTCATAAAATGAGTTTCTTTGCCAAA", "TA",  "T",  "INS:T:1:4")


Get1BPIndelFlanks("CGAGACTCCGTCTCAAAAAAATAAAAATAAAAATAAATAAAATT", "TA",  "T",  "DEL:T:1:4")
Get1BPIndelFlanks("TTTATCTATTGAGCTGTGTCTATTTTTGCATGTTTGGCAGTCAG", "A",  "DEL:T:1:4")
Get1BPIndelFlanks("AAAATAGTTGAAGAACAGGTGGTCAAGACTTTTCCAAATGTGAT", "G",  "DEL:T:1:0") # still even

load("vcfs.to.test.seq.context.function.Rdata")
source("~/monster.ID.project/seqcontext.to.deliver.R")



#################################################################
#test case for T insertion, indel vcf of HMF genome CPCT02130013T
#################################################################

# test.ins.vcf <- HMF_IDvcfs$CPCT02130013T
test.vcf <- HMF_IDvcfs$CPCT02130013T


zz <- ICAMS::VCFsToIDCatalogs(list(a = HMF_IDvcfs$CPCT02130013T), ref.genome = "hg19", return.annotated.vcfs = TRUE)
test.ins.vcf <- zz$annotated.vcfs[[1]]
test.ins.vcf <- test.ins.vcf[,c(1,2,4,5,12,13)]
test.ins.vcf <- test.ins.vcf[test.ins.vcf$ID.class %in% ICAMS::catalog.row.order$ID[c(1:5,7:11,13:17,19:23)],]


test.ins.vcf$CHROM <- "1"
test.ins.vcf$POS <- 1:nrow(test.ins.vcf)

extend.seq.context <- ExtendSeqContextForOnebpINDEL(annotated.vcf = test.ins.vcf,indel.class = "DEL:T:1:4")
#Generate a matrix for plotting
positions<-c(paste0("-",5:1),"T",paste0("+",1:5))

classes<-c("A","C","G","T")
output<-matrix(data=NA,nrow=length(positions),ncol=length(classes),
               dimnames = list(positions,classes))

for(row in 1:nrow(output)){
  tmp<-substr(extend.seq.context,row,row)
  output[row,"A"]<-sum(tmp == "A")
  output[row,"C"]<-sum(tmp == "C")
  output[row,"G"]<-sum(tmp == "G")
  output[row,"T"]<-sum(tmp == "T")
}
plotPWM(PWMmatrix = output,main = paste0("CPCT02130013T"," T>TT "," counts ",sum(output[1,])))

extend.seq.context <- ExtendSeqContextForOnebpINDEL(vcf = test.ins.vcf,indel_class = "INS:T:1:2")
#Generate a matrix for plotting
positions<-c(paste0("-",5:1),"T1","T2",paste0("+",1:5))

classes<-c("A","C","G","T")
output<-matrix(data=NA,nrow=length(positions),ncol=length(classes),
               dimnames = list(positions,classes))

for(row in 1:nrow(output)){
  tmp<-substr(extend.seq.context,row,row)
  output[row,"A"]<-sum(tmp == "A")
  output[row,"C"]<-sum(tmp == "C")
  output[row,"G"]<-sum(tmp == "G")
  output[row,"T"]<-sum(tmp == "T")
}
plotPWM(PWMmatrix = output,main = paste0("CPCT02130013T"," TT>TTT "," counts ",sum(output[1,])))


##########################################################################
#test case for de novo C insertion, indel vcf of HMF genome DRUP01340003T#
##########################################################################

extend.seq.context <- ExtendSeqContextForOnebpINDEL(vcf = HMF_IDvcfs$CPCT02050104T,indel_class = "INS:C:1:0")

#Generate a matrix for plotting

classes<-c("A","C","G","T")
positions<-c(paste0("-",5:1),paste0("+",1:5)) ##need to change manually

 output<-matrix(data=NA,nrow=length(positions),ncol=length(classes),
               dimnames = list(positions,classes))

for(row in 1:nrow(output)){
  tmp<-substr(extend.seq.context,row,row)
  output[row,"A"]<-sum(tmp == "A")
  output[row,"C"]<-sum(tmp == "C")
  output[row,"G"]<-sum(tmp == "G")
  output[row,"T"]<-sum(tmp == "T")
}
plotPWM(PWMmatrix = output,main = paste0("CPCT02050104T"," - > C "," counts ",sum(output[1,]))) ##title needs to be changed manually


