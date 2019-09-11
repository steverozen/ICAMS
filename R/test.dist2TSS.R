dist2TSS <- function(annotated.SBS.vcf) {

#file <- c(system.file("extdata",
                       #"Strelka.SBS.GRCh37.vcf",
                       #package = "ICAMS"))
#list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs(file)
#SBS.vcf <- list.of.vcfs$SBS.vcfs[[1]]
#annotated.SBS.vcf <- AnnotateSBSVCF(SBS.vcf, ref.genome = "hg19", 
                                    #trans.ranges = trans.ranges.GRCh37)
df <- annotated.SBS.vcf
dist2TSSbin<-function(df,pos="POS",TSS="start"){
  df <- data.frame(df)
  df$dist2TSS<-abs(df[,pos] - df[,TSS])
  df$distBins<-NULL
  df$distBins[df$dist2TSS <= 100000]<-1
  df$distBins[df$dist2TSS <= 200000 & df$dist2TSS > 100000]<-2
  df$distBins[df$dist2TSS <= 300000 & df$dist2TSS > 200000]<-3
  df$distBins[df$dist2TSS <= 400000 & df$dist2TSS > 300000]<-4
  df$distBins[df$dist2TSS <= 500000 & df$dist2TSS > 400000]<-5
  df$distBins[df$dist2TSS <= 600000 & df$dist2TSS > 500000]<-6
  df$distBins[df$dist2TSS <= 700000 & df$dist2TSS > 600000]<-7
  df$distBins[df$dist2TSS <= 800000 & df$dist2TSS > 700000]<-8
  df$distBins[df$dist2TSS <= 900000 & df$dist2TSS > 800000]<-9
  df$distBins[df$dist2TSS <= 1000000 & df$dist2TSS > 900000]<-10
  return(df)
}

df <- dist2TSSbin(annotated.SBS.vcf)
df <- subset(df, df$distBins!='NA')

## the ref.context and var.context columns are already oriented so the strand info is correct
df$mutation<-paste0(df$REF,">",df$ALT)

## this part 
output<-matrix(data=NA,nrow=10,ncol=12)
rownames(output)<-c("group1","group2","group3","group4","group5","group6","group7","group8","group9","group10")

i<-table(df$mutation[df$distBins == 1])
colnames(output)<-names(i)
output[1,]<-i
for(class in colnames(output)){
  for(row in 2:10){
    output[row,class]<-sum(df$mutation[df$distBins == row] == class)
  }
}
output<-as.matrix(output)
output<-cbind(output,Sense=as.numeric(margin.table(output[,c(4:6,10:12)],1)))
output<-cbind(output,antiSense=as.numeric(margin.table(output[,c(1:3,7:9)],1)))

return(output)

}

