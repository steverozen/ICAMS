## code for plotting distance to TSS
if (FALSE) {
  ## calculate dist2TSS per mutation in the vcf, and binning this
  dist2TTSbin<-function(df,pos="POS",TSS="TSS"){
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
    return(df[,-colnames(df) == "dist2TSS"])
  }
  
  ## the ref.context and var.context columns are already oriented so the strand info is correct
  df$mutation<-paste0(substr(df$REF,6,6),">",substr(df$ALT,6,6))
  
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
  
  
  ## plotting
  barplot(t(as.data.frame(output)[,c("C>A","G>T")]),beside=T,xaxt="n",
          xlab=c("distance to TSS"),ylab="mutations (N)",main="dist2TSS C>A mutations")
  axis(1,at=c(2,5,8,11,14,17,20,23,26,29),labels = c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9",">9"),cex=0.8,tick=F,line=F,padj=-1)
  legend("topright",legend=c("Transcribed","nonTranscribed"),fill=c("grey80","black"),bty="n",cex=1)
  
  ## logistic regresssion
  df$newStrand <- NA
  df$newStrand[df$mutation %in% c("A>C","A>G","A>T","G>A","G>C","G>T")] <- 'sense'
  df$newStrand[df$mutation %in% c("C>A","C>G","C>T","T>A","T>C","T>G")] <- 'antisense'
  
  ## for all mutations together 
  fit<-glm(as.factor(df$newStrand) ~ df$TSS, family = 'binomial')
}

