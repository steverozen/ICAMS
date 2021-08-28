test.vcf <- HMF_IDvcfs$CPCT02130013T


zz <- ICAMS::VCFsToIDCatalogs(list(a = HMF_IDvcfs$CPCT02130013T), ref.genome = "hg19", return.annotated.vcfs = TRUE)
test.ins.vcf <- zz$annotated.vcfs[[1]]
test.ins.vcf <- test.ins.vcf[,c(1,2,4,5,12,13)]
test.ins.vcf <- test.ins.vcf[test.ins.vcf$ID.class %in% ICAMS::catalog.row.order$ID[c(1:5,7:11,13:17,19:23)],]
write.table(test.ins.vcf,"~/ICAMS/data-raw/test.extendseqcontext.vcf",sep=",",col.names=T,row.names=F,quote=F)

