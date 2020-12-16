file <- "D:/data/consensus-vcfs/DeIdentified-Sample.consensus.20160830.somatic.snv_mnv.vcf"

df <- ICAMS:::ReadVCF(file)
vcfs <- ICAMS:::SplitOneVCF(df)


foo <- abs(dt2$VAF.x - dt2$VAF.y)
