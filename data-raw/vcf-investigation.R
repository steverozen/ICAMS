
df <- MakeDataFrameFromMutectVCF(
  "../OTA_vcf/0-5mM-OTA-S9-1_S7.bam.recal_calls.vcf")
dt <- data.table(df)
dt1 <- dt[nchar(REF) >= 2 & nchar(ALT) >= 2]
View(dt1)

dt$VAF <- GetMutectVAF(dt)
dt.SBS <- dt[nchar(REF) == 1 & nchar(ALT) == 1]
split.vcf <- SplitStrelkaSBSVCF(dt.SBS)
View(split.vcf$DBS.vcf)
