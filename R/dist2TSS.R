
# Check the current working directory
if (grepl("nat", getwd())) {
  vcfs <- ReadAndSplitStrelkaSBSVCFs("/Users/nat/Documents/ICAMS-dev/HepG2_Cis_1_SNVresult_rmDup2.vcf")
  sbs.vcf <- vcfs$SBS.vcfs[[1]]
  annotated.vcf <- AnnotateSBSVCF(SBS.vcf = sbs.vcf, ref.genome = "hg19", 
                                  trans.ranges = trans.ranges.GRCh37)
}


