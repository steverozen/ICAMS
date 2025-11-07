xx = readRDS("../CPCT02070023T.indel.vcf.rds")
yy = rbind(xx[1:3, ], xx[xx$POS == 36972193, ])
devtools::load_all()
zz = ICAMS::AnnotateIDVCF(yy, "hg19", flag.mismatches = T, explain_indels = T)
