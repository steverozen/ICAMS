# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

list.of.vcf <- ReadAndSplitStrelkaSBSVCFs("data-raw/VCF/Strelka.SBS.GRCh37.vcf")
sbs.vcf <- list.of.vcf$SBS.vcfs[[1]]
strelka.SBS.vcf.GRCh37 <- AddSeqContext(sbs.vcf, ref.genome = "hg19")

list.of.vcf <- ReadAndSplitStrelkaSBSVCFs("data-raw/VCF/Strelka.SBS.GRCh38.vcf")
sbs.vcf <- list.of.vcf$SBS.vcfs[[1]]
strelka.SBS.vcf.GRCh38 <- AddSeqContext(sbs.vcf, ref.genome = "hg38")

list.of.vcf <- ReadAndSplitStrelkaSBSVCFs("data-raw/VCF/Strelka.SBS.GRCm38.vcf")
sbs.vcf <- list.of.vcf$SBS.vcfs[[1]]
strelka.SBS.vcf.GRCm38 <- AddSeqContext(sbs.vcf, ref.genome = "mm10")

save(strelka.SBS.vcf.GRCh37, strelka.SBS.vcf.GRCh38, strelka.SBS.vcf.GRCm38,
     file = "tests/testthat/testdata/test_AddSeqContext.Rdata")

id.vcf <- ReadStrelkaIDVCF("data-raw/VCF/Strelka.ID.GRCh37.vcf")
strelka.ID.vcf.GRCh37 <- AnnotateIDVCF(id.vcf, ref.genome = "hg19")[[1]]

id.vcf <- ReadStrelkaIDVCF("data-raw/VCF/Strelka.ID.GRCh38.vcf")
strelka.ID.vcf.GRCh38 <- AnnotateIDVCF(id.vcf, ref.genome = "hg38")[[1]]

id.vcf <- ReadStrelkaIDVCF("data-raw/VCF/Strelka.ID.GRCm38.vcf")
strelka.ID.vcf.GRCm38 <- AnnotateIDVCF(id.vcf, ref.genome = "mm10")[[1]]

save(strelka.ID.vcf.GRCh37, strelka.ID.vcf.GRCh38, strelka.ID.vcf.GRCm38,
     file = "tests/testthat/testdata/test_AnnotateIDVCF.Rdata")
