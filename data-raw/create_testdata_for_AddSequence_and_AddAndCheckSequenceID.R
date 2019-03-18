# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

list.of.vcf <- ReadAndSplitStrelkaSNSVCFs("data-raw/VCF/Strelka.SNS.GRCh37.vcf")
sns.vcf <- list.of.vcf$SNS.vcfs[[1]]
strelka.SNS.vcf.GRCh37 <-
  AddSequence(sns.vcf, genome = BSgenome.Hsapiens.1000genomes.hs37d5)

list.of.vcf <- ReadAndSplitStrelkaSNSVCFs("data-raw/VCF/Strelka.SNS.GRCh38.vcf")
sns.vcf <- list.of.vcf$SNS.vcfs[[1]]
strelka.SNS.vcf.GRCh38 <-
  AddSequence(sns.vcf, genome = BSgenome.Hsapiens.UCSC.hg38)

save(strelka.SNS.vcf.GRCh37, strelka.SNS.vcf.GRCh38,
     file = "tests/testthat/testdata/test_AddSequence.Rdata")

id.vcf <- ReadStrelkaIDVCF("data-raw/VCF/Strelka.ID.GRCh37.vcf")
strelka.ID.vcf.GRCh37 <-
  AddAndCheckSequenceID(id.vcf, genome = BSgenome.Hsapiens.1000genomes.hs37d5)

id.vcf <- ReadStrelkaIDVCF("data-raw/VCF/Strelka.ID.GRCh38.vcf")
strelka.ID.vcf.GRCh38 <-
  AddAndCheckSequenceID(id.vcf, genome = BSgenome.Hsapiens.UCSC.hg38)

save(strelka.ID.vcf.GRCh37, strelka.ID.vcf.GRCh38,
     file = "tests/testthat/testdata/test_AddSequence.Rdata")
