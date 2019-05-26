# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

list.of.vcf <- ReadAndSplitStrelkaSBSVCFs("data-raw/VCF/Strelka.SBS.GRCh37.vcf")
sbs.vcf <- list.of.vcf$SBS.vcfs[[1]]
strelka.SBS.vcf.GRCh37 <-
  AddSequence(sbs.vcf, ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5)

list.of.vcf <- ReadAndSplitStrelkaSBSVCFs("data-raw/VCF/Strelka.SBS.GRCh38.vcf")
sbs.vcf <- list.of.vcf$SBS.vcfs[[1]]
strelka.SBS.vcf.GRCh38 <-
  AddSequence(sbs.vcf, ref.genome = BSgenome.Hsapiens.UCSC.hg38)

save(strelka.SBS.vcf.GRCh37, strelka.SBS.vcf.GRCh38,
     file = "tests/testthat/testdata/test_AddSequence.Rdata")

id.vcf <- ReadStrelkaIDVCF("data-raw/VCF/Strelka.ID.GRCh37.vcf")
strelka.ID.vcf.GRCh37 <-
  AddAndCheckSequenceID(id.vcf, ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5)

id.vcf <- ReadStrelkaIDVCF("data-raw/VCF/Strelka.ID.GRCh38.vcf")
strelka.ID.vcf.GRCh38 <-
  AddAndCheckSequenceID(id.vcf, ref.genome = BSgenome.Hsapiens.UCSC.hg38)

save(strelka.ID.vcf.GRCh37, strelka.ID.vcf.GRCh38,
     file = "tests/testthat/testdata/test_AddAndCheckSequenceID.Rdata")
