list.of.vcf <- ReadAndSplitStrelkaSNSVCFs("data-raw/Strelka.SNS.GRCh37.vcf")
sns.vcf <- list.of.vcf$SNS.vcfs[[1]]
strelka.SNS.vcf.GRCh37 <-
  AddSequence(sns.vcf, genome = BSgenome.Hsapiens.1000genomes.hs37d5)

list.of.vcf <- ReadAndSplitStrelkaSNSVCFs("data-raw/Strelka.SNS.GRCh38.vcf")
sns.vcf <- list.of.vcf$SNS.vcfs[[1]]
strelka.SNS.vcf.GRCh38 <-
  AddSequence(sns.vcf, genome = BSgenome.Hsapiens.UCSC.hg38)

save(strelka.SNS.vcf.GRCh37, strelka.SNS.vcf.GRCh38, file = "data-raw/test_AddSequence.Rdata")

id.vcf <- ReadStrelkaIDVCF("data-raw/Strelka.ID.GRCh37.vcf")
strelka.ID.vcf.GRCh37 <-
  AddAndCheckSequenceID(id.vcf, genome = BSgenome.Hsapiens.1000genomes.hs37d5)

id.vcf <- ReadStrelkaIDVCF("data-raw/Strelka.ID.GRCh38.vcf")
strelka.ID.vcf.GRCh38 <-
  AddAndCheckSequenceID(id.vcf, genome = BSgenome.Hsapiens.UCSC.hg38)

save(strelka.ID.vcf.GRCh37, strelka.ID.vcf.GRCh38,
     file = "data-raw/test_AddAndCheckSequenceID.Rdata")
