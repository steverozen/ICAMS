context("Transcription strand bias functions")

test_that("PlotTransBiasGeneExpToPdf function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh37.vcf")
  annotated.SBS.vcf <- 
    AnnotateSBSVCF(list.of.vcfs$SBS.vcfs[[1]], "hg19", trans.ranges.GRCh37)
  out <- 
    PlotTransBiasGeneExpToPdf(annotated.SBS.vcf = annotated.SBS.vcf,
                          file = file.path(tempdir(), "test1.pdf"),
                          expression.level = gene.expression.level.example.GRCh37, 
                          Ensembl.gene.ID.col = "Ensembl.gene.ID", TPM.col = "TPM",
                          num.of.bins = 4
                          )
  expect_equal(out, TRUE)
  unlink(file.path(tempdir(), "test1.pdf"))
})

test_that("PlotTransBiasDist2TSSToPdf function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh37.vcf")
  annotated.SBS.vcf <- 
    AnnotateSBSVCF(list.of.vcfs$SBS.vcfs[[1]], "hg19", trans.ranges.GRCh37)
  out <- PlotTransBiasDist2TSSToPdf(annotated.SBS.vcf, 
                                    plot.type = c("C>A", "C>G"),
                                    file = file.path(tempdir(), "test2.pdf"))
  expect_equal(out, TRUE)
  unlink(file.path(tempdir(), "test2.pdf"))
  
})

test_that("PlotTransBiasGeneExp function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh37.vcf")
  annotated.SBS.vcf <- 
    AnnotateSBSVCF(list.of.vcfs$SBS.vcfs[[1]], "hg19", trans.ranges.GRCh37)
  vcf1 <- annotated.SBS.vcf[REF == "C",]
  vcf2 <- annotated.SBS.vcf[REF %in% c("C", "G"), ]
  vcf3 <- annotated.SBS.vcf[REF == "T",]
  vcf4 <- annotated.SBS.vcf[REF %in% c("A", "T"), ]
  out <- 
    PlotTransBiasGeneExp(annotated.SBS.vcf = annotated.SBS.vcf, 
                         expression.level = gene.expression.level.example.GRCh37, 
                         Ensembl.gene.ID.col = "Ensembl.gene.ID", TPM.col = "TPM",
                         num.of.bins = 4, plot.type = "C>A")
  out1 <- 
    PlotTransBiasGeneExp(annotated.SBS.vcf = vcf1, 
                         expression.level = gene.expression.level.example.GRCh37, 
                         Ensembl.gene.ID.col = "Ensembl.gene.ID", TPM.col = "TPM",
                         num.of.bins = 4, plot.type = "C>A")
  
  out2 <- 
    PlotTransBiasGeneExp(annotated.SBS.vcf = vcf2, 
                         expression.level = gene.expression.level.example.GRCh37, 
                         Ensembl.gene.ID.col = "Ensembl.gene.ID", TPM.col = "TPM",
                         num.of.bins = 4, plot.type = "C>A")
  
  out3 <- 
    PlotTransBiasGeneExp(annotated.SBS.vcf = vcf3, 
                         expression.level = gene.expression.level.example.GRCh37, 
                         Ensembl.gene.ID.col = "Ensembl.gene.ID", TPM.col = "TPM",
                         num.of.bins = 4, plot.type = "C>A")
  
  out4 <- 
    PlotTransBiasGeneExp(annotated.SBS.vcf = vcf4, 
                         expression.level = gene.expression.level.example.GRCh37, 
                         Ensembl.gene.ID.col = "Ensembl.gene.ID", TPM.col = "TPM",
                         num.of.bins = 4, plot.type = "C>A")
  
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_equal(out3, TRUE)
  expect_equal(out4, TRUE)
})

test_that("PlotTransBiasDist2TSS function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh37.vcf")
  annotated.SBS.vcf <- 
    AnnotateSBSVCF(list.of.vcfs$SBS.vcfs[[1]], "hg19", trans.ranges.GRCh37)
  out <- PlotTransBiasDist2TSS(annotated.SBS.vcf, "C>A")
  expect_equal(out, TRUE)
  unlink("Rplots.pdf")
})
