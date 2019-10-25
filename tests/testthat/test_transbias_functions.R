context("Transcription strand bias functions")

test_that("PlotTransBiasExpToPdf function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh37.vcf")
  annotated.SBS.vcf <- 
    AnnotateSBSVCF(list.of.vcfs$SBS.vcfs[[1]], "hg19", trans.ranges.GRCh37)
  out <- 
    PlotTransBiasExpToPdf(annotated.SBS.vcf = annotated.SBS.vcf,
                          file = file.path(tempdir(), "test1.pdf"),
                          expression.level = gene.expression.level.example.GRCh37, 
                          Ensembl.gene.ID.col = "Ensembl.gene.ID", TPM.col = "TPM",
                          num.of.bins = 4
                          )
  expect_equal(out, TRUE)
  
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
  
})

test_that("PlotTransBiasExp function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh37.vcf")
  annotated.SBS.vcf <- 
    AnnotateSBSVCF(list.of.vcfs$SBS.vcfs[[1]], "hg19", trans.ranges.GRCh37)
  out <- 
    PlotTransBiasExp(annotated.SBS.vcf = annotated.SBS.vcf, 
                     expression.level = gene.expression.level.example.GRCh37, 
                     Ensembl.gene.ID.col = "Ensembl.gene.ID", TPM.col = "TPM",
                     num.of.bins = 4, plot.type = "C>A")
  expect_equal(out, TRUE)
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
unlink(file.path(tempdir(), "test1.pdf"))
unlink(file.path(tempdir(), "test2.pdf"))
