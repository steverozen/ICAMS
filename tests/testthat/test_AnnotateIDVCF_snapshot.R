context("AnnotateIDVCF snapshots")
helper_path <- "../../inst/xcategorize_1_justified_indel.R"
if (!file.exists(helper_path)) {
  helper_path <- "../../inst/prev_categorize_1_justified_indel.R"
}
if (file.exists(helper_path)) {
  source(helper_path)
}
library(ICAMS)

snapshot_data_frame <- function(dt, cols) {
  # Ensure snapshots capture plain data frames for stable printing
  as.data.frame(dt[, cols, with = FALSE])
}

test_that("AnnotateIDVCF snapshot with hg19", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))

  load("testdata/test_AnnotateIDVCF.Rdata")
  id.vcf <- ICAMS:::ReadStrelkaIDVCF(
    "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf"
  )

  result <- AnnotateIDVCF(
    id.vcf,
    ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
  )
  via_grch37 <- AnnotateIDVCF(id.vcf, ref.genome = "GRCh37")
  via_hg19 <- AnnotateIDVCF(id.vcf, ref.genome = "hg19")

  cols <- seq_len(29)
  snapshot_payload <- list(
    bsgenome = snapshot_data_frame(result$annotated.vcf, cols),
    grch37 = snapshot_data_frame(via_grch37$annotated.vcf, cols),
    hg19 = snapshot_data_frame(via_hg19$annotated.vcf, cols)
  )

  expect_snapshot(snapshot_payload)
})

test_that("AnnotateIDVCF snapshot with hg38", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))

  load("testdata/test_AnnotateIDVCF.Rdata")
  id.vcf <- ICAMS:::ReadStrelkaIDVCF("testdata/Strelka.ID.GRCh38.vcf")

  result <- AnnotateIDVCF(
    id.vcf,
    ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  )
  via_grch38 <- AnnotateIDVCF(id.vcf, ref.genome = "GRCh38")
  via_hg38 <- AnnotateIDVCF(id.vcf, ref.genome = "hg38")

  cols <- seq_len(ncol(result$annotated.vcf))
  snapshot_payload <- list(
    bsgenome = snapshot_data_frame(result$annotated.vcf, cols),
    grch38 = snapshot_data_frame(via_grch38$annotated.vcf, cols),
    hg38 = snapshot_data_frame(via_hg38$annotated.vcf, cols)
  )

  expect_snapshot(snapshot_payload)
})
