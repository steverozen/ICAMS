test_that("annot_vcf_to_83_catalog produces correct output", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))

  id.vcf <- ICAMS:::ReadStrelkaIDVCF(
    testthat::test_path("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf")
  )

  annotated <- AnnotateIDVCF(
    id.vcf,
    ref.genome = "hg19",
    explain_indels = 0
  )

  result <- annot_vcf_to_83_catalog(
    annotated$annotated.vcf,
    sample_id = "s1"
  )

  # Check dimensions: 83 rows, 1 column
  expect_equal(nrow(result), 83)
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "s1")

  # Row names should match the canonical ID row order
  expect_equal(rownames(result), ICAMS::catalog.row.order$ID)

  # All counts should be non-negative integers
  expect_true(all(result >= 0))

  # Total counts should equal the number of unique positions in the VCF
  annot <- annotated$annotated.vcf
  n_unique <- length(unique(paste0(annot$CHROM, "-", annot$POS)))
  expect_equal(sum(result), n_unique)
  expect_snapshot(result)
})

test_that("annot_vcf_to_83_catalog clip_le_9 produces <= counts", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))

  id.vcf <- ICAMS:::ReadStrelkaIDVCF(
    testthat::test_path("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf")
  )

  annotated <- AnnotateIDVCF(
    id.vcf,
    ref.genome = "hg19",
    explain_indels = 0
  )

  full <- annot_vcf_to_83_catalog(
    annotated$annotated.vcf, sample_id = "s1"
  )
  clipped <- annot_vcf_to_83_catalog(
    annotated$annotated.vcf, sample_id = "s1", clip_le_9 = TRUE
  )

  expect_equal(nrow(clipped), 83)
  expect_true(all(clipped >= 0))
  expect_true(sum(clipped) <= sum(full))
})

test_that("annot_vcf_to_83_catalog produces correct output -- test 2", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))

  id.vcf <- ICAMS:::ReadMutectVCF(
    testthat::test_path("testdata/Mutect.GRCh38.vcf")
  )

  annotated <- AnnotateIDVCF(
    id.vcf,
    ref.genome = "hg38",
    explain_indels = 0
  )

  result <- annot_vcf_to_83_catalog(
    annotated$annotated.vcf,
    sample_id = "s1"
  )

  # Check dimensions: 83 rows, 1 column
  expect_equal(nrow(result), 83)
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "s1")

  # Row names should match the canonical ID row order
  expect_equal(rownames(result), ICAMS::catalog.row.order$ID)

  # All counts should be non-negative integers
  expect_true(all(result >= 0))

  # Total counts should equal the number of unique positions in the VCF
  annot <- annotated$annotated.vcf
  n_unique <- length(unique(paste0(annot$CHROM, "-", annot$POS)))
  expect_equal(sum(result), n_unique)
  expect_snapshot(result)
})
