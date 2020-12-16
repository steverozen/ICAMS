# Test example for user

vcf.files <-
  dir(devtools::package_file("data-raw/VCF"),
      pattern = "Strelka.*SBS.*\\.vcf",
      full.names = TRUE)


StrelkaSBSVCFFilesToCatalogAndPlotToPdf(
  vcf.files,
  ref.genome = "hg19",
  output.file = devtools::package_file("data-raw/X.pdf")
)
