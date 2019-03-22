#' This function makes catalogs from the sample Mutect VCF file
#' and compares it with the expected catalog information.
#' @keywords internal
TestMakeCatalogFromMutectVCFs <- function() {
  files <- c(system.file("extdata",
                         "MCF10A_Carb_Low_cl2_Mutect.vcf",
                         package = "ICAMS",
                         mustWork = TRUE))

  cats <-
    MutectVCFFilesToCatalog(
      files,
      genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
      # Use default transcript ranges
      trans.ranges = trans.ranges.GRCh37)

  prev.catalog.192 <-
    ReadCatalog(system.file("extdata",
                            "mutect.regress.cat.sns.192.csv",
                            package = "ICAMS",
                            mustWork = TRUE),
                ref.genome = "GRCh37", region = "genome", type = "counts")

  stopifnot(cats$catSNS192 == prev.catalog.192$catalog)

  prev.catalog.96 <-
    ReadCatalog(
      system.file("extdata",
                  "mutect.regress.cat.sns.96.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", type = "counts")
  stopifnot(cats$catSNS96 == prev.catalog.96$catalog)

  prev.catalog.1536 <-
    ReadCatalog(
      system.file("extdata",
                  "mutect.regress.cat.sns.1536.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", type = "counts")
  stopifnot(cats$catSNS1536 == prev.catalog.1536$catalog)

  prev.catalog.DNS.78<-
    ReadCatalog(
      system.file("extdata",
                  "mutect.regress.cat.dns.78.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", type = "counts")
  stopifnot(cats$catDNS78 == prev.catalog.DNS.78$catalog)

  prev.catalog.DNS.136<-
    ReadCatalog(
      system.file("extdata",
                  "mutect.regress.cat.dns.136.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", type = "counts")
  stopifnot(cats$catQUAD136 == prev.catalog.DNS.136$catalog)

  prev.catalog.DNS.144<-
    ReadCatalog(
      system.file("extdata",
                  "mutect.regress.cat.dns.144.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", type = "counts")
  stopifnot(cats$catDNS144 == prev.catalog.DNS.144$catalog)

  prev.catalog.indels<-
    ReadCatalog(
      system.file("extdata",
                  "mutect.regress.cat.indels.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", type = "counts")
  stopifnot(cats$catID == prev.catalog.indels$catalog)

  cat("ok\n")

  invisible(cats)
}

#' This function is to make catalogs from the sample Strelka SNS VCF files
#' to compare with the expected catalog information.
#' @keywords internal
TestMakeCatalogFromStrelkaSNSVCFs <- function() {
  files <- c(system.file("extdata",
                         "MCF10A_Carb_Low_cl2_Strelka_SNS.vcf",
                         package = "ICAMS",
                         mustWork = TRUE))

  cats <-
    StrelkaSNSVCFFilesToCatalog(
      files,
      genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
      # Use default transcript ranges
      trans.ranges = trans.ranges.GRCh37)

  prev.catalog.192 <-
    ReadCatalog(system.file("extdata",
                            "new.regress.cat.sns.192.csv",
                            package = "ICAMS",
                            mustWork = TRUE),
                ref.genome = "GRCh37", region = "genome", type = "counts")

  stopifnot(cats$catSNS192 == prev.catalog.192$catalog)

  prev.catalog.96 <-
    ReadCatalog(
      system.file("extdata",
                  "new.regress.cat.sns.96.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", type = "counts")
  stopifnot(cats$catSNS96 == prev.catalog.96$catalog)

  prev.catalog.1536 <-
    ReadCatalog(
      system.file("extdata",
                  "new.regress.cat.sns.1536.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", type = "counts")
  stopifnot(cats$catSNS1536 == prev.catalog.1536$catalog)

  prev.catalog.DNS.78<-
    ReadCatalog(
      system.file("extdata",
                  "new.regress.cat.dns.78.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", type = "counts")
  stopifnot(cats$catDNS78 == prev.catalog.DNS.78$catalog)

  prev.catalog.DNS.136<-
    ReadCatalog(
      system.file("extdata",
                  "new.regress.cat.dns.136.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", type = "counts")
  stopifnot(cats$catQUAD136 == prev.catalog.DNS.136$catalog)

  prev.catalog.DNS.144<-
    ReadCatalog(
      system.file("extdata",
                  "new.regress.cat.dns.144.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", type = "counts")
  stopifnot(cats$catDNS144 == prev.catalog.DNS.144$catalog)

  cat("ok\n")

  invisible(cats)
}

#' This function is to make catalogs from the sample Strelka ID VCF files
#' to compare with the expected catalog information.
#' @keywords internal
TestMakeCatalogFromStrelkaIDVCFs <- function() {
  files <- c(system.file("extdata",
                         "MCF10A_Carb_Low_cl2_INDELresult.vcf",
                         package = "ICAMS",
                         mustWork = TRUE))

  cat.ID <-
    StrelkaIDVCFFilesToCatalog(
      files,
      genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5)

  prev.catalog.indels<-
    ReadCatID(
      system.file("extdata",
                  "new.regress.cat.indels.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cat.ID == prev.catalog.indels)

  cat("ok\n")

  invisible(cat.ID)
}
