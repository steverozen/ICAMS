#' This function makes catalogs from the sample Mutect VCF file
#' and compares it with the expected catalog information.
#' @keywords internal
TestMakeCatalogFromMutectVCFs <- function() {
  file <-
    "https://raw.githubusercontent.com/steverozen/ICAMS/master/data-raw/VCF/MCF10A_Carb_Low_cl2_Mutect.vcf"

  cats <-
    MutectVCFFilesToCatalog(file, ref.genome = "GRCh37",
                            # Use default transcript ranges
                            trans.ranges = trans.ranges.GRCh37,
                            region = "genome")

  prev.catalog.192 <-
    ReadCatalog(system.file("extdata",
                            "mutect.regress.cat.sbs.192.csv",
                            package = "ICAMS",
                            mustWork = TRUE),
                ref.genome = "GRCh37", region = "transcript", catalog.type = "counts")

  stopifnot(cats$catSBS192 == prev.catalog.192)

  prev.catalog.96 <-
    ReadCatalog(
      system.file("extdata",
                  "mutect.regress.cat.sbs.96.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", catalog.type = "counts")
  stopifnot(cats$catSBS96 == prev.catalog.96)

  prev.catalog.1536 <-
    ReadCatalog(
      system.file("extdata",
                  "mutect.regress.cat.sbs.1536.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", catalog.type = "counts")
  stopifnot(cats$catSBS1536 == prev.catalog.1536)

  prev.catalog.DBS.78<-
    ReadCatalog(
      system.file("extdata",
                  "mutect.regress.cat.dbs.78.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", catalog.type = "counts")
  stopifnot(cats$catDBS78 == prev.catalog.DBS.78)

  prev.catalog.DBS.136<-
    ReadCatalog(
      system.file("extdata",
                  "mutect.regress.cat.dbs.136.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", catalog.type = "counts")
  stopifnot(cats$catDBS136 == prev.catalog.DBS.136)

  prev.catalog.DBS.144<-
    ReadCatalog(
      system.file("extdata",
                  "mutect.regress.cat.dbs.144.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "transcript", catalog.type = "counts")
  stopifnot(cats$catDBS144 == prev.catalog.DBS.144)

  prev.catalog.indels<-
    ReadCatalog(
      system.file("extdata",
                  "mutect.regress.cat.indels.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", catalog.type = "counts")
  stopifnot(cats$catID == prev.catalog.indels)

  message("ok\n")

  invisible(cats)
}

#' This function is to make catalogs from the sample Strelka SBS VCF files
#' to compare with the expected catalog information.
#' @keywords internal
TestMakeCatalogFromStrelkaSBSVCFs <- function() {
  file <-
    "https://raw.githubusercontent.com/steverozen/ICAMS/master/data-raw/VCF/MCF10A_Carb_Low_cl2_Strelka_SBS.vcf"

  cats <-
    StrelkaSBSVCFFilesToCatalog(file, ref.genome = "GRCh37",
                                # Use default transcript ranges
                                trans.ranges = trans.ranges.GRCh37,
                                region = "genome")

  prev.catalog.192 <-
    ReadCatalog(system.file("extdata",
                            "strelka.regress.cat.sbs.192.csv",
                            package = "ICAMS",
                            mustWork = TRUE),
                ref.genome = "GRCh37", region = "transcript", catalog.type = "counts")

  stopifnot(cats$catSBS192 == prev.catalog.192)

  prev.catalog.96 <-
    ReadCatalog(
      system.file("extdata",
                  "strelka.regress.cat.sbs.96.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", catalog.type = "counts")
  stopifnot(cats$catSBS96 == prev.catalog.96)

  prev.catalog.1536 <-
    ReadCatalog(
      system.file("extdata",
                  "strelka.regress.cat.sbs.1536.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", catalog.type = "counts")
  stopifnot(cats$catSBS1536 == prev.catalog.1536)

  prev.catalog.DBS.78<-
    ReadCatalog(
      system.file("extdata",
                  "strelka.regress.cat.dbs.78.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", catalog.type = "counts")
  stopifnot(cats$catDBS78 == prev.catalog.DBS.78)

  prev.catalog.DBS.136<-
    ReadCatalog(
      system.file("extdata",
                  "strelka.regress.cat.dbs.136.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", catalog.type = "counts")
  stopifnot(cats$catDBS136 == prev.catalog.DBS.136)

  prev.catalog.DBS.144<-
    ReadCatalog(
      system.file("extdata",
                  "strelka.regress.cat.dbs.144.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "transcript", catalog.type = "counts")
  stopifnot(cats$catDBS144 == prev.catalog.DBS.144)

  message("ok\n")

  invisible(cats)
}

#' This function is to make catalogs from the sample Strelka ID VCF files
#' to compare with the expected catalog information.
#' @keywords internal
TestMakeCatalogFromStrelkaIDVCFs <- function() {
  file <-
    "https://raw.githubusercontent.com/steverozen/ICAMS/master/data-raw/VCF/MCF10A_Carb_Low_cl2_Strelka_ID.vcf"

  cat.ID <- StrelkaIDVCFFilesToCatalog(file, ref.genome = "GRCh37",
                                       region = "genome")

  prev.catalog.indels<-
    ReadCatalog(
      system.file("extdata",
                  "strelka.regress.cat.indels.csv",
                  package = "ICAMS",
                  mustWork = TRUE),
      ref.genome = "GRCh37", region = "genome", catalog.type = "counts")
  stopifnot(cat.ID == prev.catalog.indels)

  message("ok\n")

  invisible(cat.ID)
}
