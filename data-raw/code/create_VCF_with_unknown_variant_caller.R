# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

df <- 
  ICAMS:::MakeDataFrameFromVCF("tests/testthat/testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf")

df1 <- df[, -(9:19)]

tmp <- sapply(1:nrow(df1), FUN = function(x){
  as.character(df1[x])
})

tmp1 <- apply(X = tmp, MARGIN = 2, FUN = function(x) {
  paste(x, collapse = "\t")
})

line1 <- "##fileformat=VCFv4.1"
line2 <- "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"
line3 <- "13	55736216	.	T	A	.	PASS	NT=ref;QSS=80;QSS_NT=80;SGT=TT->AT;SOMATIC;TQSS=1;TQSS_NT=1"
line4 <- "13	55736217	.	C	A	.	PASS	NT=ref;QSS=75;QSS_NT=75;SGT=CC->AC;SOMATIC;TQSS=1;TQSS_NT=1"
line5 <- "13	55736218	.	T	A	.	PASS	NT=ref;QSS=75;QSS_NT=75;SGT=TT->AT;SOMATIC;TQSS=1;TQSS_NT=1"
tmpfile <- tempfile(fileext = ".vcf")
write(line1, file = tmpfile, append = TRUE)
write(line2, file = tmpfile, append = TRUE)
write(tmp1, file = tmpfile, append = TRUE)
write(line3, file = tmpfile, append = TRUE)
write(line4, file = tmpfile, append = TRUE)
write(line5, file = tmpfile, append = TRUE)

file.copy(from = tmpfile, 
          to = "tests/testthat/testdata/SBS.GRCh37.variantcaller.unknown.vcf", 
          overwrite = TRUE)

file.remove(tmpfile)
