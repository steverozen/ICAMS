# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

file <- "tests/testthat/testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.vcf"
file1 <- "tests/testthat/testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.vcf"

df.SBS <- read.csv(file, header = FALSE, sep = "\t", quote = "",
               col.names = paste0("c", 1:100), as.is = TRUE)
# Delete the columns which are totally empty
df.SBS <- df[!sapply(df, function(x) all(is.na(x)))]

df.ID <- MakeDataFrameFromVCF(file1)
colnames(df.ID) <- paste0("c", 1:19)

df.mixed <- rbind(df.SBS, df.ID)

write.table(x = df.mixed, file = "tests/testthat/testdata/Strelka.mixed.GRCh37.vcf", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
