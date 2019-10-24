# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

gene.expression.level.example.GRCh37 <- 
  read.table("data-raw/gene.expression.level.example.GRCh37.txt", header = TRUE)
colnames(gene.expression.level.example.GRCh37) <-
  c("Ensembl.gene.ID", "gene.symbol", "counts", "TPM")