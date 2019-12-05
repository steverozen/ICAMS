# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

gene.expression.data.HepG2 <- 
  fread("data-raw/HepG2_SC1.txt", header = TRUE)
colnames(gene.expression.data.HepG2) <-
  c("Ensembl.gene.ID", "gene.symbol", "counts", "TPM")

gene.expression.data.MCF10A <- 
  fread("data-raw/MCF10A_SC1.txt", header = TRUE)
colnames(gene.expression.data.MCF10A) <-
  c("Ensembl.gene.ID", "gene.symbol", "counts", "TPM")