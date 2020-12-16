# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

gene.expression.data.HepG2 <- 
  fread("data-raw/data/HepG2_SC1.txt", header = TRUE)
colnames(gene.expression.data.HepG2) <-
  c("Ensembl.gene.ID", "gene.symbol", "counts", "TPM")

gene.expression.data.MCF10A <- 
  fread("data-raw/data/MCF10A_SC1.txt", header = TRUE)
colnames(gene.expression.data.MCF10A) <-
  c("Ensembl.gene.ID", "gene.symbol", "counts", "TPM")