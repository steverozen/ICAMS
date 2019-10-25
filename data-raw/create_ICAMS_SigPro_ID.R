# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

SigPro.to.ICAMS.ID <- read.csv("data-raw/SigPro.to.ICAMS.ID.csv", row.names = 1)
SigPro.to.ICAMS.ID <- as.matrix(SigPro.to.ICAMS.ID)


ICAMS.to.SigPro.ID <- read.csv("data-raw/ICAMS.to.SigPro.ID.csv", row.names = 1)
ICAMS.to.SigPro.ID <- as.matrix(ICAMS.to.SigPro.ID)