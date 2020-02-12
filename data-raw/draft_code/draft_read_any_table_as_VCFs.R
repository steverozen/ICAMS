### ICAMS related functions needed for arbitrary? VCF / MAF like data formats

# Example

library(data.table)

d1 <- fread("data-raw/draft_code/sup-table-5.csv")

# assume microbiopsy_id is the ID we want

new.col.heads <- c("SampleID",
                   "data_access_id",
                   "histological_feature",
                   "CHROM",
                   "POS",
                   "REF",
                   "ALT",
                   "mut_depth",
                   "total_depth",
                   "VAF")

colnames(d1) <- new.col.heads

indel.rows <- which(nchar(d1$REF) != nchar(d1$ALT))
# Note, these data contain context for indels
# 
del.rows <- which(d1$ALT == "-") # none

di <- d1[indel.rows, ]

d1 <- as.data.frame(d1)

# vcfs <- split(d1, d1$SampleID)

# Since this contains indels will try

# split.vcfs <- ICAMS:::SplitListOfMutectVCFs(vcfs)

# Note, no DBSs -- possible problem dealing with VAFs -- not checked

d2 <- d1[-indel.rows, ]

vcfs <- split(d2, d2$SampleID)
id.vcfs <- split(di, di$SampleID)

split.vcfs2 <- ICAMS:::SplitListOfStrelkaSBSVCFs(vcfs)

# There _are_ DBS.
# 

library(BSgenome.Hsapiens.1000genomes.hs37d5)

SBS.cats <-
  ICAMS::VCFsToSBSCatalogs(
    split.vcfs2$SBS.vcfs, 
    ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
    region = "genome")


ID.cats <-
  ICAMS:::VCFsToIDCatalogs(
    id.vcfs,
    ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
    region = "genome",
    flag.mismatches = 1
  )
ICAMS::PlotCatalogToPdf(ID.cats$catalog, "bladder-indel.pdf")



SBS.ppms <- 
  ICAMS:::CreatePPMFromSBSVCFs(
    split.vcfs2$SBS.vcfs,
    ref.genome = 
      BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
    seq.context.width = 4)
ICAMS:::PlotPPMToPdf(list.of.ppm = SBS.ppms, 
                    file = "extended.pdf",
                     titles = names(split.vcfs2$SBS.vcfs))


d3 <- d2[d2$REF == "A" & d2$ALT != "C",  ]
avcfs <- split(d3, d3$SampleID)
a.ppmm <- ICAMS:::CreatePPMFromSBSVCFs(
  avcfs,
  ref.genome = 
    BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
  seq.context.width = 4)
ICAMS:::PlotPPMToPdf(list.of.ppm = a.ppmm, 
                     file = "A.to.A.or.Gextended.pdf",
                     titles = names(a.ppmm))



ee <- fread("data-raw/draft_code/bladder-exome.csv")
colnames(ee) <- new.col.heads
ee <- ee[nchar(ee$REF) == nchar(ee$ALT), ]
eevcfs <- vcfs <- split(ee, ee$SampleID)

