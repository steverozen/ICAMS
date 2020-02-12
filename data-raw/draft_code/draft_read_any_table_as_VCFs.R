### ICAMS related functions needed for arbitrary? VCF / MAF like data formats

# Example

library(data.table)



# General workflow for a file with multiple samples
# 
# 
# Read the file

d1 <- fread("....path to file here.....csv")


# Make sure the essential columns are named as correctly
# (manual, will vary from file to file)


new.col.heads <- c("SampleID",
                   "data_access_id",
                   "field_we_dont_care_about",
                   "CHROM",
                   "POS",
                   "REF",
                   "ALT",
                   "field2",
                   "field3",
                   "VAF")

colnames(d1) <- new.col.heads


# check for presence of indels
# If they exist, need to check if there is context for the insertions
# If not, can only try one's best
indel.rows <- which(nchar(d1$REF) != nchar(d1$ALT))
# Note, these data contain context for indels
# 
del.rows <- which(d1$ALT == "-") # none

# di <- d1[indel.rows, ]

d1 <- as.data.frame(d1)

# vcfs <- split(d1, d1$SampleID)

# Since this contains indels will try

# Check if there are DBSs etc (not shown)
# 
# Check for complex indels
# 
# split.vcfs <- ICAMS:::SplitListOfMutectVCFs(vcfs)

# Note, no DBSs -- possible problem dealing with VAFs -- not checked

d2 <- d1[-indel.rows, ]


# get  a list of VCFs
vcfs <- split(d2, d2$SampleID)

# id.vcfs <- split(di, di$SampleID)


# This one will find DBSs if necessary
# Probably the VAF column has to be present, so
# if one wants this and there is no VAF probably just set all varfs to 50%
split.vcfs2 <- ICAMS:::SplitListOfStrelkaSBSVCFs(vcfs)

# There _are_ DBS, so this will find them
# 

library(BSgenome.Hsapiens.1000genomes.hs37d5)


# generate the SBS catalogs
SBS.cats <-
  ICAMS::VCFsToSBSCatalogs(
    split.vcfs2$SBS.vcfs, 
    ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
    region = "genome")


# not shown, generate the DBS catalogs


# Generate the id catalogs
# ID.cats <-
#  ICAMS:::VCFsToIDCatalogs(
#    id.vcfs,
#    ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
#   region = "genome",
#    flag.mismatches = 1
#  )
# ICAMS::PlotCatalogToPdf(ID.cats$catalog, "indel.pdf")


# get ppms
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
                     file = "A.to.A.or.G.extended.pdf",
                     titles = names(a.ppmm))


# Plot everything
# 