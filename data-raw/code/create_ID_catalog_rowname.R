# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

# This function creates the ICAMS ID catalog rowname that is added as an extra
# column to the Excel file
# "./data-raw/PCAWG7_indel_classification_2021_09_03.xlsx"

types <- c("DEL", "INS")
bases <- c("T", "C")
repeat.sizes <- c(0:4, "5+")
deletion.sizes <- c(2:4, "5+")
homology.sizes <- c(1:4, "5+")

one.bp.indel.rownames <- character()
for (type in types) {
  for (size in repeat.sizes) {
    for (base in bases) {
      one.bp.indel.rownames <- 
          c(one.bp.indel.rownames, paste0(type, ":", base, ":", "1", ":", size))
    }
  }
}

repeats.indel.rownames <- character()
for (type in types) {
  for (deletion.size in deletion.sizes) {
    for (repeat.size in repeat.sizes) {
      repeats.indel.rownames <-
        c(repeats.indel.rownames, 
          paste0(type, ":repeats:", deletion.size, ":", repeat.size))
    }
  }
}
 
MH.indel.rownames <- character()                   
for (deletion.size in deletion.sizes) {
  for (homology.size in homology.sizes) {
    if (deletion.size != "5+") {
      if (homology.size < deletion.size) {
        MH.indel.rownames <-
          c(MH.indel.rownames, paste0("DEL:MH:", deletion.size, ":", homology.size))
      }
    } else {
      MH.indel.rownames <-
        c(MH.indel.rownames, paste0("DEL:MH:", deletion.size, ":", homology.size))
    }
  }
}

all.indel.rownames <- 
  c(one.bp.indel.rownames, repeats.indel.rownames, MH.indel.rownames)
setdiff(all.indel.rownames, ICAMS::catalog.row.order$ID)
setdiff(ICAMS::catalog.row.order$ID, all.indel.rownames)

df <- data.frame(ICAMS.ID.catalog.rownames = all.indel.rownames)
write.csv(df, file = "./data-raw/data/ICAMS.ID.catalog.rownames.csv")
