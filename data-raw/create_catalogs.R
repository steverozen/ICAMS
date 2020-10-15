make.row.order1536 <- function() {
  # Make the row order of 1536 pentanucleotide mutation types.

  all4 <- c("A", "C", "G", "T")

  retval1 <- character(1536)
  retval2 <- character(1536)
  i <- 1

  for (center in c("C", 'T')) {
    for (alt in setdiff(all4, center)) {
      for (minus2 in all4) {
        for (minus1 in all4) {
          for (plus1 in all4) {
            for (plus2 in all4) {
              retval1[i] <-
                paste0(minus2, minus1, center, plus1, plus2, alt)
              retval2[i] <-
                paste0(center, alt, minus2, minus1, plus1, plus2)
              i <- i + 1
            }
          }
        }
      }
    }
  }
  return(list(standard = retval1, mini = retval2))
}

make.row.order.sp1536 <- function() {
  # Make the row order of 1536 pentanucleotide mutation types in SigProfiler
  # format
  all4 <- c("A", "C", "G", "T")
  retval <- character(1536)
  i <- 1
  
  for (minus2 in all4) {
    for (minus1 in all4) {
      for (center in c("C", 'T')) {
        for (alt in setdiff(all4, center)) {
          for (plus1 in all4) {
            for (plus2 in all4) {
              retval[i] <-
                paste0(minus2, minus1, "[", center, ">", alt, "]", plus1, plus2)
              i <- i + 1
            }
          }
        }
      }
    }
  }
  return(retval)
}

catalog.row.order.SBS.1536 <- make.row.order1536()$standard
catalog.row.order.sp.SBS.1536 <- make.row.order.sp1536()

catalog.row.order.SBS.192 <-
  c("AAAC","AACC","AAGC","AATC","CAAC","CACC","CAGC","CATC","GAAC",
    "GACC","GAGC","GATC","TAAC","TACC","TAGC","TATC","AAAG","AACG",
    "AAGG","AATG","CAAG","CACG","CAGG","CATG","GAAG","GACG","GAGG",
    "GATG","TAAG","TACG","TAGG","TATG","AAAT","AACT","AAGT","AATT",
    "CAAT","CACT","CAGT","CATT","GAAT","GACT","GAGT","GATT","TAAT",
    "TACT","TAGT","TATT","ACAA","ACCA","ACGA","ACTA","CCAA","CCCA",
    "CCGA","CCTA","GCAA","GCCA","GCGA","GCTA","TCAA","TCCA","TCGA",
    "TCTA","ACAG","ACCG","ACGG","ACTG","CCAG","CCCG","CCGG","CCTG",
    "GCAG","GCCG","GCGG","GCTG","TCAG","TCCG","TCGG","TCTG","ACAT",
    "ACCT","ACGT","ACTT","CCAT","CCCT","CCGT","CCTT","GCAT","GCCT",
    "GCGT","GCTT","TCAT","TCCT","TCGT","TCTT","AGAA","AGCA","AGGA",
    "AGTA","CGAA","CGCA","CGGA","CGTA","GGAA","GGCA","GGGA","GGTA",
    "TGAA","TGCA","TGGA","TGTA","AGAC","AGCC","AGGC","AGTC","CGAC",
    "CGCC","CGGC","CGTC","GGAC","GGCC","GGGC","GGTC","TGAC","TGCC",
    "TGGC","TGTC","AGAT","AGCT","AGGT","AGTT","CGAT","CGCT","CGGT",
    "CGTT","GGAT","GGCT","GGGT","GGTT","TGAT","TGCT","TGGT","TGTT",
    "ATAA","ATCA","ATGA","ATTA","CTAA","CTCA","CTGA","CTTA","GTAA",
    "GTCA","GTGA","GTTA","TTAA","TTCA","TTGA","TTTA","ATAC","ATCC",
    "ATGC","ATTC","CTAC","CTCC","CTGC","CTTC","GTAC","GTCC","GTGC",
    "GTTC","TTAC","TTCC","TTGC","TTTC","ATAG","ATCG","ATGG","ATTG",
    "CTAG","CTCG","CTGG","CTTG","GTAG","GTCG","GTGG","GTTG","TTAG",
    "TTCG","TTGG","TTTG")

to.reorder.SBS.192.for.plotting <-
  c("TGTT", "ACAA", "GGTT", "ACCA", "CGTT", "ACGA", "AGTT", "ACTA",
    "TGGT", "CCAA", "GGGT", "CCCA", "CGGT", "CCGA", "AGGT", "CCTA",
    "TGCT", "GCAA", "GGCT", "GCCA", "CGCT", "GCGA", "AGCT", "GCTA",
    "TGAT", "TCAA", "GGAT", "TCCA", "CGAT", "TCGA", "AGAT", "TCTA",
    "TGTC", "ACAG", "GGTC", "ACCG", "CGTC", "ACGG", "AGTC", "ACTG",
    "TGGC", "CCAG", "GGGC", "CCCG", "CGGC", "CCGG", "AGGC", "CCTG",
    "TGCC", "GCAG", "GGCC", "GCCG", "CGCC", "GCGG", "AGCC", "GCTG",
    "TGAC", "TCAG", "GGAC", "TCCG", "CGAC", "TCGG", "AGAC", "TCTG",
    "TGTA", "ACAT", "GGTA", "ACCT", "CGTA", "ACGT", "AGTA", "ACTT",
    "TGGA", "CCAT", "GGGA", "CCCT", "CGGA", "CCGT", "AGGA", "CCTT",
    "TGCA", "GCAT", "GGCA", "GCCT", "CGCA", "GCGT", "AGCA", "GCTT",
    "TGAA", "TCAT", "GGAA", "TCCT", "CGAA", "TCGT", "AGAA", "TCTT",
    "TATT", "ATAA", "GATT", "ATCA", "CATT", "ATGA", "AATT", "ATTA",
    "TAGT", "CTAA", "GAGT", "CTCA", "CAGT", "CTGA", "AAGT", "CTTA",
    "TACT", "GTAA", "GACT", "GTCA", "CACT", "GTGA", "AACT", "GTTA",
    "TAAT", "TTAA", "GAAT", "TTCA", "CAAT", "TTGA", "AAAT", "TTTA",
    "TATG", "ATAC", "GATG", "ATCC", "CATG", "ATGC", "AATG", "ATTC",
    "TAGG", "CTAC", "GAGG", "CTCC", "CAGG", "CTGC", "AAGG", "CTTC",
    "TACG", "GTAC", "GACG", "GTCC", "CACG", "GTGC", "AACG", "GTTC",
    "TAAG", "TTAC", "GAAG", "TTCC", "CAAG", "TTGC", "AAAG", "TTTC",
    "TATC", "ATAG", "GATC", "ATCG", "CATC", "ATGG", "AATC", "ATTG",
    "TAGC", "CTAG", "GAGC", "CTCG", "CAGC", "CTGG", "AAGC", "CTTG",
    "TACC", "GTAG", "GACC", "GTCG", "CACC", "GTGG", "AACC", "GTTG",
    "TAAC", "TTAG", "GAAC", "TTCG", "CAAC", "TTGG", "AAAC", "TTTG"
  )

catalog.row.order.SBS.96 <-
  c("ACAA", "ACCA", "ACGA", "ACTA", "CCAA", "CCCA", "CCGA", "CCTA",
    "GCAA", "GCCA", "GCGA", "GCTA", "TCAA", "TCCA", "TCGA", "TCTA",
    "ACAG", "ACCG", "ACGG", "ACTG", "CCAG", "CCCG", "CCGG", "CCTG",
    "GCAG", "GCCG", "GCGG", "GCTG", "TCAG", "TCCG", "TCGG", "TCTG",
    "ACAT", "ACCT", "ACGT", "ACTT", "CCAT", "CCCT", "CCGT", "CCTT",
    "GCAT", "GCCT", "GCGT", "GCTT", "TCAT", "TCCT", "TCGT", "TCTT",
    "ATAA", "ATCA", "ATGA", "ATTA", "CTAA", "CTCA", "CTGA", "CTTA",
    "GTAA", "GTCA", "GTGA", "GTTA", "TTAA", "TTCA", "TTGA", "TTTA",
    "ATAC", "ATCC", "ATGC", "ATTC", "CTAC", "CTCC", "CTGC", "CTTC",
    "GTAC", "GTCC", "GTGC", "GTTC", "TTAC", "TTCC", "TTGC", "TTTC",
    "ATAG", "ATCG", "ATGG", "ATTG", "CTAG", "CTCG", "CTGG", "CTTG",
    "GTAG", "GTCG", "GTGG", "GTTG", "TTAG", "TTCG", "TTGG", "TTTG"
  )

catalog.row.order.sp.SBS.96 <- 
  c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G",
    "A[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","A[T>A]A","A[T>A]C",
    "A[T>A]G","A[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","A[T>G]A",
    "A[T>G]C","A[T>G]G","A[T>G]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T",
    "C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","C[C>T]A","C[C>T]C","C[C>T]G",
    "C[C>T]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","C[T>C]A","C[T>C]C",
    "C[T>C]G","C[T>C]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[C>A]A",
    "G[C>A]C","G[C>A]G","G[C>A]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T",
    "G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","G[T>A]A","G[T>A]C","G[T>A]G",
    "G[T>A]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","G[T>G]A","G[T>G]C",
    "G[T>G]G","G[T>G]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","T[C>G]A",
    "T[C>G]C","T[C>G]G","T[C>G]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T",
    "T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","T[T>C]A","T[T>C]C","T[T>C]G",
    "T[T>C]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")

# The first two letters are reference bases while the
# last two letters are the altered bases after mutation.
catalog.row.order.DBS.78 <-
  c("ACCA", "ACCG", "ACCT", "ACGA", "ACGG", "ACGT", "ACTA", "ACTG",
    "ACTT", "ATCA", "ATCC", "ATCG", "ATGA", "ATGC", "ATTA", "CCAA",
    "CCAG", "CCAT", "CCGA", "CCGG", "CCGT", "CCTA", "CCTG", "CCTT",
    "CGAT", "CGGC", "CGGT", "CGTA", "CGTC", "CGTT", "CTAA", "CTAC",
    "CTAG", "CTGA", "CTGC", "CTGG", "CTTA", "CTTC", "CTTG", "GCAA",
    "GCAG", "GCAT", "GCCA", "GCCG", "GCTA", "TAAT", "TACG", "TACT",
    "TAGC", "TAGG", "TAGT", "TCAA", "TCAG", "TCAT", "TCCA", "TCCG",
    "TCCT", "TCGA", "TCGG", "TCGT", "TGAA", "TGAC", "TGAT", "TGCA",
    "TGCC", "TGCT", "TGGA", "TGGC", "TGGT", "TTAA", "TTAC", "TTAG",
    "TTCA", "TTCC", "TTCG", "TTGA", "TTGC", "TTGG")

catalog.row.order.sp.DBS.78 <-
  paste0(catalog.row.headers.DBS.78$Ref, ">", catalog.row.headers.DBS.78$Var)

# There are 144 stranded DBSs: 4 X 4 sources and 3 X 3 alternates;
# 4 x 4 x 3 x 3 = 144.
catalog.row.order.DBS.144 <-
  c("AACC", "AACG", "AACT", "AAGC", "AAGG", "AAGT", "AATC", "AATG",
    "AATT", "ACCA", "ACCG", "ACCT", "ACGA", "ACGG", "ACGT", "ACTA",
    "ACTG", "ACTT", "AGCA", "AGCC", "AGCT", "AGGA", "AGGC", "AGGT",
    "AGTA", "AGTC", "AGTT", "ATCA", "ATCC", "ATCG", "ATGA", "ATGC",
    "ATGG", "ATTA", "ATTC", "ATTG", "CAAC", "CAAG", "CAAT", "CAGC",
    "CAGG", "CAGT", "CATC", "CATG", "CATT", "CCAA", "CCAG", "CCAT",
    "CCGA", "CCGG", "CCGT", "CCTA", "CCTG", "CCTT", "CGAA", "CGAC",
    "CGAT", "CGGA", "CGGC", "CGGT", "CGTA", "CGTC", "CGTT", "CTAA",
    "CTAC", "CTAG", "CTGA", "CTGC", "CTGG", "CTTA", "CTTC", "CTTG",
    "GAAC", "GAAG", "GAAT", "GACC", "GACG", "GACT", "GATC", "GATG",
    "GATT", "GCAA", "GCAG", "GCAT", "GCCA", "GCCG", "GCCT", "GCTA",
    "GCTG", "GCTT", "GGAA", "GGAC", "GGAT", "GGCA", "GGCC", "GGCT",
    "GGTA", "GGTC", "GGTT", "GTAA", "GTAC", "GTAG", "GTCA", "GTCC",
    "GTCG", "GTTA", "GTTC", "GTTG", "TAAC", "TAAG", "TAAT", "TACC",
    "TACG", "TACT", "TAGC", "TAGG", "TAGT", "TCAA", "TCAG", "TCAT",
    "TCCA", "TCCG", "TCCT", "TCGA", "TCGG", "TCGT", "TGAA", "TGAC",
    "TGAT", "TGCA", "TGCC", "TGCT", "TGGA", "TGGC", "TGGT", "TTAA",
    "TTAC", "TTAG", "TTCA", "TTCC", "TTCG", "TTGA", "TTGC", "TTGG"
  )

# Create the order for plotting transcription strand bias graph of the
# 144 DBS catalog. There are 12 dinucleotide mutation types which are not
# used in the transcription analysis because the REF and ALT bases are both
# symmetrical, hence not able to determine the strand information.
# ("ATCG", "ATGC", "ATTA", "CGAT", "CGGC", "CGTA", "GCAT", "GCCG",
# "GCTA", "TAAT", "TACG", "TAGC")
to.reorder.DBS.144.for.plotting <-
  c("GTTG", "ACCA", "GTCG", "ACCG", "GTAG", "ACCT", "GTTC", "ACGA",
    "GTCC", "ACGG", "GTAC", "ACGT", "GTTA", "ACTA", "GTCA", "ACTG",
    "GTAA", "ACTT", "ATTG", "ATCA", "ATGG", "ATCC", "ATTC", "ATGA",
    "GGTT", "CCAA", "GGCT", "CCAG", "GGAT", "CCAT", "GGTC", "CCGA",
    "GGCC", "CCGG", "GGAC", "CCGT", "GGTA", "CCTA", "GGCA", "CCTG",
    "GGAA", "CCTT", "CGAC", "CGGT", "CGGA", "CGTC", "CGAA", "CGTT",
    "AGTT", "CTAA", "AGGT", "CTAC", "AGCT", "CTAG", "AGTC", "CTGA",
    "AGGC", "CTGC", "AGCC", "CTGG", "AGTA", "CTTA", "AGGA", "CTTC",
    "AGCA", "CTTG", "GCTT", "GCAA", "GCCT", "GCAG", "GCTG", "GCCA",
    "TAAG", "TACT", "TACC", "TAGG", "TAAC", "TAGT", "GATT", "TCAA",
    "GACT", "TCAG", "GAAT", "TCAT", "GATG", "TCCA", "GACG", "TCCG",
    "GAAG", "TCCT", "GATC", "TCGA", "GACC", "TCGG", "GAAC", "TCGT",
    "CATT", "TGAA", "CAGT", "TGAC", "CAAT", "TGAT", "CATG", "TGCA",
    "CAGG", "TGCC", "CAAG", "TGCT", "CATC", "TGGA", "CAGC", "TGGC",
    "CAAC", "TGGT", "AATT", "TTAA", "AAGT", "TTAC", "AACT", "TTAG",
    "AATG", "TTCA", "AAGG", "TTCC", "AACG", "TTCG", "AATC", "TTGA",
    "AAGC", "TTGC", "AACC", "TTGG")

catalog.row.order.DBS.136 <-
  c("AACA", "AACC", "AACG", "AACT", "AATA",
    "AATC", "AATG", "AATT", "ACCA", "ACCC", "ACCG", "ACCT", "ACGA",
    "ACGC", "ACGG", "ACGT", "ACTA", "ACTC", "ACTG", "ACTT", "AGCA",
    "AGCC", "AGCG", "AGCT", "ATAA", "ATAC", "ATAG", "ATAT", "ATCA",
    "ATCC", "ATCG", "ATCT", "ATGA", "ATGC", "ATGG", "ATGT", "ATTA",
    "ATTC", "ATTG", "ATTT", "CACA", "CACC", "CACG", "CACT", "CATA",
    "CATC", "CATG", "CCCA", "CCCC", "CCCG", "CCCT", "CCGA", "CCGC",
    "CCGG", "CCTA", "CCTC", "CCTG", "CCTT", "CGCA", "CGCC", "CGCG",
    "CTAA", "CTAC", "CTAG", "CTCA", "CTCC", "CTCG", "CTCT", "CTGA",
    "CTGC", "CTGG", "CTGT", "CTTA", "CTTC", "CTTG", "CTTT", "GACA",
    "GACC", "GACG", "GACT", "GATA", "GATC", "GCCA", "GCCC", "GCCG",
    "GCCT", "GCGA", "GCGC", "GCTA", "GCTC", "GCTG", "GCTT", "GGCA",
    "GGCC", "GTAA", "GTAC", "GTCA", "GTCC", "GTCG", "GTCT", "GTGA",
    "GTGC", "GTGG", "GTGT", "GTTA", "GTTC", "GTTG", "GTTT", "TACA",
    "TACC", "TACG", "TACT", "TATA", "TCCA", "TCCC", "TCCG", "TCCT",
    "TCGA", "TCTA", "TCTC", "TCTG", "TCTT", "TGCA", "TTAA", "TTCA",
    "TTCC", "TTCG", "TTCT", "TTGA", "TTGC", "TTGG", "TTGT", "TTTA",
    "TTTC", "TTTG", "TTTT")

catalog.row.order.ID <-
  c("DEL:C:1:0", "DEL:C:1:1", "DEL:C:1:2", "DEL:C:1:3", "DEL:C:1:4",
    "DEL:C:1:5+", "DEL:T:1:0", "DEL:T:1:1", "DEL:T:1:2", "DEL:T:1:3",
    "DEL:T:1:4", "DEL:T:1:5+", "INS:C:1:0", "INS:C:1:1", "INS:C:1:2",
    "INS:C:1:3", "INS:C:1:4", "INS:C:1:5+", "INS:T:1:0", "INS:T:1:1",
    "INS:T:1:2", "INS:T:1:3", "INS:T:1:4", "INS:T:1:5+", "DEL:repeats:2:0",
    "DEL:repeats:2:1", "DEL:repeats:2:2", "DEL:repeats:2:3", "DEL:repeats:2:4",
    "DEL:repeats:2:5+", "DEL:repeats:3:0", "DEL:repeats:3:1", "DEL:repeats:3:2",
    "DEL:repeats:3:3", "DEL:repeats:3:4", "DEL:repeats:3:5+", "DEL:repeats:4:0",
    "DEL:repeats:4:1", "DEL:repeats:4:2", "DEL:repeats:4:3", "DEL:repeats:4:4",
    "DEL:repeats:4:5+", "DEL:repeats:5+:0", "DEL:repeats:5+:1", "DEL:repeats:5+:2",
    "DEL:repeats:5+:3", "DEL:repeats:5+:4", "DEL:repeats:5+:5+",
    "INS:repeats:2:0", "INS:repeats:2:1", "INS:repeats:2:2", "INS:repeats:2:3",
    "INS:repeats:2:4", "INS:repeats:2:5+", "INS:repeats:3:0", "INS:repeats:3:1",
    "INS:repeats:3:2", "INS:repeats:3:3", "INS:repeats:3:4", "INS:repeats:3:5+",
    "INS:repeats:4:0", "INS:repeats:4:1", "INS:repeats:4:2", "INS:repeats:4:3",
    "INS:repeats:4:4", "INS:repeats:4:5+", "INS:repeats:5+:0", "INS:repeats:5+:1",
    "INS:repeats:5+:2", "INS:repeats:5+:3", "INS:repeats:5+:4", "INS:repeats:5+:5+",
    "DEL:MH:2:1", "DEL:MH:3:1", "DEL:MH:3:2", "DEL:MH:4:1", "DEL:MH:4:2",
    "DEL:MH:4:3", "DEL:MH:5+:1", "DEL:MH:5+:2", "DEL:MH:5+:3", "DEL:MH:5+:4",
    "DEL:MH:5+:5+")

catalog.row.order.sp.ID.83 <-
  c("1:Del:C:0", "1:Del:C:1", "1:Del:C:2", "1:Del:C:3", "1:Del:C:4", 
    "1:Del:C:5", "1:Del:T:0", "1:Del:T:1", "1:Del:T:2", "1:Del:T:3", 
    "1:Del:T:4", "1:Del:T:5", "1:Ins:C:0", "1:Ins:C:1", "1:Ins:C:2", 
    "1:Ins:C:3", "1:Ins:C:4", "1:Ins:C:5", "1:Ins:T:0", "1:Ins:T:1", 
    "1:Ins:T:2", "1:Ins:T:3", "1:Ins:T:4", "1:Ins:T:5", "2:Del:M:1", 
    "2:Del:R:0", "2:Del:R:1", "2:Del:R:2", "2:Del:R:3", "2:Del:R:4", 
    "2:Del:R:5", "2:Ins:R:0", "2:Ins:R:1", "2:Ins:R:2", "2:Ins:R:3", 
    "2:Ins:R:4", "2:Ins:R:5", "3:Del:M:1", "3:Del:M:2", "3:Del:R:0", 
    "3:Del:R:1", "3:Del:R:2", "3:Del:R:3", "3:Del:R:4", "3:Del:R:5", 
    "3:Ins:R:0", "3:Ins:R:1", "3:Ins:R:2", "3:Ins:R:3", "3:Ins:R:4", 
    "3:Ins:R:5", "4:Del:M:1", "4:Del:M:2", "4:Del:M:3", "4:Del:R:0", 
    "4:Del:R:1", "4:Del:R:2", "4:Del:R:3", "4:Del:R:4", "4:Del:R:5", 
    "4:Ins:R:0", "4:Ins:R:1", "4:Ins:R:2", "4:Ins:R:3", "4:Ins:R:4", 
    "4:Ins:R:5", "5:Del:M:1", "5:Del:M:2", "5:Del:M:3", "5:Del:M:4", 
    "5:Del:M:5", "5:Del:R:0", "5:Del:R:1", "5:Del:R:2", "5:Del:R:3", 
    "5:Del:R:4", "5:Del:R:5", "5:Ins:R:0", "5:Ins:R:1", "5:Ins:R:2", 
    "5:Ins:R:3", "5:Ins:R:4", "5:Ins:R:5")

##############################################################
# Next are "empty row headers", data frames with the columns
# need to write catalogs to disk.
##############################################################

catalog.row.headers.SBS.96 <-
  structure(
    list(
      `Mutation type` =
        c("C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G"),
      Trinucleotide =
        c("ACA", "ACC", "ACG", "ACT",
          "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA",
          "TCC", "TCG", "TCT", "ACA", "ACC", "ACG", "ACT", "CCA", "CCC",
          "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG",
          "TCT", "ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT",
          "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT", "ATA",
          "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC",
          "GTG", "GTT", "TTA", "TTC", "TTG", "TTT", "ATA", "ATC", "ATG",
          "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT",
          "TTA", "TTC", "TTG", "TTT", "ATA", "ATC", "ATG", "ATT", "CTA",
          "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT", "TTA", "TTC",
          "TTG", "TTT")),
    class = c("data.table", "data.frame"),
    row.names = c(NA, -96L)
    # , .internal.selfref = <pointer: 0x0000000002551ef0>
  )

catalog.row.headers.SBS.96.v1 <-
  structure(
    list(
      Bef = c("A", "A", "A", "A", "C", "C", "C", "C", "G", "G", "G", "G", "T", 
              "T", "T", "T", "A", "A", "A", "A", "C", "C", "C", "C", "G", "G", 
              "G", "G", "T", "T", "T", "T", "A", "A", "A", "A", "C", "C", "C", 
              "C", "G", "G", "G", "G", "T", "T", "T", "T", "A", "A", "A", "A", 
              "C", "C", "C", "C", "G", "G", "G", "G", "T", "T", "T", "T", "A", 
              "A", "A", "A", "C", "C", "C", "C", "G", "G", "G", "G", "T", "T", 
              "T", "T", "A", "A", "A", "A", "C", "C", "C", "C", "G", "G", "G", 
              "G", "T", "T", "T", "T" ), 
      Ref = c("C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", 
              "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C",
              "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", 
              "C", "C", "C", "C", "C", "C", "C", "C", "C", "T", "T", "T", "T", 
              "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", 
              "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", 
              "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", 
              "T", "T", "T", "T", "T"), 
      After = c("A",  "C", "G", "T", "A", "C", "G", "T", "A", "C", "G", "T", "A", 
                "C", "G", "T", "A", "C", "G", "T", "A", "C", "G", "T", "A", "C", 
                "G", "T", "A", "C", "G", "T", "A", "C", "G", "T", "A", "C", "G", 
                "T", "A", "C", "G", "T", "A", "C", "G", "T", "A", "C", "G", "T", 
                "A", "C", "G", "T", "A", "C", "G", "T", "A", "C", "G", "T", "A", 
                "C", "G", "T", "A", "C", "G", "T", "A", "C", "G", "T", "A", "C", 
                "G", "T", "A", "C", "G", "T", "A", "C", "G", "T", "A", "C", "G", 
                "T", "A", "C", "G", "T"), 
      Var = c("A", "A", "A", "A",  "A", "A", "A", "A", "A", "A", "A", "A", "A", 
              "A", "A", "A", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", 
              "G", "G", "G", "G", "G", "G", "T", "T", "T", "T", "T", "T", "T", 
              "T", "T", "T", "T", "T", "T", "T", "T", "T", "A", "A", "A", "A", 
              "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "C",
              "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", 
              "C", "C", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", 
              "G", "G", "G", "G", "G")), 
    row.names = c(NA, -96L), 
    class = c("data.table",  "data.frame")
    #.internal.selfref = <pointer: 0x00000200ca341ef0>
  )

catalog.row.headers.SBS.192 <-
  structure(
    list(
      Strand =
        c("T", "T", "T", "T", "T", "T", "T",
          "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T",
          "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T",
          "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T",
          "T", "T", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U",
          "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U",
          "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U",
          "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "T", "T",
          "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T",
          "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T",
          "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T",
          "T", "T", "T", "T", "T", "T", "T", "U", "U", "U", "U", "U", "U",
          "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U",
          "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U",
          "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U", "U",
          "U", "U", "U"),
      `Mutation type` =
        c("T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G"),
      Trinucleotide =
        c("TTT",
          "GTT", "CTT", "ATT", "TTG", "GTG", "CTG", "ATG", "TTC", "GTC",
          "CTC", "ATC", "TTA", "GTA", "CTA", "ATA", "TTT", "GTT", "CTT",
          "ATT", "TTG", "GTG", "CTG", "ATG", "TTC", "GTC", "CTC", "ATC",
          "TTA", "GTA", "CTA", "ATA", "TTT", "GTT", "CTT", "ATT", "TTG",
          "GTG", "CTG", "ATG", "TTC", "GTC", "CTC", "ATC", "TTA", "GTA",
          "CTA", "ATA", "ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG",
          "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT",
          "ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA",
          "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
          "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG",
          "GCT", "TCA", "TCC", "TCG", "TCT", "TCT", "GCT", "CCT", "ACT",
          "TCG", "GCG", "CCG", "ACG", "TCC", "GCC", "CCC", "ACC", "TCA",
          "GCA", "CCA", "ACA", "TCT", "GCT", "CCT", "ACT", "TCG", "GCG",
          "CCG", "ACG", "TCC", "GCC", "CCC", "ACC", "TCA", "GCA", "CCA",
          "ACA", "TCT", "GCT", "CCT", "ACT", "TCG", "GCG", "CCG", "ACG",
          "TCC", "GCC", "CCC", "ACC", "TCA", "GCA", "CCA", "ACA", "ATA",
          "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC",
          "GTG", "GTT", "TTA", "TTC", "TTG", "TTT", "ATA", "ATC", "ATG",
          "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT",
          "TTA", "TTC", "TTG", "TTT", "ATA", "ATC", "ATG", "ATT", "CTA",
          "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT", "TTA", "TTC",
          "TTG", "TTT")),
    class = c("data.table", "data.frame"),
    row.names = c(NA, -192L)
    # .internal.selfref = <pointer: 0x0000000002551ef0>
  )

catalog.row.headers.SBS.1536 <-
  structure(
    list(
      `Mutation type` =
        c("C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A", "C>A",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G", "C>G",
          "C>G", "C>G", "C>G", "C>G", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T",
          "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "C>T", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A", "T>A",
          "T>A", "T>A", "T>A", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C",
          "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>C", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G", "T>G",
          "T>G", "T>G"),
      Pentanucleotide =
        c("AACAA", "AACAC", "AACAG",
          "AACAT", "AACCA", "AACCC", "AACCG", "AACCT", "AACGA", "AACGC",
          "AACGG", "AACGT", "AACTA", "AACTC", "AACTG", "AACTT", "ACCAA",
          "ACCAC", "ACCAG", "ACCAT", "ACCCA", "ACCCC", "ACCCG", "ACCCT",
          "ACCGA", "ACCGC", "ACCGG", "ACCGT", "ACCTA", "ACCTC", "ACCTG",
          "ACCTT", "AGCAA", "AGCAC", "AGCAG", "AGCAT", "AGCCA", "AGCCC",
          "AGCCG", "AGCCT", "AGCGA", "AGCGC", "AGCGG", "AGCGT", "AGCTA",
          "AGCTC", "AGCTG", "AGCTT", "ATCAA", "ATCAC", "ATCAG", "ATCAT",
          "ATCCA", "ATCCC", "ATCCG", "ATCCT", "ATCGA", "ATCGC", "ATCGG",
          "ATCGT", "ATCTA", "ATCTC", "ATCTG", "ATCTT", "CACAA", "CACAC",
          "CACAG", "CACAT", "CACCA", "CACCC", "CACCG", "CACCT", "CACGA",
          "CACGC", "CACGG", "CACGT", "CACTA", "CACTC", "CACTG", "CACTT",
          "CCCAA", "CCCAC", "CCCAG", "CCCAT", "CCCCA", "CCCCC", "CCCCG",
          "CCCCT", "CCCGA", "CCCGC", "CCCGG", "CCCGT", "CCCTA", "CCCTC",
          "CCCTG", "CCCTT", "CGCAA", "CGCAC", "CGCAG", "CGCAT", "CGCCA",
          "CGCCC", "CGCCG", "CGCCT", "CGCGA", "CGCGC", "CGCGG", "CGCGT",
          "CGCTA", "CGCTC", "CGCTG", "CGCTT", "CTCAA", "CTCAC", "CTCAG",
          "CTCAT", "CTCCA", "CTCCC", "CTCCG", "CTCCT", "CTCGA", "CTCGC",
          "CTCGG", "CTCGT", "CTCTA", "CTCTC", "CTCTG", "CTCTT", "GACAA",
          "GACAC", "GACAG", "GACAT", "GACCA", "GACCC", "GACCG", "GACCT",
          "GACGA", "GACGC", "GACGG", "GACGT", "GACTA", "GACTC", "GACTG",
          "GACTT", "GCCAA", "GCCAC", "GCCAG", "GCCAT", "GCCCA", "GCCCC",
          "GCCCG", "GCCCT", "GCCGA", "GCCGC", "GCCGG", "GCCGT", "GCCTA",
          "GCCTC", "GCCTG", "GCCTT", "GGCAA", "GGCAC", "GGCAG", "GGCAT",
          "GGCCA", "GGCCC", "GGCCG", "GGCCT", "GGCGA", "GGCGC", "GGCGG",
          "GGCGT", "GGCTA", "GGCTC", "GGCTG", "GGCTT", "GTCAA", "GTCAC",
          "GTCAG", "GTCAT", "GTCCA", "GTCCC", "GTCCG", "GTCCT", "GTCGA",
          "GTCGC", "GTCGG", "GTCGT", "GTCTA", "GTCTC", "GTCTG", "GTCTT",
          "TACAA", "TACAC", "TACAG", "TACAT", "TACCA", "TACCC", "TACCG",
          "TACCT", "TACGA", "TACGC", "TACGG", "TACGT", "TACTA", "TACTC",
          "TACTG", "TACTT", "TCCAA", "TCCAC", "TCCAG", "TCCAT", "TCCCA",
          "TCCCC", "TCCCG", "TCCCT", "TCCGA", "TCCGC", "TCCGG", "TCCGT",
          "TCCTA", "TCCTC", "TCCTG", "TCCTT", "TGCAA", "TGCAC", "TGCAG",
          "TGCAT", "TGCCA", "TGCCC", "TGCCG", "TGCCT", "TGCGA", "TGCGC",
          "TGCGG", "TGCGT", "TGCTA", "TGCTC", "TGCTG", "TGCTT", "TTCAA",
          "TTCAC", "TTCAG", "TTCAT", "TTCCA", "TTCCC", "TTCCG", "TTCCT",
          "TTCGA", "TTCGC", "TTCGG", "TTCGT", "TTCTA", "TTCTC", "TTCTG",
          "TTCTT", "AACAA", "AACAC", "AACAG", "AACAT", "AACCA", "AACCC",
          "AACCG", "AACCT", "AACGA", "AACGC", "AACGG", "AACGT", "AACTA",
          "AACTC", "AACTG", "AACTT", "ACCAA", "ACCAC", "ACCAG", "ACCAT",
          "ACCCA", "ACCCC", "ACCCG", "ACCCT", "ACCGA", "ACCGC", "ACCGG",
          "ACCGT", "ACCTA", "ACCTC", "ACCTG", "ACCTT", "AGCAA", "AGCAC",
          "AGCAG", "AGCAT", "AGCCA", "AGCCC", "AGCCG", "AGCCT", "AGCGA",
          "AGCGC", "AGCGG", "AGCGT", "AGCTA", "AGCTC", "AGCTG", "AGCTT",
          "ATCAA", "ATCAC", "ATCAG", "ATCAT", "ATCCA", "ATCCC", "ATCCG",
          "ATCCT", "ATCGA", "ATCGC", "ATCGG", "ATCGT", "ATCTA", "ATCTC",
          "ATCTG", "ATCTT", "CACAA", "CACAC", "CACAG", "CACAT", "CACCA",
          "CACCC", "CACCG", "CACCT", "CACGA", "CACGC", "CACGG", "CACGT",
          "CACTA", "CACTC", "CACTG", "CACTT", "CCCAA", "CCCAC", "CCCAG",
          "CCCAT", "CCCCA", "CCCCC", "CCCCG", "CCCCT", "CCCGA", "CCCGC",
          "CCCGG", "CCCGT", "CCCTA", "CCCTC", "CCCTG", "CCCTT", "CGCAA",
          "CGCAC", "CGCAG", "CGCAT", "CGCCA", "CGCCC", "CGCCG", "CGCCT",
          "CGCGA", "CGCGC", "CGCGG", "CGCGT", "CGCTA", "CGCTC", "CGCTG",
          "CGCTT", "CTCAA", "CTCAC", "CTCAG", "CTCAT", "CTCCA", "CTCCC",
          "CTCCG", "CTCCT", "CTCGA", "CTCGC", "CTCGG", "CTCGT", "CTCTA",
          "CTCTC", "CTCTG", "CTCTT", "GACAA", "GACAC", "GACAG", "GACAT",
          "GACCA", "GACCC", "GACCG", "GACCT", "GACGA", "GACGC", "GACGG",
          "GACGT", "GACTA", "GACTC", "GACTG", "GACTT", "GCCAA", "GCCAC",
          "GCCAG", "GCCAT", "GCCCA", "GCCCC", "GCCCG", "GCCCT", "GCCGA",
          "GCCGC", "GCCGG", "GCCGT", "GCCTA", "GCCTC", "GCCTG", "GCCTT",
          "GGCAA", "GGCAC", "GGCAG", "GGCAT", "GGCCA", "GGCCC", "GGCCG",
          "GGCCT", "GGCGA", "GGCGC", "GGCGG", "GGCGT", "GGCTA", "GGCTC",
          "GGCTG", "GGCTT", "GTCAA", "GTCAC", "GTCAG", "GTCAT", "GTCCA",
          "GTCCC", "GTCCG", "GTCCT", "GTCGA", "GTCGC", "GTCGG", "GTCGT",
          "GTCTA", "GTCTC", "GTCTG", "GTCTT", "TACAA", "TACAC", "TACAG",
          "TACAT", "TACCA", "TACCC", "TACCG", "TACCT", "TACGA", "TACGC",
          "TACGG", "TACGT", "TACTA", "TACTC", "TACTG", "TACTT", "TCCAA",
          "TCCAC", "TCCAG", "TCCAT", "TCCCA", "TCCCC", "TCCCG", "TCCCT",
          "TCCGA", "TCCGC", "TCCGG", "TCCGT", "TCCTA", "TCCTC", "TCCTG",
          "TCCTT", "TGCAA", "TGCAC", "TGCAG", "TGCAT", "TGCCA", "TGCCC",
          "TGCCG", "TGCCT", "TGCGA", "TGCGC", "TGCGG", "TGCGT", "TGCTA",
          "TGCTC", "TGCTG", "TGCTT", "TTCAA", "TTCAC", "TTCAG", "TTCAT",
          "TTCCA", "TTCCC", "TTCCG", "TTCCT", "TTCGA", "TTCGC", "TTCGG",
          "TTCGT", "TTCTA", "TTCTC", "TTCTG", "TTCTT", "AACAA", "AACAC",
          "AACAG", "AACAT", "AACCA", "AACCC", "AACCG", "AACCT", "AACGA",
          "AACGC", "AACGG", "AACGT", "AACTA", "AACTC", "AACTG", "AACTT",
          "ACCAA", "ACCAC", "ACCAG", "ACCAT", "ACCCA", "ACCCC", "ACCCG",
          "ACCCT", "ACCGA", "ACCGC", "ACCGG", "ACCGT", "ACCTA", "ACCTC",
          "ACCTG", "ACCTT", "AGCAA", "AGCAC", "AGCAG", "AGCAT", "AGCCA",
          "AGCCC", "AGCCG", "AGCCT", "AGCGA", "AGCGC", "AGCGG", "AGCGT",
          "AGCTA", "AGCTC", "AGCTG", "AGCTT", "ATCAA", "ATCAC", "ATCAG",
          "ATCAT", "ATCCA", "ATCCC", "ATCCG", "ATCCT", "ATCGA", "ATCGC",
          "ATCGG", "ATCGT", "ATCTA", "ATCTC", "ATCTG", "ATCTT", "CACAA",
          "CACAC", "CACAG", "CACAT", "CACCA", "CACCC", "CACCG", "CACCT",
          "CACGA", "CACGC", "CACGG", "CACGT", "CACTA", "CACTC", "CACTG",
          "CACTT", "CCCAA", "CCCAC", "CCCAG", "CCCAT", "CCCCA", "CCCCC",
          "CCCCG", "CCCCT", "CCCGA", "CCCGC", "CCCGG", "CCCGT", "CCCTA",
          "CCCTC", "CCCTG", "CCCTT", "CGCAA", "CGCAC", "CGCAG", "CGCAT",
          "CGCCA", "CGCCC", "CGCCG", "CGCCT", "CGCGA", "CGCGC", "CGCGG",
          "CGCGT", "CGCTA", "CGCTC", "CGCTG", "CGCTT", "CTCAA", "CTCAC",
          "CTCAG", "CTCAT", "CTCCA", "CTCCC", "CTCCG", "CTCCT", "CTCGA",
          "CTCGC", "CTCGG", "CTCGT", "CTCTA", "CTCTC", "CTCTG", "CTCTT",
          "GACAA", "GACAC", "GACAG", "GACAT", "GACCA", "GACCC", "GACCG",
          "GACCT", "GACGA", "GACGC", "GACGG", "GACGT", "GACTA", "GACTC",
          "GACTG", "GACTT", "GCCAA", "GCCAC", "GCCAG", "GCCAT", "GCCCA",
          "GCCCC", "GCCCG", "GCCCT", "GCCGA", "GCCGC", "GCCGG", "GCCGT",
          "GCCTA", "GCCTC", "GCCTG", "GCCTT", "GGCAA", "GGCAC", "GGCAG",
          "GGCAT", "GGCCA", "GGCCC", "GGCCG", "GGCCT", "GGCGA", "GGCGC",
          "GGCGG", "GGCGT", "GGCTA", "GGCTC", "GGCTG", "GGCTT", "GTCAA",
          "GTCAC", "GTCAG", "GTCAT", "GTCCA", "GTCCC", "GTCCG", "GTCCT",
          "GTCGA", "GTCGC", "GTCGG", "GTCGT", "GTCTA", "GTCTC", "GTCTG",
          "GTCTT", "TACAA", "TACAC", "TACAG", "TACAT", "TACCA", "TACCC",
          "TACCG", "TACCT", "TACGA", "TACGC", "TACGG", "TACGT", "TACTA",
          "TACTC", "TACTG", "TACTT", "TCCAA", "TCCAC", "TCCAG", "TCCAT",
          "TCCCA", "TCCCC", "TCCCG", "TCCCT", "TCCGA", "TCCGC", "TCCGG",
          "TCCGT", "TCCTA", "TCCTC", "TCCTG", "TCCTT", "TGCAA", "TGCAC",
          "TGCAG", "TGCAT", "TGCCA", "TGCCC", "TGCCG", "TGCCT", "TGCGA",
          "TGCGC", "TGCGG", "TGCGT", "TGCTA", "TGCTC", "TGCTG", "TGCTT",
          "TTCAA", "TTCAC", "TTCAG", "TTCAT", "TTCCA", "TTCCC", "TTCCG",
          "TTCCT", "TTCGA", "TTCGC", "TTCGG", "TTCGT", "TTCTA", "TTCTC",
          "TTCTG", "TTCTT", "AATAA", "AATAC", "AATAG", "AATAT", "AATCA",
          "AATCC", "AATCG", "AATCT", "AATGA", "AATGC", "AATGG", "AATGT",
          "AATTA", "AATTC", "AATTG", "AATTT", "ACTAA", "ACTAC", "ACTAG",
          "ACTAT", "ACTCA", "ACTCC", "ACTCG", "ACTCT", "ACTGA", "ACTGC",
          "ACTGG", "ACTGT", "ACTTA", "ACTTC", "ACTTG", "ACTTT", "AGTAA",
          "AGTAC", "AGTAG", "AGTAT", "AGTCA", "AGTCC", "AGTCG", "AGTCT",
          "AGTGA", "AGTGC", "AGTGG", "AGTGT", "AGTTA", "AGTTC", "AGTTG",
          "AGTTT", "ATTAA", "ATTAC", "ATTAG", "ATTAT", "ATTCA", "ATTCC",
          "ATTCG", "ATTCT", "ATTGA", "ATTGC", "ATTGG", "ATTGT", "ATTTA",
          "ATTTC", "ATTTG", "ATTTT", "CATAA", "CATAC", "CATAG", "CATAT",
          "CATCA", "CATCC", "CATCG", "CATCT", "CATGA", "CATGC", "CATGG",
          "CATGT", "CATTA", "CATTC", "CATTG", "CATTT", "CCTAA", "CCTAC",
          "CCTAG", "CCTAT", "CCTCA", "CCTCC", "CCTCG", "CCTCT", "CCTGA",
          "CCTGC", "CCTGG", "CCTGT", "CCTTA", "CCTTC", "CCTTG", "CCTTT",
          "CGTAA", "CGTAC", "CGTAG", "CGTAT", "CGTCA", "CGTCC", "CGTCG",
          "CGTCT", "CGTGA", "CGTGC", "CGTGG", "CGTGT", "CGTTA", "CGTTC",
          "CGTTG", "CGTTT", "CTTAA", "CTTAC", "CTTAG", "CTTAT", "CTTCA",
          "CTTCC", "CTTCG", "CTTCT", "CTTGA", "CTTGC", "CTTGG", "CTTGT",
          "CTTTA", "CTTTC", "CTTTG", "CTTTT", "GATAA", "GATAC", "GATAG",
          "GATAT", "GATCA", "GATCC", "GATCG", "GATCT", "GATGA", "GATGC",
          "GATGG", "GATGT", "GATTA", "GATTC", "GATTG", "GATTT", "GCTAA",
          "GCTAC", "GCTAG", "GCTAT", "GCTCA", "GCTCC", "GCTCG", "GCTCT",
          "GCTGA", "GCTGC", "GCTGG", "GCTGT", "GCTTA", "GCTTC", "GCTTG",
          "GCTTT", "GGTAA", "GGTAC", "GGTAG", "GGTAT", "GGTCA", "GGTCC",
          "GGTCG", "GGTCT", "GGTGA", "GGTGC", "GGTGG", "GGTGT", "GGTTA",
          "GGTTC", "GGTTG", "GGTTT", "GTTAA", "GTTAC", "GTTAG", "GTTAT",
          "GTTCA", "GTTCC", "GTTCG", "GTTCT", "GTTGA", "GTTGC", "GTTGG",
          "GTTGT", "GTTTA", "GTTTC", "GTTTG", "GTTTT", "TATAA", "TATAC",
          "TATAG", "TATAT", "TATCA", "TATCC", "TATCG", "TATCT", "TATGA",
          "TATGC", "TATGG", "TATGT", "TATTA", "TATTC", "TATTG", "TATTT",
          "TCTAA", "TCTAC", "TCTAG", "TCTAT", "TCTCA", "TCTCC", "TCTCG",
          "TCTCT", "TCTGA", "TCTGC", "TCTGG", "TCTGT", "TCTTA", "TCTTC",
          "TCTTG", "TCTTT", "TGTAA", "TGTAC", "TGTAG", "TGTAT", "TGTCA",
          "TGTCC", "TGTCG", "TGTCT", "TGTGA", "TGTGC", "TGTGG", "TGTGT",
          "TGTTA", "TGTTC", "TGTTG", "TGTTT", "TTTAA", "TTTAC", "TTTAG",
          "TTTAT", "TTTCA", "TTTCC", "TTTCG", "TTTCT", "TTTGA", "TTTGC",
          "TTTGG", "TTTGT", "TTTTA", "TTTTC", "TTTTG", "TTTTT", "AATAA",
          "AATAC", "AATAG", "AATAT", "AATCA", "AATCC", "AATCG", "AATCT",
          "AATGA", "AATGC", "AATGG", "AATGT", "AATTA", "AATTC", "AATTG",
          "AATTT", "ACTAA", "ACTAC", "ACTAG", "ACTAT", "ACTCA", "ACTCC",
          "ACTCG", "ACTCT", "ACTGA", "ACTGC", "ACTGG", "ACTGT", "ACTTA",
          "ACTTC", "ACTTG", "ACTTT", "AGTAA", "AGTAC", "AGTAG", "AGTAT",
          "AGTCA", "AGTCC", "AGTCG", "AGTCT", "AGTGA", "AGTGC", "AGTGG",
          "AGTGT", "AGTTA", "AGTTC", "AGTTG", "AGTTT", "ATTAA", "ATTAC",
          "ATTAG", "ATTAT", "ATTCA", "ATTCC", "ATTCG", "ATTCT", "ATTGA",
          "ATTGC", "ATTGG", "ATTGT", "ATTTA", "ATTTC", "ATTTG", "ATTTT",
          "CATAA", "CATAC", "CATAG", "CATAT", "CATCA", "CATCC", "CATCG",
          "CATCT", "CATGA", "CATGC", "CATGG", "CATGT", "CATTA", "CATTC",
          "CATTG", "CATTT", "CCTAA", "CCTAC", "CCTAG", "CCTAT", "CCTCA",
          "CCTCC", "CCTCG", "CCTCT", "CCTGA", "CCTGC", "CCTGG", "CCTGT",
          "CCTTA", "CCTTC", "CCTTG", "CCTTT", "CGTAA", "CGTAC", "CGTAG",
          "CGTAT", "CGTCA", "CGTCC", "CGTCG", "CGTCT", "CGTGA", "CGTGC",
          "CGTGG", "CGTGT", "CGTTA", "CGTTC", "CGTTG", "CGTTT", "CTTAA",
          "CTTAC", "CTTAG", "CTTAT", "CTTCA", "CTTCC", "CTTCG", "CTTCT",
          "CTTGA", "CTTGC", "CTTGG", "CTTGT", "CTTTA", "CTTTC", "CTTTG",
          "CTTTT", "GATAA", "GATAC", "GATAG", "GATAT", "GATCA", "GATCC",
          "GATCG", "GATCT", "GATGA", "GATGC", "GATGG", "GATGT", "GATTA",
          "GATTC", "GATTG", "GATTT", "GCTAA", "GCTAC", "GCTAG", "GCTAT",
          "GCTCA", "GCTCC", "GCTCG", "GCTCT", "GCTGA", "GCTGC", "GCTGG",
          "GCTGT", "GCTTA", "GCTTC", "GCTTG", "GCTTT", "GGTAA", "GGTAC",
          "GGTAG", "GGTAT", "GGTCA", "GGTCC", "GGTCG", "GGTCT", "GGTGA",
          "GGTGC", "GGTGG", "GGTGT", "GGTTA", "GGTTC", "GGTTG", "GGTTT",
          "GTTAA", "GTTAC", "GTTAG", "GTTAT", "GTTCA", "GTTCC", "GTTCG",
          "GTTCT", "GTTGA", "GTTGC", "GTTGG", "GTTGT", "GTTTA", "GTTTC",
          "GTTTG", "GTTTT", "TATAA", "TATAC", "TATAG", "TATAT", "TATCA",
          "TATCC", "TATCG", "TATCT", "TATGA", "TATGC", "TATGG", "TATGT",
          "TATTA", "TATTC", "TATTG", "TATTT", "TCTAA", "TCTAC", "TCTAG",
          "TCTAT", "TCTCA", "TCTCC", "TCTCG", "TCTCT", "TCTGA", "TCTGC",
          "TCTGG", "TCTGT", "TCTTA", "TCTTC", "TCTTG", "TCTTT", "TGTAA",
          "TGTAC", "TGTAG", "TGTAT", "TGTCA", "TGTCC", "TGTCG", "TGTCT",
          "TGTGA", "TGTGC", "TGTGG", "TGTGT", "TGTTA", "TGTTC", "TGTTG",
          "TGTTT", "TTTAA", "TTTAC", "TTTAG", "TTTAT", "TTTCA", "TTTCC",
          "TTTCG", "TTTCT", "TTTGA", "TTTGC", "TTTGG", "TTTGT", "TTTTA",
          "TTTTC", "TTTTG", "TTTTT", "AATAA", "AATAC", "AATAG", "AATAT",
          "AATCA", "AATCC", "AATCG", "AATCT", "AATGA", "AATGC", "AATGG",
          "AATGT", "AATTA", "AATTC", "AATTG", "AATTT", "ACTAA", "ACTAC",
          "ACTAG", "ACTAT", "ACTCA", "ACTCC", "ACTCG", "ACTCT", "ACTGA",
          "ACTGC", "ACTGG", "ACTGT", "ACTTA", "ACTTC", "ACTTG", "ACTTT",
          "AGTAA", "AGTAC", "AGTAG", "AGTAT", "AGTCA", "AGTCC", "AGTCG",
          "AGTCT", "AGTGA", "AGTGC", "AGTGG", "AGTGT", "AGTTA", "AGTTC",
          "AGTTG", "AGTTT", "ATTAA", "ATTAC", "ATTAG", "ATTAT", "ATTCA",
          "ATTCC", "ATTCG", "ATTCT", "ATTGA", "ATTGC", "ATTGG", "ATTGT",
          "ATTTA", "ATTTC", "ATTTG", "ATTTT", "CATAA", "CATAC", "CATAG",
          "CATAT", "CATCA", "CATCC", "CATCG", "CATCT", "CATGA", "CATGC",
          "CATGG", "CATGT", "CATTA", "CATTC", "CATTG", "CATTT", "CCTAA",
          "CCTAC", "CCTAG", "CCTAT", "CCTCA", "CCTCC", "CCTCG", "CCTCT",
          "CCTGA", "CCTGC", "CCTGG", "CCTGT", "CCTTA", "CCTTC", "CCTTG",
          "CCTTT", "CGTAA", "CGTAC", "CGTAG", "CGTAT", "CGTCA", "CGTCC",
          "CGTCG", "CGTCT", "CGTGA", "CGTGC", "CGTGG", "CGTGT", "CGTTA",
          "CGTTC", "CGTTG", "CGTTT", "CTTAA", "CTTAC", "CTTAG", "CTTAT",
          "CTTCA", "CTTCC", "CTTCG", "CTTCT", "CTTGA", "CTTGC", "CTTGG",
          "CTTGT", "CTTTA", "CTTTC", "CTTTG", "CTTTT", "GATAA", "GATAC",
          "GATAG", "GATAT", "GATCA", "GATCC", "GATCG", "GATCT", "GATGA",
          "GATGC", "GATGG", "GATGT", "GATTA", "GATTC", "GATTG", "GATTT",
          "GCTAA", "GCTAC", "GCTAG", "GCTAT", "GCTCA", "GCTCC", "GCTCG",
          "GCTCT", "GCTGA", "GCTGC", "GCTGG", "GCTGT", "GCTTA", "GCTTC",
          "GCTTG", "GCTTT", "GGTAA", "GGTAC", "GGTAG", "GGTAT", "GGTCA",
          "GGTCC", "GGTCG", "GGTCT", "GGTGA", "GGTGC", "GGTGG", "GGTGT",
          "GGTTA", "GGTTC", "GGTTG", "GGTTT", "GTTAA", "GTTAC", "GTTAG",
          "GTTAT", "GTTCA", "GTTCC", "GTTCG", "GTTCT", "GTTGA", "GTTGC",
          "GTTGG", "GTTGT", "GTTTA", "GTTTC", "GTTTG", "GTTTT", "TATAA",
          "TATAC", "TATAG", "TATAT", "TATCA", "TATCC", "TATCG", "TATCT",
          "TATGA", "TATGC", "TATGG", "TATGT", "TATTA", "TATTC", "TATTG",
          "TATTT", "TCTAA", "TCTAC", "TCTAG", "TCTAT", "TCTCA", "TCTCC",
          "TCTCG", "TCTCT", "TCTGA", "TCTGC", "TCTGG", "TCTGT", "TCTTA",
          "TCTTC", "TCTTG", "TCTTT", "TGTAA", "TGTAC", "TGTAG", "TGTAT",
          "TGTCA", "TGTCC", "TGTCG", "TGTCT", "TGTGA", "TGTGC", "TGTGG",
          "TGTGT", "TGTTA", "TGTTC", "TGTTG", "TGTTT", "TTTAA", "TTTAC",
          "TTTAG", "TTTAT", "TTTCA", "TTTCC", "TTTCG", "TTTCT", "TTTGA",
          "TTTGC", "TTTGG", "TTTGT", "TTTTA", "TTTTC", "TTTTG", "TTTTT"
        )),
    class = c("data.table", "data.frame"),
    row.names = c(NA, -1536L)
    #, .internal.selfref = <pointer: 0x0000000006221ef0>
  )

catalog.row.headers.DBS.78 <-
  structure(
    list(
      Ref =
        c("AC", "AC", "AC", "AC", "AC", "AC", "AC",
          "AC", "AC", "AT", "AT", "AT", "AT", "AT", "AT", "CC", "CC", "CC",
          "CC", "CC", "CC", "CC", "CC", "CC", "CG", "CG", "CG", "CG", "CG",
          "CG", "CT", "CT", "CT", "CT", "CT", "CT", "CT", "CT", "CT", "GC",
          "GC", "GC", "GC", "GC", "GC", "TA", "TA", "TA", "TA", "TA", "TA",
          "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TG", "TG",
          "TG", "TG", "TG", "TG", "TG", "TG", "TG", "TT", "TT", "TT", "TT",
          "TT", "TT", "TT", "TT", "TT"),
      Var = c("CA", "CG", "CT", "GA",
              "GG", "GT", "TA", "TG", "TT", "CA", "CC", "CG", "GA", "GC", "TA",
              "AA", "AG", "AT", "GA", "GG", "GT", "TA", "TG", "TT", "AT", "GC",
              "GT", "TA", "TC", "TT", "AA", "AC", "AG", "GA", "GC", "GG", "TA",
              "TC", "TG", "AA", "AG", "AT", "CA", "CG", "TA", "AT", "CG", "CT",
              "GC", "GG", "GT", "AA", "AG", "AT", "CA", "CG", "CT", "GA", "GG",
              "GT", "AA", "AC", "AT", "CA", "CC", "CT", "GA", "GC", "GT", "AA",
              "AC", "AG", "CA", "CC", "CG", "GA", "GC", "GG")),
    class = c("data.table", "data.frame"),
    row.names = c(NA, -78L)
    # ,  .internal.selfref = <pointer: 0x0000000006221ef0>
  )

catalog.row.headers.DBS.144 <-
  structure(
    list(
      Ref =
        c("AA", "AA", "AA", "AA", "AA", "AA", "AA",
          "AA", "AA", "AC", "AC", "AC", "AC", "AC", "AC", "AC", "AC", "AC",
          "AG", "AG", "AG", "AG", "AG", "AG", "AG", "AG", "AG", "AT", "AT",
          "AT", "AT", "AT", "AT", "AT", "AT", "AT", "CA", "CA", "CA", "CA",
          "CA", "CA", "CA", "CA", "CA", "CC", "CC", "CC", "CC", "CC", "CC",
          "CC", "CC", "CC", "CG", "CG", "CG", "CG", "CG", "CG", "CG", "CG",
          "CG", "CT", "CT", "CT", "CT", "CT", "CT", "CT", "CT", "CT", "GA",
          "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GC", "GC", "GC",
          "GC", "GC", "GC", "GC", "GC", "GC", "GG", "GG", "GG", "GG", "GG",
          "GG", "GG", "GG", "GG", "GT", "GT", "GT", "GT", "GT", "GT", "GT",
          "GT", "GT", "TA", "TA", "TA", "TA", "TA", "TA", "TA", "TA", "TA",
          "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TG", "TG",
          "TG", "TG", "TG", "TG", "TG", "TG", "TG", "TT", "TT", "TT", "TT",
          "TT", "TT", "TT", "TT", "TT"),
      Var = c("CC", "CG", "CT", "GC",
              "GG", "GT", "TC", "TG", "TT", "CA", "CG", "CT", "GA", "GG", "GT",
              "TA", "TG", "TT", "CA", "CC", "CT", "GA", "GC", "GT", "TA", "TC",
              "TT", "CA", "CC", "CG", "GA", "GC", "GG", "TA", "TC", "TG", "AC",
              "AG", "AT", "GC", "GG", "GT", "TC", "TG", "TT", "AA", "AG", "AT",
              "GA", "GG", "GT", "TA", "TG", "TT", "AA", "AC", "AT", "GA", "GC",
              "GT", "TA", "TC", "TT", "AA", "AC", "AG", "GA", "GC", "GG", "TA",
              "TC", "TG", "AC", "AG", "AT", "CC", "CG", "CT", "TC", "TG", "TT",
              "AA", "AG", "AT", "CA", "CG", "CT", "TA", "TG", "TT", "AA", "AC",
              "AT", "CA", "CC", "CT", "TA", "TC", "TT", "AA", "AC", "AG", "CA",
              "CC", "CG", "TA", "TC", "TG", "AC", "AG", "AT", "CC", "CG", "CT",
              "GC", "GG", "GT", "AA", "AG", "AT", "CA", "CG", "CT", "GA", "GG",
              "GT", "AA", "AC", "AT", "CA", "CC", "CT", "GA", "GC", "GT", "AA",
              "AC", "AG", "CA", "CC", "CG", "GA", "GC", "GG")),
    class = c("data.table", "data.frame"),
    row.names = c(NA, -144L))
# .internal.selfref = <pointer: 0x0000000002611ef0>)

catalog.row.headers.DBS.136 <-
  structure(
    list(
      Quad =
        c("AACA", "AACC", "AACG", "AACT", "AATA",
          "AATC", "AATG", "AATT", "ACCA", "ACCC", "ACCG", "ACCT", "ACGA",
          "ACGC", "ACGG", "ACGT", "ACTA", "ACTC", "ACTG", "ACTT", "AGCA",
          "AGCC", "AGCG", "AGCT", "ATAA", "ATAC", "ATAG", "ATAT", "ATCA",
          "ATCC", "ATCG", "ATCT", "ATGA", "ATGC", "ATGG", "ATGT", "ATTA",
          "ATTC", "ATTG", "ATTT", "CACA", "CACC", "CACG", "CACT", "CATA",
          "CATC", "CATG", "CCCA", "CCCC", "CCCG", "CCCT", "CCGA", "CCGC",
          "CCGG", "CCTA", "CCTC", "CCTG", "CCTT", "CGCA", "CGCC", "CGCG",
          "CTAA", "CTAC", "CTAG", "CTCA", "CTCC", "CTCG", "CTCT", "CTGA",
          "CTGC", "CTGG", "CTGT", "CTTA", "CTTC", "CTTG", "CTTT", "GACA",
          "GACC", "GACG", "GACT", "GATA", "GATC", "GCCA", "GCCC", "GCCG",
          "GCCT", "GCGA", "GCGC", "GCTA", "GCTC", "GCTG", "GCTT", "GGCA",
          "GGCC", "GTAA", "GTAC", "GTCA", "GTCC", "GTCG", "GTCT", "GTGA",
          "GTGC", "GTGG", "GTGT", "GTTA", "GTTC", "GTTG", "GTTT", "TACA",
          "TACC", "TACG", "TACT", "TATA", "TCCA", "TCCC", "TCCG", "TCCT",
          "TCGA", "TCTA", "TCTC", "TCTG", "TCTT", "TGCA", "TTAA", "TTCA",
          "TTCC", "TTCG", "TTCT", "TTGA", "TTGC", "TTGG", "TTGT", "TTTA",
          "TTTC", "TTTG", "TTTT")),
    class = c("data.table", "data.frame"),
    row.names = c(NA, -136L))
# .internal.selfref = <pointer: 0x0000000002611ef0>)

catalog.row.headers.ID <-
  structure(
    list(
      Type =
        c("DEL", "DEL", "DEL", "DEL", "DEL", "DEL",
          "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "INS", "INS", "INS",
          "INS", "INS", "INS", "INS", "INS", "INS", "INS", "INS", "INS",
          "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL",
          "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL",
          "DEL", "DEL", "DEL", "DEL", "DEL", "DEL", "INS", "INS", "INS",
          "INS", "INS", "INS", "INS", "INS", "INS", "INS", "INS", "INS",
          "INS", "INS", "INS", "INS", "INS", "INS", "INS", "INS", "INS",
          "INS", "INS", "INS", "DEL", "DEL", "DEL", "DEL", "DEL", "DEL",
          "DEL", "DEL", "DEL", "DEL", "DEL"),
      Subtype =
        c("C", "C", "C",
          "C", "C", "C", "T", "T", "T", "T", "T", "T", "C", "C", "C", "C",
          "C", "C", "T", "T", "T", "T", "T", "T", "repeats", "repeats",
          "repeats", "repeats", "repeats", "repeats", "repeats", "repeats",
          "repeats", "repeats", "repeats", "repeats", "repeats", "repeats",
          "repeats", "repeats", "repeats", "repeats", "repeats", "repeats",
          "repeats", "repeats", "repeats", "repeats", "repeats", "repeats",
          "repeats", "repeats", "repeats", "repeats", "repeats", "repeats",
          "repeats", "repeats", "repeats", "repeats", "repeats", "repeats",
          "repeats", "repeats", "repeats", "repeats", "repeats", "repeats",
          "repeats", "repeats", "repeats", "repeats", "MH", "MH", "MH",
          "MH", "MH", "MH", "MH", "MH", "MH", "MH", "MH"),
      Indel_size =
        c("1",
          "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1",
          "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2",
          "2", "2", "2", "3", "3", "3", "3", "3", "3", "4", "4", "4", "4",
          "4", "4", "5+", "5+", "5+", "5+", "5+", "5+", "2", "2", "2",
          "2", "2", "2", "3", "3", "3", "3", "3", "3", "4", "4", "4", "4",
          "4", "4", "5+", "5+", "5+", "5+", "5+", "5+", "2", "3", "3",
          "4", "4", "4", "5+", "5+", "5+", "5+", "5+"),
      Repeat_MH_size
      = c("0",
          "1", "2", "3", "4", "5+", "0", "1", "2", "3", "4", "5+", "0",
          "1", "2", "3", "4", "5+", "0", "1", "2", "3", "4", "5+", "0",
          "1", "2", "3", "4", "5+", "0", "1", "2", "3", "4", "5+", "0",
          "1", "2", "3", "4", "5+", "0", "1", "2", "3", "4", "5+", "0",
          "1", "2", "3", "4", "5+", "0", "1", "2", "3", "4", "5+", "0",
          "1", "2", "3", "4", "5+", "0", "1", "2", "3", "4", "5+", "1",
          "1", "2", "1", "2", "3", "1", "2", "3", "4", "5+")),
    class =
      c("data.table", "data.frame"), row.names = c(NA, -83L)
    #, .internal.selfref = <pointer: 0x0000000008dc1ef0>
  )
# Create a list of empty matrices
emptySBS96 <- matrix(0, nrow = 96, ncol = 0)
rownames(emptySBS96) <- catalog.row.order.SBS.96

emptySBS192 <- matrix(0, nrow = 192, ncol = 0)
rownames(emptySBS192) <- catalog.row.order.SBS.192

emptySBS1536 <- matrix(0, nrow = 1536, ncol = 0)
rownames(emptySBS1536) <- catalog.row.order.SBS.1536

emptyDBS78 <- matrix(0, nrow = 78, ncol = 0)
rownames(emptyDBS78) <- catalog.row.order.DBS.78

emptyDBS144 <- matrix(0, nrow = 144, ncol = 0)
rownames(emptyDBS144) <- catalog.row.order.DBS.144

emptyDBS136 <- matrix(0, nrow = 136, ncol = 0)
rownames(emptyDBS136) <- catalog.row.order.DBS.136

empty.cats <- list(catSBS96 = emptySBS96,
                   catSBS192 = emptySBS192,
                   catSBS1536 = emptySBS1536,
                   catDBS78 = emptyDBS78,
                   catDBS144 = emptyDBS144,
                   catDBS136 = emptyDBS136)

catalog.row.order <- list(SBS96 = catalog.row.order.SBS.96,
                          SBS192 = catalog.row.order.SBS.192,
                          SBS1536 = catalog.row.order.SBS.1536,
                          DBS78 = catalog.row.order.DBS.78,
                          DBS136 = catalog.row.order.DBS.136,
                          DBS144 = catalog.row.order.DBS.144,
                          ID = catalog.row.order.ID,
                          # NOT TESTED
                          COMPOSITE = c(catalog.row.order.SBS.1536,
                                        catalog.row.order.DBS.78,
                                        catalog.row.order.ID))

catalog.row.order.sp <- list(SBS96 = catalog.row.order.sp.SBS.96,
                             SBS1536 = catalog.row.order.sp.SBS.1536,
                             DBS78 = catalog.row.order.sp.DBS.78,
                             ID83 = catalog.row.order.sp.ID.83)

catalog.row.headers.COMPOSITE <-
  data.frame("Mutation type" = catalog.row.order[["COMPOSITE"]])
colnames(catalog.row.headers.COMPOSITE) <- "Mutation type"

#Create regex pattern for FilterWithHomopolymerMS
homopolymer.ms.regex.pattern <-
  c("A", "C", "G", "T", "(AC)", "(AG)", "(AT)", "(CA)",
    "(CG)", "(CT)", "(GA)", "(GC)", "(GT)", "(TA)", "(TC)", "(TG)")

homopolymer.ms.regex.pattern <-
  unlist(lapply(homopolymer.ms.regex.pattern, function(x) {
    return(paste(paste(rep(x, each = 5), collapse = ""), "+", sep = ""))
  }))

homopolymer.ms.regex.pattern <-
  paste(homopolymer.ms.regex.pattern, collapse = "|")
