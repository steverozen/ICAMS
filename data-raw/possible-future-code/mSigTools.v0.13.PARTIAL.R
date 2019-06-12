# Source this file -- creates some gobal varaibles.
#
# Partial version of misigtools that can be sourced without
# causeing too many problems.
#
# mSigTools
#
# v0.11
#
# An alpha version
#
# 2019 06 12
#
# Copyright 2017 by Alvin Wei Tian Ng, Steven G. Rozen
#
# The code is released under GPL-3
# https://www.gnu.org/licenses/gpl-3.0.en.html

# Dependencies
library(stringi) # Needed?
library('lsa') # for cosine()
library(ICAMS)

# Not needed, ICAMS revc works the same
# revc <- function(seq) {
#
#  rq <- stri_reverse(seq)
#  rq1 <- stri_trans_char(rq, 'ACGT', '1234')
#  # We have to do this in two steps because
#  # stri_trans_char computes something like the transitive
#  # closure of the substitutions.
# stri_trans_char(rq1, '1234', 'TGCA')
# }

xtract.col <- function(spec, sample.name.or.index) {
  tmp <- as.matrix(spec[ , sample.name.or.index])
  if (ncol(tmp) > 1) return(tmp)
  rownames(tmp) <- rownames(spec)
  colnames(tmp) <-
    ifelse(mode(sample.name.or.index) == 'numeric',
           colnames(spec)[sample.name.or.index],
           sample.name.or.index)
  tmp
}

######################################
# Funtions for reading catalogs, etc.
######################################

.canonical.96.row.order <-
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

.canonical.192.row.order <-
  
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
    "TTCG","TTGG","TTTG"
    
  )

# Read 96-channel spectrum or signatures in Ludmil format
# I thinks this function is the same as ICAMS::ReadCatalog
# I.e. Ludmil's format is the same as that used in 
# ICAMS

read.96.ludmil.format <- function(path) {
  cos <- read.table(path,
                    stringsAsFactors = F,
                    as.is = T,
                    header = T,
                    sep=',', 
                    check.names = F)
  stopifnot(nrow(cos)==96)
  ref.gt.var <- cos[ ,1]
  before.ref.after <- cos[ ,2]
  var <- substring(ref.gt.var, 3, 3)
  tmp <- paste(before.ref.after, var, sep='')
  rownames(cos) <- tmp
  out <- cos[ ,-(1:2), drop=F]
  out <- as.matrix(out)
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out[.canonical.96.row.order, drop=F,]
}

### Read 96-channel spectrum or signatures in Duke-NUS format
### This is a tab separated format with mutations information in
### the columns Before, Ref, After, Var

read.96.duke.nus.format <- function(path) {
  cos <- read.table(path,
                    stringsAsFactors = F,
                    as.is = T,
                    header = T,
                    sep='\t', 
                    check.names = F)
  stopifnot(nrow(cos)==96)
  # Not sure yet what will the form of the row labels
  # But Somatic.Mutation.Type has all the necessary information,
  # after which the first 3 columns are not needed
  tmp <- paste(cos$Before, cos$Ref, cos$After, cos$Var, sep='')
  rownames(cos) <- tmp
  out <- as.matrix(cos[ ,-(1:4), drop=F])
  out[.canonical.96.row.order, drop=F,]
}

### Read 192 channel spectra in Duke-NUS format
read.and.prep.192.duke.nus.catalog <- function(path) {
  df <- read.table(path,
                       stringsAsFactors = F,
                       as.is=T,
                       header=T, 
                   check.names = F)
  # Careful, df has 192 row
  stopifnot(nrow(df)==192)
  aa.must.complement <- df$Ref %in% c('A', 'G')

  xc <- df[aa.must.complement, ]
  rn1 <- revc(paste(xc$Before, xc$Ref, xc$After, sep = ''))
  rownames(xc) <- paste(rn1, revc(xc$Var), sep='')

  x1 <- df[!aa.must.complement, ]
  rownames(x1) <- paste(x1$Before, x1$Ref, x1$After, x1$Var, sep='')
  sort.string <- paste(x1$Ref, x1$Var, x1$Before, x1$After, sep = '')
  df.nocomp <- x1[order(sort.string), ]

  df.comp <- xc[rownames(df.nocomp), ]

  # For debugging
  # tmp.out <- cbind(df.comp[ ,1:4], df.nocomp[ , 1:4])

  # We no longer need the first 4 columns
  df.comp <- df.comp[ , -(1:4)]
  df.nocomp <- df.nocomp[ , -(1:4)]
  out <- df.comp + df.nocomp
  out <- as.matrix(out)
  out.order <- order(margin.table(out, 2), colnames(out), decreasing = T)
  out2 <- out[ , out.order]
  stopifnot(rownames(out2) == .canonical.96.row.order)
  rownames(df) <- paste(df$Before, df$Ref, df$After, df$Var, sep='')

  list(channel.96=out2, channel.192=df[.canonical.192.row.order,])
}

