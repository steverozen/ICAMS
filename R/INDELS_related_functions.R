#' Add sequence context to a data frame with mutation records
#'
#' @param df A data frame storing mutation records of a VCF file. IMPORTANT: The
#'   representation of indels in df must have been canonicalized, so that
#'   context bases (which are added by some indel callers) are placed in a
#'   column "Left.context.base" and so that, for deletions, ALT is the empty
#'   string, and, for insertions, REF is the empty string.
#' @param seq A particular reference genome.
#' @importFrom methods as
#' @import BSgenome.Hsapiens.1000genomes.hs37d5
#' @return A data frame with 2 new columns added to the input data frame. One
#'   column contains sequence context information and the other column contains
#'   the length of the "context" string to the left of the site of the variant.
#' @export
AddSequenceID <- function(df, seq = BSgenome.Hsapiens.1000genomes.hs37d5) {

  # Create a GRanges object with the needed width.

  stopifnot(nchar(df$REF) != nchar(dr$ALT))

  # First, figure out how much sequence context is needed.
  df$pad.width <- -99

  df[nchar(df$REF) > nchar(df$ALT),
     # This is a deletion
     "pad.width"] <- nchar(df$REF) * 5
  # Need sequence context of 5+ for repeats of any size

  df[nchar(df$REF) < nchar(df$ALT),
     # This is an insertion
     "pad.width"] <- nchar(df$ALT) * 5

  # Not sure this will work if the pad.widths are different.
  Ranges <-
    as(data.frame(chrom = df$CHROM,
                  start = df$POS - df$pad.width, # 10,
                  end = df$POS + diff + df$pad.width # 10
    ),
    "GRanges")

  # Extract sequence context from the reference genome
  df <- dplyr::mutate(df,
                      seq.context = getSeq(seq, Ranges, as.character = TRUE))
  df <- dplyr::mutate(df, seq.pad.width = df$pad.width)
  return(df)
}

#' FindMaxRepeatDel
#'
#' q is a substring of context at pos to pos + len(q) - 1
#'
#'
#' e.g. q = ac
#' pos = 3
#' context = xyaczt
#' pos         ^
#'
#' for deletion, if substr(context, pos, pos + len(q) - 1) != q
#' there is an error.
#' @param context TODO
#' @param q TODO
#' @param pos TODO
#'
#' @return TODO
#' @export
FindMaxRepeatDel <- function(context, q, pos) {
  stopifnot(substring(context, pos, pos + n - 1) == q)

  # Look left
  i <- pos - n
  left.count <- 0
  while (i > 0) {
    if (substr(context, i, i + n - 1) == q) {
      left.count <- left.count + 1
      i <- i - n
    } else break
  }

  right.count <- 0
  i <- pos + n
  tot.len <- nchar(context)
  while ((i + n - 1) <= tot.len) {
    if (substr(context, i, i + n -1) == q) {
      right.count <- right.count + 1
      i <- i + n
    } else break
  }

  return(left.count + 1 + right.count)
}

#' FindDelMH
#'
#' Microhomology can be alligned in multiple equivalent ways.
#' Example:
#'
#' GGCTAGTT aligned to
#'
#' GGCTAGAACTAGTT
#' GG------CTAGTT GGCTAGTT GG[CTAGAA]CTAGTT
#'                            ----   ----
#' GGC------TAGTT GGCTAGTT GGC[TAGAAC]TAGTT
#'                           * ---  * ---
#' GGCT------AGTT GGCTAGTT
#' GGCTA------GTT GGCTAGTT
#' GGCTAG------TT GGCTAGTT
#'
#' All the same pairs of sequence, aligned 5 different ways.
#' 4 bp of microhomology.
#'
#' Need to find:
#'
#' (1) The maxium match of undeleted sequence on left that is
#' identical to the right end of deleted sequence, and
#'
#' (2) The maxium match of undeleted sequence on right that is
#' identical to the left end of deleted sequence.
#'
#' The microhomology sequence is the concatenation of items
#' (1) and (2).
#'
#' @param context TODO
#' @param q TODO
#' @param pos TODO
#'
#' @return TODO
#' @export
FindDelMH <- function(context, q, pos) {
  stopifnot(substr(context, pos, pos + nchar(q) - 1) == q)
  stopifnot(pos > nchar(q) + 1) # The context on the left has to be
  # longer than q

  i <- 0
  while (substr(context, pos - (i + 1), pos - 1) ==
         substr(q, nchar(q) - i, nchar(q))) {
    i <- i + 1
  }
  left.len <- i

  i <- 0
  right.context.start <- pos + nchar(q)
  while (substr(q, 1, 1 + i ) ==
         substr(context, right.context.start, right.context.start + i)) {
    i <- i + 1
  }
  return(left.len + i)
}

#' FindMaxRepeatIns
#'
#' If q is an insertion into context between pos and pos+1
#' if q is repeated in context it might start at pos+1:
#'
#' e.g. q = ac
#' pos = 4
#' context = abxyac
#'  pos         ^
#'  start        ^
#'
#' or q might start at pos + 1 - len(q)
#'
#' e.g. q = ac
#' pos = 4
#' context = xyaczz
#'  pos         ^
#'  start      ^
#'
#' @param context TODO
#' @param q TODO
#' @param pos TODO
#'
#' @return TODO
#' @export
FindMaxRepeatIns <- function(context, q, pos) {
  stop() # this code is complete broken
  n <- nchar(q)

  extra.count <- 0
  if (substring(context, pos, pos + n - 1) == q) {
    extra.count <- 1
  }

  # Look left
  i <- pos - n
  left.count <- 0
  while (i > 0) {
    if (substr(context, i, i + n - 1) == q) {
      left.count <- left.count + 1
      i <- i - n
    } else break
  }

  right.count <- 0
  i <- pos + n
  tot.len <- nchar(context)
  while ((i + n - 1) <= tot.len) {
    if (substr(context, i, i + n -1) == q) {
      right.count <- right.count + 1
      i <- i + n
    } else break
  }

  return(left.count + extra.count + right.count)
}

#' Canonicalize1DEL
#'
#' @param ref TODO
#' @param alt TODO
#' @param context TODO
#'
#' @return TODO
#' @export
Canonicalize1DEL <- function(ref, alt, context) {
  if ("-" == alt) {
    alt <- ""
  }
  # TODO(steve): insert code to deal with the case that ref and alt shae a 1-base prefix
  if (nchar(alt) > 0) {
    cat("possible complex indel:", ref, alt, context, "\n")
    stop()
  }
  del.len <- nchar(ref)
  stopifnot(substr(context, 11, del.len) == ref)
  # Is the deletion involved in a repeat?
  rep.count <- FindMaxRepeatDel(context, ref, 11)
}

#' Canonicalize1INS
#'
#' @param ref TODO
#' @param alt TODO
#' @param context TODO
#'
#' @return TODO
#' @export
Canonicalize1INS <- function(ref, alt, context) {
  stop() # not done
}

#' Canonicalize1ID
#'
#' @param ref TODO
#' @param alt TODO
#' @param context TODO
#'
#' @return TODO
#' @export
Canonicalize1ID <- function(ref, alt, context) {
  if (nchar(alt) < nchar(ref) || "-" == alt) {
    # A deletion
    Canonicalize1INS(ref, alt, context)
  } else if (nchar(alt) > nchar(ref) ||  "-" == ref) {
    # An insertion
    Canonicalize1INS(ref, alt, context)
  } else {
    cat("Non-insertion / non-deletion found:", ref, alt, context, "\n")
    stop()
  }
}

#' CanonicalizeID
#'
#' @param ref TODO
#' @param alt TODO
#' @param context TODO
#'
#' @return TODO
#' @export
CanonicalizeID <- function(ref, alt, context) {
  ret <- mapply(Canonicalize1ID, ref, alt, context)
  return(ret)
}

#' Create mutation catalog for *one* sample from a Variant Call Format (VCF)
#' file
#'
#' @param ID.vcf An in-memory VCF as a data.frame annotated by the AddSequence
#'   and AddTranscript functions. It must only contain indels and must *not*
#'   contain SBS (single base substituions), DBS or triplet base substituions
#'   etc.
#'
#'   * Sequence must already have been added to ID.vcf
#'
#'   One design decision for variant callers is the representation of "complex
#'   indels", e.g. mutations e.g. CAT > GC. Some callers represent this as C>G,
#'   A>C, and T>_. Others might represent it as CAT > CG. Multiple issues can
#'   arise. In PCAWG, overlapping indel/SBS calls from different callers were
#'   included in the indel VCFs.
#' @param SBS.vcf An in-memory VCF as a data frame. Because we have to work with
#'   some PCAWG data, we will look for neigboring indels and indels adjoining
#'   SBS. That means this functions takes an SBS VCF and an ID VCF from the same
#'   sample.
#'
#' @return Returns a list with two elemsents:
#'   ID.cat:   A matrix containing the mutation catalog information.
#'   Problems: Locations of neighboring indels or indels neighboring SBS.
#'             In the future we might handle these depending on what we
#'             find in the indel calls from different variant callers.
#' @export
CreateOneColIDCatalog <- function(ID.vcf, SBS.vcf) {
  # TODO(steve): more checking of the ID VCF here
  stopifnot(nchar(SBS.vcf$ALT) == 1)
  stopifnot(nchar(SBS.vcf$REF) == 1)
  stopifnot("seq.21context" %in% names(ID$vcf))


  # START here

  canon.ID <- CanonicalizeID(vcf$REF, vcf$ALT, vcf$seq.21context)




  # Create the DBS catalog matrix
  tab.DBS <- table(canon.DBS)

  row.order <- data.table(rn = .catalog.row.order.DBS)
  DBS.dt <- as.data.table(tab.DBS)
  # DBS.dt has two columns, names cannon.dt (from the table() function
  # and N (the count)

  DBS.dt2 <-
    merge(row.order, DBS.dt, by.x="rn", by.y="canon.DBS", all = TRUE)
  DBS.dt2[ is.na(N) , N := 0]
  stopifnot(unlist(DBS.dt2$rn) == .catalog.row.order.DBS)

  DBS.mat <- as.matrix(DBS.dt2[ , 2])
  rownames(DBS.mat) <- DBS.dt2$rn
  return(DBS.mat)
}
