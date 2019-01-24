#' Add sequence context to a data frame with ID (insertion/deletion) mutation records
#'
#' @param df A data frame storing mutation records of a VCF file
#'   containing only insertions and deletions. IMPORTANT: The
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
#' @keywords internal
AddSequenceID <- function(df, seq = BSgenome.Hsapiens.1000genomes.hs37d5) {

  stopifnot(nchar(df$REF) != nchar(df$ALT))
  stopifnot((df$ref == "") | (df$ALT == "")) # The representation has to be canonicalized

  # First, figure out how much sequence context is needed.
  df$pad.width <- -99

  df[nchar(df$REF) > nchar(df$ALT),
     # This is a deletion, so need sequence context of 5+ for repeats of any
     # size
     "pad.width"] <- nchar(df$REF) * 5

  df[nchar(df$REF) < nchar(df$ALT),
     # This is an insertion
     "pad.width"] <- nchar(df$ALT) * 5

  # Create a GRanges object with the needed width.
  # TODO(steve): diff is not defined; should it be in df$
  # e.g. df[   , "diff"] <- max(nchar(df$REF), nchar(df$ALT))
  Ranges <-
    as(data.frame(chrom = df$CHROM,
                  start = df$POS - df$pad.width, # 10,
                  end = df$POS + diff + df$pad.width # 10
    ),
    "GRanges")

  # Extract sequence context from the reference genome and add to df
  df <- dplyr::mutate(df,
                      seq.context = getSeq(seq, Ranges, as.character = TRUE))
  df <- dplyr::mutate(df, seq.pad.width = df$pad.width)
  return(df)
}

#' @title Return the number of repeat units in which a deletion
#' is embedded.
#'
#' @param context A string that embeds \code{rep.unit.seq} at position
#'  \code{pos}
#'
#' @param rep.unit.seq A substring of \code{context} at \code{pos}
#'  to \code{pos + nchar(rep.unit.seq) - 1}, which is the repeat
#'   unit sequence.
#'
#' @param pos The position of \code{rep.unit.seq} in \code{context}.
#'
#' @return The number of repeat units in which \code{rep.unit.seq} is
#' embedded, not including
#' the input rep.unit.seq in the count.
#'
#' @details
#'
#' For example \code{FindMaxRepeatDel("xyaczt", "ac", 3)}
#' returns 0.
#'
#' If
#' \code{substr(context, pos, pos + nchar(rep.unit.seq) - 1) != rep.unit.seq}
#'  then stop.
#'
#' If this functions returns 0, then it is necessary to
#'   look for microhomology.
#'
#' @keywords internal
FindMaxRepeatDel <- function(context, rep.unit.seq, pos) {
  n <- nchar(rep.unit.seq)
  stopifnot(substring(context, pos, pos + n - 1) == rep.unit.seq)

  # Look left
  i <- pos - n
  left.count <- 0
  while (i > 0) {
    if (substr(context, i, i + n - 1) == rep.unit.seq) {
      left.count <- left.count + 1
      i <- i - n
    } else break
  }

  # Look right
  right.count <- 0
  i <- pos + n
  tot.len <- nchar(context)
  while ((i + n - 1) <= tot.len) {
    if (substr(context, i, i + n -1) == rep.unit.seq) {
      right.count <- right.count + 1
      i <- i + n
    } else break
  }

  # IMPORTANT - in the catalog version of the the mutation class we exclude the
  # deleted unit from the count, but in plotting we include it.
  return(left.count + right.count)
}

#' @title Return the length of microhomology at a deletion
#'
#' TODO(steve):not finished
#'
#' Microhomology can be alligned in multiple equivalent ways.
#' Example:
#'
#' GGCTAGTT aligned to GGCTAGAACTAGTT
#'
#' \preformatted{
#' GGCTAGAACTAGTT
#' GG------CTAGTT GGCTAGTT GG[CTAGAA]CTAGTT
#'                            ----   ----
#'
#' GGC------TAGTT GGCTAGTT GGC[TAGAAC]TAGTT
#'                           * ---  * ---
#'
#' GGCT------AGTT GGCTAGTT GGCT[AGAACT]AGTT
#'                           ** --  ** --
#'
#' GGCTA------GTT GGCTAGTT GGCTA[GAACTA]GTT
#'                           *** -  *** -
#'
#' GGCTAG------TT GGCTAGTT GGCTAG[AACTAG]TT
#'                           ****   ****
#' GAC------TAGTT GACTAGTT GAC[TAGAAC]TAGTT
#'                          ** --- ** ---
#' }
#'
#'
#' \preformatted{
#' GACTAGCTAGTT
#' GACTA----GTT GACTAGTT GACTA[GCTA]GTT
#'                         *** -*** -
#'
#' ** This is really repeat TODO(steve): add check in code
#' GACTAG----TT GACTAGTT GACTAG[CTAG]TT
#'                         **** ----
#'
#' GACT----AGTT GACTAGTT GACT[AGCT]AGTT
#'                         ** --** --
#'
#'
#' }
#'
#' All the same pairs of sequence, aligned 5 different ways.
#' 4 bp of microhomology.
#'
#' Need to find:
#'
#' (1) The maxium match of undeleted sequence on left that is
#' identical to the right end of deleted sequence, and
#'
#' The microhomology sequence is the concatenation of items
#' (1) and (2).
#'
#' @param context The deleted sequence plus ample surround
#'   sequence on each side (at least as long as \code{del.sequence}).
#'
#' @param deleted.seq The deleted sequence in \code{context}.
#' #'
#' @param pos The position of \code{del.sequence} in \code{context}.
#'
#' @param trace If > 0, cat various messages.
#'
#' @return The length of the maxium microhomology of \code{del.sequence}
#'   in \code{context}.
#'
#' @keywords internal
#'
FindDelMH <- function(context, deleted.seq, pos, trace = 0) {
  n <- nchar(deleted.seq)

  stopifnot(substr(context, pos, pos + n - 1) == deleted.seq)
  # The context on the left has to be longer then deleted.seq
  stopifnot((pos - 1) > n)
  # The context on the right as to be longer than deleted.seq
  stopifnot(nchar(context) - (pos + n - 1) > n)

  ds <- unlist(strsplit(deleted.seq, ""))

  # Look for microhomology to the left in context.
  left.context <- substr(context, pos - n, pos - 1)
  left <- unlist(strsplit(x = left.context, ""))
  for (i in n:1) {
    if (ds[i] != left[i]) break
    if (i == 1) stop("Thre is a repeat to the left of ", deleted.seq)
  }
  left.len <- n - i
  if (trace > 0 ) {
    cat("Left break", i, "\nleft.len =", left.len, "\n")
  }

  # Look for microhomology to the right in context.
  right.context <- substr(context, pos + n, pos + 2 * n)
  right <- unlist(strsplit(x = right.context, ""))
  for (i2 in 1:n) {
    if (ds[i2] != right[i2]) break
    if (i2 == n) stop("There is repeat to the right of ", deleted.seq)
  }
  right.len <- i2 - 1
  if (trace > 0) {
    cat("Right break", i2, "\nright.len =", right.len, "\n")
    cat(paste0(left.context, "[",
               deleted.seq, "]",
               right.context, "\n"))
    # left.context and right.context are the same length as deleted.seq

  }
  return (left.len + right.len)
}

if (FALSE) {

  # GAGAGG[CTAGAA]CTAGTT
  #        ----   ----
  FindDelMH("GGAGAGGCTAGAACTAGTTAAAAA", "CTAGAA", 8, trace = 1)

  # GAGAGGC[TAGAAC]TAGTT
  #       * ---  * ---
  FindDelMH("GGAGAGGCTAGAACTAGTTAAAAA", "TAGAAC", 9, trace = 1)


  # TGACTA[GCTA]GTTAA
  #    *** -*** -
  FindDelMH("TGACTAGCTAGTTAA", "GCTA", 7)

  # AGATA[GATA]CCCCA
  #  **** ----
  FindDelMH("AGATAGATACCCCA", "GATA", 6)

  # ACCCCC[GATA]GATACCCCA
  #        **** ----
  FindDelMH("ACCCCCGATAGATACCCCA", "GATA", 7)




  # AAGATA[GATAG]CCCCAA
  #   **** ----
  FindDelMH("AAGATAGATAGCCCCAA", "GATAG", 7)

  # AAGATA[GGATA]CCCCAAA
  #   ****  ----
  FindDelMH("AAGATAGGATACCCCAAA", "GGATA", 7)


  }

#' @title Return the number of repeat units in which an insertion
#' is embedded.
#'
#' @param context A string into which \code{rep.unit.seq} was
#'  inseted at position \code{pos}
#'
#' @param rep.unit.seq The inserted sequence and potention repeat unit
#'
#' @param pos \code{rep.unit.seq} is understood to be inserted between
#'   positions \code{pos} amd \code{pos + 1}.
#'
#' @return If same sequence as \code{rep.unit.seq} occurs ending at
#'   \code{pos} or starting at \code{pos + 1} then the number of
#'   repeat units before the insertion, otherwise 0.
#'
#'
#' @details
#' For example
#'
#' \preformatted{
#'
#' rep.unit.seq = ac
#' pos = 2
#' context = xyaczt
#' return 1
#'
#' rep.unit.seq = ac
#' pos = 4
#' context = xyaczt
#' return 1
#'
#' rep.unit.seq = cgct
#' pos = 2
#' rep.unit.seq = at
#' return 0
#'
#' context = gacacacacg
#' rep.unit.seq = ac
#' pos = any of 1, 3, 5, 7, 9
#' return 4
#' }
#'
#' If \code{substr(context, pos, pos + nchar(rep.unit.seq) - 1) != rep.unit.seq} then stop.
#'
#' @keywords internal
#'
FindMaxRepeatIns <- function(context, rep.unit.seq, pos) {

  n <- nchar(rep.unit.seq)

  # If rep.unit.seq is in context adjacent to pos, it might start at pos + 1 -
  # len(rep.unit.seq), so look left
  left.count <- 0
  p <- pos + 1 - n
  while (p > 0) {
    if (substring(context, p, p + n - 1) == rep.unit.seq) {
      left.count <- left.count + 1
    } else break
    p <- p - n
  }

  # If rep.unit.seq is repeated in context it might start at pos + 1,
  # so look right.
  right.count <- 0
  p <- pos + 1
  tot.len <- nchar(context)
  while ((p + n - 1) <= tot.len) {
    if (substr(context, p, p + n - 1) == rep.unit.seq) {
      right.count <- right.count + 1
    } else break
    p <- p + n
  }

  return(left.count + right.count)
}

#' Canonicalize1DEL
#'
#' @param ref TODO
#' @param alt TODO
#' @param context TODO
#' @keywords internal
#' @return TODO
#' @export
Canonicalize1DEL <- function(ref, alt, context) {
  if ("-" == alt) {
    alt <- ""
  }
  # TODO(steve): insert code to deal with the case that ref and alt share a 1-base prefix
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
#' @keywords internal
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
#' @keywords internal
#' @return TODO
#' @export
#'
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
#' @keywords internal
#' @return TODO

CanonicalizeID <- function(ref, alt, context) {
  ret <- mapply(Canonicalize1ID, ref, alt, context)
  return(ret)
}

#' Create an indel (ID) mutation catalog for *one* sample from a Variant Call Format (VCF)
#' file
#'
#' @param ID.vcf An in-memory VCF as a data.frame annotated by the AddSequence
#'   and AddTranscript functions. It must only contain indels and must *not*
#'   contain SBS (single base substituions), DBS, or triplet base substituions
#'   etc.
#'
#'   * Sequence must already have been added to ID.vcf
#'
#'   One design decision for variant callers is the representation of "complex
#'   indels", e.g. mutations e.g. CAT > GC. Some callers represent this as C>G,
#'   A>C, and T>_. Others might represent it as CAT > CG. Multiple issues can
#'   arise. In PCAWG, overlapping indel/SBS calls from different callers were
#'   included in the indel VCFs.
#'
#' @param SBS.vcf An in-memory VCF as a data frame. Because we have to work with
#'   some PCAWG data, we will look for neigboring indels and indels adjoining
#'   SBS. That means this functions takes an SBS VCF and an ID VCF from the same
#'   sample.
#'
#' @return A list with two elemsents:
#'   ID.cat:   A 1-column matrix containing the mutation catalog information.
#'   problems: Locations of neighboring indels or indels neighboring SBS.
#'             In the future we might handle these depending on what we
#'             find in the indel calls from different variant callers.
#' TODO(steve) Is problems implemented?
CreateOneColIDCatalog <- function(ID.vcf, SBS.vcf) {
  # TODO(steve): more checking of the ID VCF here
  stopifnot(nchar(SBS.vcf$ALT) == 1)
  stopifnot(nchar(SBS.vcf$REF) == 1)
  stopifnot("seq.21context" %in% names(ID.vcf))

  canon.ID <- CanonicalizeID(ID.vcf$REF, ID.vcf$ALT, ID.vcf$seq.21context)

  # Create the ID catalog matrix
  tab.ID <- table(canon.ID)

  row.order <- data.table(rn = .catalog.row.order.ID) # TODO(steve): can reduce use of data.table?

  ID.dt <- as.data.table(tab.ID)
  # ID.dt has two columns, names cannon.dt (from the table() function
  # and N (the count)

  ID.dt2 <-
    merge(row.order, ID.dt, by.x="rn", by.y="canon.ID", all = TRUE)
  ID.dt2[ is.na(N) , N := 0]
  stopifnot(unlist(ID.dt2$rn) == .catalog.row.order.ID)

  ID.mat <- as.matrix(ID.dt2[ , 2])
  rownames(ID.mat) <- ID.dt2$rn
  return(ID.mat)
}


if (FALSE) {
  TestFindMaxRepeatDel <- function() {
    FindMaxRepeatDel("abcabc", "abc", 1)
    FindMaxRepeatDel("abc", "abc", 1)
    FindMaxRepeatDel("abcabc", "abc", 4)
    FindMaxRepeatDel("abcabcabc", "abc", 4)
    FindMaxRepeatDel("abcxyx", "abc", 4)
    FindMaxRepeatDel("xyzxyz", "abc", 4)
    FindMaxRepeatDel("xyzabcxyz", "abc", 4)
    FindMaxRepeatDel("xxxaaayyy", "a", 4)
    FindMaxRepeatDel("xxxaaayyy", "a", 3)
    FindMaxRepeatDel("xxxaaayyy", "a", 2)
    FindMaxRepeatDel("xxxaaayyy", "a", 7)
    FindMaxRepeatDel("xxxaaayyy", "a", 8)

    TestFindMaxRepeatIns <- function() {
      stop() # totally broken
      FindMaxRepeatIns("abcabc", "abc", 1)
      FindMaxRepeatIns("abc", "abc", 1)
      FindMaxRepeatIns("abcabc", "abc", 4)
      FindMaxRepeatIns("abcabcabc", "abc", 4)
      FindMaxRepeatIns("abcxyx", "abc", 4)
      FindMaxRepeatIns("xyzxyz", "abc", 4)
      FindMaxRepeatIns("xyzabcxyz", "abc", 4)
      FindMaxRepeatIns("xxxaaayyy", "a", 4)
      FindMaxRepeatIns("xxxaaayyy", "a", 3)
      FindMaxRepeatIns("xxxaaayyy", "a", 2)
      FindMaxRepeatIns("xxxaaayyy", "a", 7)
      FindMaxRepeatIns("xxxaaayyy", "a", 8)
    }

  }

}
