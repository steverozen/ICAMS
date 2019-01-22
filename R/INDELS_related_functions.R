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
AddSequenceID <- function(df, seq = BSgenome.Hsapiens.1000genomes.hs37d5) {

  # Create a GRanges object with the needed width.

  stopifnot(nchar(df$REF) != nchar(dr$ALT))
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

  # Not sure this will work if the pad.widths are different.
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

#' Return the number of repeat units in which a deletion
#' is embedded. TODO(Steve): check this statement; what
#' if there is no repeat?
#'
#' e.g. rep.unit.seq = ac
#' pos = 3
#' context = xyaczt
#' pos         ^
#' Return 0
#'
#' If \code{substr(context, pos, pos + nchar(rep.unit.seq) - 1) != rep.unit.seq} then stop.
#'
#' @param context A string that embeds \code{rep.unit.seq} at position \code{pos}
#' @param rep.unit.seq A substring of \code{context} at \code{pos}
#'  to \code{pos + nchar(rep.unit.seq) - 1}, which is the repeat unit sequence.
#' @param pos The position of \code{rep.unit.seq}.
#'
#' @return The number of repeat units in which \code{rep.unit.seq} is embedded, not include
#'   the input rep.unit.seq in the count.
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

  right.count <- 0
  i <- pos + n
  tot.len <- nchar(context)
  while ((i + n - 1) <= tot.len) {
    if (substr(context, i, i + n -1) == rep.unit.seq) {
      right.count <- right.count + 1
      i <- i + n
    } else break
  }

  # IMPORTANT - in the catalog we exclude deleted unit from the count.
  return(left.count + right.count)
}

#' FindDelMH TODO(steve):not finished
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

#' FindMaxRepeatIns TODO(steve):finish this
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
#' @export
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
