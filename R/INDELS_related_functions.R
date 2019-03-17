#' @title Add sequence context to a data frame with ID (insertion/deletion) mutation records,
#'  and confirm that they match the given reference genome.
#'
#' @param df A data frame storing mutation records of a VCF file
#'   containing only insertions and deletions. This function expects that
#'   there is a "context base" to the left, for example REF = ACG, ALT = A
#'  (deletion of CG) or REF = A, ALT = ACC (insertion of CC).
#'
#' @param genome A genome argument as described in \code{\link{ICAMS}}.
#'
#' @param flag.mismatches If > 0, if there are mismatches to references, print
#' out the first \code{flag.mismatches} rows and continue.  Otherwise \code{stop}.
#'
#' @importFrom methods as
#'
#' @importFrom BSgenome getSeq seqnames
#'
#' @import BSgenome.Hsapiens.1000genomes.hs37d5
#'
#' @import BSgenome.Hsapiens.UCSC.hg38
#'
#' @return A data frame with 2 new columns added to the input data frame:
#' \enumerate{
#'  \item \code{seq.context} The sequence embedding the variant.
#'
#'  \item \code{seq.context.width} The width of \code{seq.context} to the left
#'     of the variant. Does not include the "context base". TODO(Steve): do we need
#'     to modify this function so that it can handle indel callers that do
#'     not provide the "context base"?
#' }
#' @keywords internal
AddAndCheckSequenceID <- function(df, genome, flag.mismatches = FALSE) {

  genome <- NormalizeGenomeArg(genome)

  stopifnot(nchar(df$REF) != nchar(df$ALT)) # This has to be an indel, maybe a complex indel
  if (any(df$REF == "" | df$ALT == "")) {
    # Not sure how to handle this yet; the code may work with minimal adjustment
    stop("Cannot handle VCF with indel representation with one allele the empty string")
  } else {
    # We expect either eg ref = ACG, alt = A (deletion of CG) or
    # ref = A, alt = ACC (insertion of CC)
    stopifnot(substr(df$REF, 1, 1) == substr(df$ALT, 1, 1))
    complex.indels.to.remove <- which((nchar(df$REF) > 1 & (nchar(df$ALT) > 1)))
    if (length(complex.indels.to.remove > 0)) {
      warning("Removing complex indels", complex.indels.to.remove, "\n")
      print(df[ complex.indels.to.remove, 1:5])
      df <- df[ -complex.indels.to.remove, ]
    }
  }
  # First, figure out how much sequence context is needed.
  var.width <- abs(nchar(df$ALT) - nchar(df$REF))

  is.del <- nchar(df$ALT) <= nchar(df$REF)
  var.width.in.genome <- ifelse(is.del, var.width, 0)

  df$seq.context.width <- var.width * 6
  # 6 because we need to find out if the insertion or deletion is embedded
  # in up to 5 additonal repeats of the inserted or deleted sequence.
  # Then 6 to avoid possible future issues.

  # Extract sequence context from the reference genome

  # Check if the format of sequence names in df and genome are the same.
  # Internally ICAMS uses human chromosomes labeled as "1", "2", ... "X"...
  # However, BSgenome.Hsapiens.UCSC.hg38 has chromosomes labeled
  # "chr1", "chr2", ....
  vcf.chr.names <- unique(df$CHROM)
  if (!all(vcf.chr.names %in% seqnames(genome))) {
    tmp.chr <- paste0("chr", vcf.chr.names)
    if (!all(tmp.chr %in% seqnames(genome))) {
      stop("Cannot match chromosome names:\n",
           sort(vcf.chr.names), "\nversus\n", sort(seqnames(genome)))
    }

    chr.names <- paste0("chr", df$CHROM)
  } else {
    chr.names <- df$CHROM
  }
  # Create a GRanges object with the needed width.
  Ranges <-
    as(data.frame(chrom = chr.names,
                  start = df$POS - df$seq.context.width, # 10,
                  end = df$POS + var.width.in.genome + df$seq.context.width # 10
    ),
    "GRanges")

  df$seq.context <- getSeq(genome, Ranges, as.character = TRUE)

  seq.to.check <-
    substr(df$seq.context, df$seq.context.width + 1,
           df$seq.context.width + var.width.in.genome + 1)

  mismatches <- which(seq.to.check != df$REF)

  if (length(mismatches) > 0) {
    tmp.table <-
      data.frame(
        df$CHROM, df$POS, df$REF, df$ALT, df$seq.context, seq.to.check)
    tmp.table <- tmp.table[mismatches, ]
    cat("\n\nMismatches between VCF and reference sequence:\n\n")
    rows.to.print <-
      ifelse(flag.mismatches == 0,
             min(10, nrow(tmp.table)),
             min(nrow(tmp.table), flag.mismatches))
    cat("Showing", rows.to.print, "rows out of", nrow(tmp.table), "\n\n")
    print(tmp.table[1:rows.to.print, ])
    if (flag.mismatches > 0) {
      cat("\n\nDiscarding rows with mismatches\n\n")
      df <- df[-mismatches, ]
    } else {
      stop("Mismatch to reference genome sequence")
    }
  }

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

#' @title Return the length of microhomology at a deletion.
#'
#' @details
#'
#' This function is primarily for internal use, but we export it so that the
#' logic behind it will be documented for users.
#'
#' Example:
#'
#' \code{GGCTAGTT} aligned to \code{GGCTAGAACTAGTT} with
#' a deletion represented as:
#' \preformatted{
#'
#' GGCTAGAACTAGTT
#' GG------CTAGTT GGCTAGTT GG[CTAGAA]CTAGTT
#'                            ----   ----
#' }
#'
#' Presumed repair mechanism leading to this:
#'
#' \preformatted{
#'   ....
#' GGCTAGAACTAGTT
#' CCGATCTTGATCAA
#'
#' =>
#'
#'   ....
#' GGCTAG      TT
#' CC      GATCAA
#'         ....
#'
#' =>
#'
#' GGCTAGTT
#' CCGATCAA
#'
#' }
#'
#' The deletion caller can represent the
#' same deletion in several
#' different, but completely equivalent, ways.
#'
#' \preformatted{
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
#' }
#'
#' A deletion in a \emph{repeat} can also be represented
#' in several different ways. A deletion in a repeat
#' is abstractly equivalent to microhomology that
#' spans the entire deleted sequence. For example;
#'
#' \preformatted{
#' GACTAGCTAGTT
#' GACTA----GTT GACTAGTT GACTA[GCTA]GTT
#'                         *** -*** -
#' }
#'
#' is really a repeat
#'
#' \preformatted{
#' GACTAG----TT GACTAGTT GACTAG[CTAG]TT
#'                         **** ----
#'
#' GACT----AGTT GACTAGTT GACT[AGCT]AGTT
#'                         ** --** --
#' }
#'
#' \strong{But the function only flags this with a -1 return; it does not figure
#' out the repeat extent.}
#'
#' In the implementation, the function finds:
#'
#' \enumerate{
#'
#' \item The maximum match of undeleted sequence on left that is
#' identical to the right end of the deleted sequence, and
#'
#' \item The maximum match of undeleted sequence on the right this
#' is identical to the left end of the deleted sequence.
#'}
#'
#' The microhomology sequence is the concatenation of items
#' (1) and (2).
#'
#' @param context The deleted sequence plus ample surrounding
#'   sequence on each side (at least as long as \code{del.sequence}).
#'
#' @param deleted.seq The deleted sequence in \code{context}.
#'
#' @param pos The position of \code{del.sequence} in \code{context}.
#'
#' @param trace If > 0, cat various messages.
#'
#' @return The length of the maximum microhomology of \code{del.sequence}
#'   in \code{context}.
#'
#' @export
#'
FindDelMH <- function(context, deleted.seq, pos, trace = 0) {
  n <- nchar(deleted.seq)

  if (substr(context, pos, pos + n - 1) != deleted.seq) {
    stop("substr(context, pos, pos + n - 1) != deleted.seq\n",
         substr(context, pos, pos + n -1), " ", deleted.seq, "\n")
  }
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
    if (i == 1) {
      stop("There is a repeated ", deleted.seq,
           " to the left of the deleted ",
           deleted.seq)
    }
  }
  left.len <- n - i
  if (trace > 0 ) {
    cat("Left break", i, "\nleft.len =", left.len, "\n")
  }

  # Look for microhomology to the right in context.
  right.context <- substr(context, pos + n, (pos + 2 * n) - 1)
  right <- unlist(strsplit(x = right.context, ""))
  for (i2 in 1:n) {
    if (ds[i2] != right[i2]) break
    if (i2 == n) {
      stop("There is a repeated ", deleted.seq,
           " to the right of the deleted ",
           deleted.seq)
    }
  }
  right.len <- i2 - 1
  if (trace > 0) {
    cat("Right break", i2, "\nright.len =", right.len, "\n")
    cat(paste0(left.context, "[",
               deleted.seq, "]",
               right.context, "\n"))
    # Print out strings of ** and -- to indicated the sequences involved in the
    # microhomology; left.context and right.context are the same length as
    # deleted.seq
    cat(paste(c(
      rep(" ", n - left.len),
      rep("*", left.len),
      " ",
      rep("-", right.len),
      rep(" ", n - (left.len + right.len)),
      rep("*", left.len),
      " ",
      rep("-", right.len),
      "\n"
    ), collapse = ""))

  }
  if (left.len + right.len == n) {
    warning("There is unhandled cryptic repeat, returning -1")
    return(-1)
  }
  return (left.len + right.len)
}

#' @title Return the number of repeat units in which an insertion
#' is embedded.
#'
#' @param context A string into which \code{rep.unit.seq} was
#'  inserted at position \code{pos}.
#'
#' @param rep.unit.seq The inserted sequence and candidate repeat unit
#' sequence.
#'
#' @param pos \code{rep.unit.seq} is understood to be inserted between
#'   positions \code{pos} and \code{pos + 1}.
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

#' @title Given a deletion and its sequence context, categorize it.
#'
#' @param context The deleted sequence plus ample surrounding
#'   sequence on each side (at least as long as \code{del.seq}).
#'
#' @param del.seq The deleted sequence in \code{context}.
#'
#' @param pos The position of \code{del.sequence} in \code{context}.
#'
#' @param trace If > 0 cat information how the computation is carried out.
#
#' @return A string that is the canonical representation of the given deletion type
#'
#' @keywords internal
Canonicalize1DEL <- function(context, del.seq, pos, trace = 0) {
  # Is the deletion involved in a repeat?
  rep.count <- FindMaxRepeatDel(context, del.seq, pos)

  rep.count.string <- ifelse(rep.count >= 5, "5+", as.character(rep.count))
  deletion.size <- nchar(del.seq)
  deletion.size.string <- ifelse(deletion.size >= 5, "5+", as.character(deletion.size))

  # Category is "1bp deletion"
  if (deletion.size == 1) {
    if (del.seq == "G") del.seq <- "C"
    if (del.seq == "A") del.seq <- "T"
    return(paste0("DEL:", del.seq, ":1:", rep.count.string))
  }

  # Category is ">2bp deletion"
  if (rep.count > 0) {
    return(paste0("DEL:repeats:", deletion.size.string, ":", rep.count.string))
  }

  # We have to look for microhomology
  microhomology.len <- FindDelMH(context, del.seq, pos, trace = trace)
  if (microhomology.len == -1) {
    stop("Non-normalized deleted repeat: ",
         "necessary to use a different indel caller or enhance this code")
  }

  if (microhomology.len == 0) {
    stopifnot(rep.count.string == 0)
    # Categorize and return non-repeat, non-microhomology deletion
    return(paste0("DEL:repeats:", deletion.size.string, ":0"))
  }

  microhomology.len.str <-
    ifelse(microhomology.len >= 5, "5+", as.character(microhomology.len))

  return(paste0("DEL:MH:", deletion.size.string, ":", microhomology.len.str))
}

#' @title Given an insertion and its sequence context, categorize it.
#'
#' @param context The deleted sequence plus ample surrounding
#'   sequence on each side (at least as long as \code{ins.sequence} * 6).
#'
#' @param ins.sequence The deleted sequence in \code{context}.
#'
#' @param pos The position of \code{ins.sequence} in \code{context}.
#'
#' @param trace If > 0, then cat information how the computation is carried out.
#
#' @return A string that is the canonical representation of the given insertion type
#'
#' @keywords internal

Canonicalize1INS <- function(context, ins.sequence, pos, trace = 0) {
  if (trace > 0) {
   cat("Canonicalize1ID(", context, ",", ins.sequence, ",", pos, "\n")
  }
  rep.count <- FindMaxRepeatIns(context, ins.sequence, pos)
  rep.count.string <- ifelse(rep.count >= 5, "5+", as.character(rep.count))
  insertion.size <- nchar(ins.sequence)
  insertion.size.string <-
    ifelse(insertion.size >= 5, "5+", as.character(insertion.size))

  if (insertion.size == 1) {
    if (ins.sequence == "G") ins.sequence <- "C"
    if (ins.sequence == "A") ins.sequence <- "T"
    retval <-
      paste0("INS:", ins.sequence, ":1:", rep.count.string)
    if (trace > 0) cat(retval, "\n")
    return(retval)
  }
  retval <-
    paste0("INS:repeats:", insertion.size.string, ":", rep.count.string)
  if (trace > 0) cat(retval, "\n")
  return(retval)
}

#' @title Given a single insertion or deletion in context categorize it.
#'
#' @param context Ample surrounding
#'   sequence on each side of the insertion or deletion.
#'
#' @param ref The reference allele (vector of length 1)
#'
#' @param alt The alternative allele (vector of length 1)
#'
#' @param pos The position of \code{ins.or.del.seq} in \code{context}.
#'
#' @param trace If > 0 cat information how the computation is carried out.
#
#' @return A string that is the canonical representation of the type of the given
#'  insertion or deletion.
#'
#' @keywords internal
Canonicalize1ID <- function(context, ref, alt, pos, trace = 0) {
  if (trace > 0) cat("Canonicalize1ID(", context, ",", ref, ",", alt, ",", pos, "\n")
  if (nchar(alt) < nchar(ref)) {
    # A deletion
    return(Canonicalize1DEL(context, ref, pos + 1, trace))
  } else if (nchar(alt) > nchar(ref)) {
    # An insertion
    return(Canonicalize1INS(context, alt, pos, trace))
  } else {
    cat("Non-insertion / non-deletion found:", ref, alt, context, "\n")
    stop()
  }
}

#' @title Given vectors of insertions and deletions in contexts categorize them.
#'
#' @param context A vector of ample surrounding
#'   sequence on each side the variants
#'
#' @param ref Vector of reference alleles
#'
#' @param alt Vector of alternative alleles
#'
#' @param pos Vector of the positions of the insertions and deletions in \code{context}.
#'
#' @param trace If > 0 cat information how the computation is carried out.
#
#' @return A vector of strings that are the canonical representations
#'  of the given insertions and deletions.
#'
#' @importFrom utils head
#'
#' @keywords internal
CanonicalizeID <- function(context, ref, alt, pos, trace = 0) {

  if (trace > 1) {
    print(head(data.frame(
      context, left.pad = substr(context, 1, pos), pos, ref, alt,
      stringsAsFactors = FALSE)))
  }

  if (all(substr(ref, 1, 1) == substr(alt, 1, 1))) {
    ref <- substr(ref, 2, nchar(ref))
    alt <- substr(alt, 2, nchar(alt))
  } else {
    stopifnot(ref != "" | alt != "")
  }
  # trace = 1
  ret <- mapply(Canonicalize1ID, context, ref, alt, pos, trace)
  return(ret)
}

#' @title Create one column of an indel catalog from one VCF
#'
#' @param ID.vcf An in-memory VCF as a data.frame annotated by the
#'   \code{\link{AddAndCheckSequenceID}} function. It must only
#'   contain indels and must \strong{not} contain SNSs
#'   (single nucleotide/base substitutions), DNSs, or triplet
#'   base substitutions, etc.
#'
#'   One design decision for variant callers is the representation of "complex
#'   indels", e.g. mutations e.g. CAT > GC. Some callers represent this as C>G,
#'   A>C, and T>_. Others might represent it as CAT > CG. Multiple issues can
#'   arise. In PCAWG, overlapping indel/SNS calls from different callers were
#'   included in the indel VCFs.
#'
#' @param SNS.vcf An in-memory VCF as a data frame. Because we have to work with
#'   some PCAWG data, we will look for neighboring indels and indels adjoining
#'   SNS. That means this functions takes an SNS VCF and an ID VCF from the same
#'   sample.
#'
#' @param trace If > 0, various called functions cat information
#'   useful for debugging and testing. The larger the number, the
#'   more output.
#'
#' @return A list with two elements:
#'   ID.cat:   A 1-column matrix containing the mutation catalog information.
#'   problems: Locations of neighboring indels or indels neighboring SNS.
#'             In the future we might handle these depending on what we
#'             find in the indel calls from different variant callers. This
#'             is not implemented at present.
#'
#' @keywords internal
CreateOneColIDCatalog <- function(ID.vcf, SNS.vcf, trace = 0) {
  # TODO(Steve): more checking of the ID VCF here

  canon.ID <- CanonicalizeID(ID.vcf$seq.context,
                             ID.vcf$REF,
                             ID.vcf$ALT,
                             ID.vcf$seq.context.width + 1,
                             trace = trace)

  # Create the ID catalog matrix
  tab.ID <- table(canon.ID)

  row.order <- data.table(rn = ICAMS::catalog.row.order$ID)
  # TODO(Steve): can reduce use of data.table?

  ID.dt <- as.data.table(tab.ID)
  # ID.dt has two columns, names cannon.dt (from the table() function
  # and N (the count)

  ID.dt2 <-
    merge(row.order, ID.dt, by.x="rn", by.y="canon.ID", all = TRUE)
  ID.dt2[ is.na(N) , N := 0]
  stopifnot(setequal(unlist(ID.dt2$rn), ICAMS::catalog.row.order$ID))

  ID.mat <- as.matrix(ID.dt2[ , 2])
  rownames(ID.mat) <- ID.dt2$rn
  return(ID.mat[ICAMS::catalog.row.order$ID, , drop = F])
}


#' Create ID (insertion and deletion) catalog from ID VCFs
#'
#' @param list.of.vcfs List of in-memory VCFs. The list names will be
#' the sample ids in the output catalog.
#'
#' @param genome A genome argument as described in \code{\link{ICAMS}}.
#'
#' @return An ID (indel) catalog
#'
#' @export
VCFsToIDCatalogs <- function(list.of.vcfs, genome) {

  ncol <- length(list.of.vcfs)

  # Create a 0-column matrix with the correct row labels.
  catID <- matrix(0, nrow = length(ICAMS::catalog.row.order$ID), ncol = 0)
  rownames(catID) <- ICAMS::catalog.row.order$ID

  for (i in 1 : ncol) {
    ID <- list.of.vcfs[[i]]
    ID <- AddAndCheckSequenceID(ID, genome = genome)
    # Unlike the case for SNS and DNS, we do not
    # add transcript information.
    one.ID.column <- CreateOneColIDCatalog(ID)
    rm(ID)
    catID <- cbind(catID, one.ID.column)
  }

  colnames(catID) <- names(list.of.vcfs)
  return(catID)
}

