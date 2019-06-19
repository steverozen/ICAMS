#' "Collapse" a catalog.
#' 
#' @description 
#' \enumerate{
#' \item Take a mutational spectrum or signature catalog
#' that is based on a fined-grained set of features (for example, single-nucleotide
#' substitutions in the context of the preceding and following 2 bases).
#'
#' \item Collapse it to a catalog based on a coarser-grained set of features
#' (for example, single-nucleotide substitutions in the context of the
#' immediately preceding and following bases).
#' }
#'
#' \code{Collapse192CatalogTo96} Collapse an SBS 192 catalog
#' to an SBS 96 catalog.
#'
#' \code{Collapse1536CatalogTo96} Collapse an SBS 1536 catalog
#'  to an SBS 96 catalog.
#'
#' \code{Collapse144CatalogTo78} Collapse a DBS 144 catalog
#' to a DBS 78 catalog.
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}}.
#'
#' @return A catalog as defined in \code{\link{ICAMS}}.
#'
#' @name CollapseCatalog
NULL

#' @rdname CollapseCatalog
#' @export
Collapse192CatalogTo96 <- function(catalog) {
  dt192 <- data.table(catalog)
  dt192$rn <- PyrTri(rownames(catalog))
  dt96 <- dt192[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat96 <- as.matrix(dt96[, -1])
  rownames(mat96) <- dt96$rn
  mat96 <- mat96[ICAMS::catalog.row.order$SBS96, , drop = FALSE]

  cat96 <-
    as.catalog(
      object = mat96,
      ref.genome = attributes(catalog)$ref.genome,
      region = attributes(catalog)$region,
      catalog.type = attributes(catalog)$catalog.type,
      abundance = Collapse192AbundanceTo96(attributes(catalog)$abundance)
    )
  return(cat96)
}

#' @keywords internal
Collapse192AbundanceTo96 <- function(abundance192) {
  PyrTri <- function(string) {
    stopifnot(nchar(string) == rep(3, length(string)))
    output <-
      ifelse(substr(string, 2, 2) %in% c("A", "G"), revc(string), string)
    return(output)
  }
  dt <- data.table(abundance192)
  rownames(dt) <- names(abundance192)
  dt$rn <- PyrTri(rownames(dt))
  dt1 <- dt[, lapply(.SD, sum), by = rn, .SDcols = ]
  abundance96 <- unlist(dt1[, -1])
  names(abundance96) <- dt1$rn
  return(abundance96)
}

#' @rdname CollapseCatalog
#' @export
Collapse1536CatalogTo96 <- function(catalog) {
  dt <- data.table(catalog)
  rn <- rownames(catalog)

  # The next gsub replaces the string representing a
  # single-base mutation in pentanucleotide with the corresponding
  # sring for that mutation in a trinucleotide context.
  dt$rn <- gsub(".(...).(.)", "\\1\\2", rn, perl = TRUE)
  dt96 <- dt[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat96 <- as.matrix(dt96[, -1])
  rownames(mat96) <- dt96$rn
  mat96 <- mat96[ICAMS::catalog.row.order$SBS96, , drop = FALSE]
  
  cat96 <-
    as.catalog(
      object = mat96,
      ref.genome = attributes(catalog)$ref.genome,
      region = attributes(catalog)$region,
      catalog.type = attributes(catalog)$catalog.type,
      abundance = Collapse1536AbundanceTo96(attributes(catalog)$abundance)
    )
  return(cat96)
}

#' @keywords internal
Collapse1536AbundanceTo96 <- function(abundance1536) {
  dt <- data.table(abundance1536)
  rownames(dt) <- names(abundance1536)
  dt$rn <- substr(rownames(dt), 2, 4)
  dt1 <- dt[, lapply(.SD, sum), by = rn, .SDcols = ]
  abundance96 <- unlist(dt1[, -1])
  names(abundance96) <- dt1$rn
  return(abundance96)
}

#' @rdname CollapseCatalog
#' @export
Collapse144CatalogTo78 <- function(catalog) {
  dt144 <- data.table(catalog)
  ref <- substr(rownames(catalog), 1, 2)
  alt <- substr(rownames(catalog), 3, 4)
  dt144$rn <- CanonicalizeDBS(ref, alt)
  dt78 <- dt144[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat78 <- as.matrix(dt78[ , -1])
  rownames(mat78) <- dt78$rn
  mat78 <- mat78[ICAMS::catalog.row.order$DBS78, , drop = FALSE]

  cat78 <-
    as.catalog(
      object = mat78,
      ref.genome = attributes(catalog)$ref.genome,
      region = attributes(catalog)$region,
      catalog.type = attributes(catalog)$catalog.type,
      abundance = Collapse1536AbundanceTo96(attributes(catalog)$abundance)
    )
  return(cat78)
}

#' @keywords internal
IsDensity <- function(x) {
  return(x %in% c("density", "density.signature"))
}

#' @keywords internal
IsCounts <- function(x) {
  return(x %in% c("counts", "counts.signature"))
}

#' @keywords internal
IsSignature <- function(x) {
  return(x %in% c("counts.signature", "density.signature"))
}

#' @keywords internal
Collapse144AbundanceTo78 <- function(abundance144) {
  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  dt <- data.table(abundance144)
  rownames(dt) <- names(abundance144)
  dt$rn <- ifelse(rownames(dt) %in% canonical.ref, rownames(dt), revc(rownames(dt)))
  dt1 <- dt[, lapply(.SD, sum), by = rn, .SDcols = ]
  abundance78 <- unlist(dt1[, -1])
  names(abundance78) <- dt1$rn
  return(abundance78)
}

#' @keywords internal
CheckAndNormalizeTranCatArgs <-
  function(catalog, 
           target.ref.genome, 
           target.region, 
           target.catalog.type, 
           target.abundance) {
    
    StopIfNrowIllegal(catalog)
    
    s <- list(
      ref.genome    = attr(catalog, "ref.genome",   exact = TRUE),
      catalog.type = attr(catalog,  "catalog.type", exact = TRUE),
      abundance    = attr(catalog,  "abundance",    exact = TRUE),
      region       = attr(catalog,  "region",       exact = TRUE))
    
    if (is.null(target.ref.genome))   target.ref.genome   <- s$ref.genome
    if (is.null(target.catalog.type)) target.catalog.type <- s$catalog.type 
    StopIfCatalogTypeIllegal(target.catalog.type)
    if (is.null(target.region))       target.region       <-s$region
    StopIfRegionIllegal(target.region)

    inferred.abundance <- 
        InferAbundance(catalog,
                       target.ref.genome,
                       target.region,
                       target.catalog.type)
                   
    if (!is.null(target.abundance) &&
        !is.null(inferred.abundance)) {
      if (!all.equal(target.abundance,
                     inferred.abundance)) {
        stop("Caller supplied abundance is different from inferred abundance")
      }
    }
    if (is.null(target.abundance)) target.abundance <- inferred.abundance

    if (s$catalog.type == target.catalog.type) {
      if (all(s[["abundance"]] == target.abundance)) { 
        warning("Input and target catalog.type and abundance are the same\n",
                "no transformation")
      }
    }
    return(list(s = s, 
                t = 
                  list(ref.genome   = target.ref.genome,
                       catalog.type = target.catalog.type,
                       region       = target.region,
                       abundance    = target.abundance)))
  }


#' @keywords internal
IsTransformationLegal <- function(s, t) {
 
  if (IsSignature(s$catalog.type) &&
      !IsSignature(t$catalog.type)) {
    stop("cannot transform from a signture to a counts or density catalog")
  }
  
  if (IsDensity(s$catalog.type))
    return(TCFromDenSigDen(s, t))
  else
    return(TCFromCouSigCou(s, t))
}


#' Source catalog type is counts or counts.signature
#' 
#' counts.signature -> density.signature, counts.signature
#' counts -> anything
#' 
#' @keywords internal
TCFromCouSigCou <- function(s, t) {
  
  stopifnot(IsCounts(s$catalog.type))
  
  if (IsSignature(s$catalog.type)) {
    # counts.signature -> density.signatyre, counts.signature
    return(TCFromCouSig(s, t))
  } else {
    # counts -> anything
    return(TCFromCou(s, t))
  }
}

#' @keywords internal
#' 
#' counts.signature -> counts.signature, density.signature
TCFromCouSig <- function(s, t) {
  
  if (t$catalog.type == "counts.signature") {
    # counts.signature -> counts.signature
    if (is.null(s[["abundance"]])) {
      stop("Cannot transform from counts.signature -> counts.signature ",
           "source abundance is NULL")
    }
    if (is.null(t[["abundance"]])) {
      stop("Cannot transform from counts.signature -> counts.signature ",
           "target abundance is NULL")
    }
    if (all(s[["abundance"]] == t[["abundance"]])) {
      warning("Transformation from counts.signature -> counts.sinature",
              "with equal abundance is a null operation")
      return(NULL)
    }
    return(TRUE)
  } else if (t$catalog.type == "counts") {
    # counts.signature -> counts
    stop("Cannot transform from counts.signature -> counts")
  } else if (t$catalog.type == "density.signature") {
    # counts.signature -> density.signature
     if (is.null(s[["abundance"]])) {
       stop("Cannot transform from counts.signature -> denisity.signature ",
            "if target abundance is NULL")
     } else {
       return(TRUE)
     }

  } else if (t$catalog.type == "density") {
    # counts.signature -> density
    stop("Cannnot tranform from counts.signature -> density")
  }
  stop("programming error, unexpected catalog.type ", t$catalog.type)
}


#' @keywords internal
#' 
#' counts -> <anything>
TCFromCou <- function(s, t) {
  if (t$catalog.type == "counts.signature") {
    # counts -> counts.signature
    if (is.null(s[["abundance"]])) {
      if (!is.null(t[["abundance"]]))
      stop("Cannot transform from counts -> counts.signature ",
           "source abundance is NULL and target abundance is not NULL")
    }
    if (is.null(t[["abundance"]])) {
      if (!is.null(s[["abundance"]]))
      stop("Cannot transform from counts.signature -> counts.signature ",
           "target abundance is NULL")
    }
    return(TRUE)
  } else if (t$catalog.type == "counts") {
    # counts -> counts
    if (is.null(s[["abundance"]]) || is.null(t[["abundance"]])) {
      stop("Cannot transform from counts -> counts if either ",
           "source or target abundance is null")
    }
    if (all(s[["abundance"]] == t[["abundance"]])) {
      warning("Trasformation from counts -> counts with equal abundances ",
              "is a null operation")
      return(NULL)
    }
    return(TRUE)
        
  } else if (t$catalog.type == "density.signature") {
    # counts -> density.signature
    if (is.null(s[["abundance"]])) {
      stop("Cannot transform from counts -> denisity.signature ",
           "if source abundance is NULL")
    } 
    return(TRUE)
  } else if (t$catalog.type == "density") {
    # counts -> density
    if (is.null(s[["abundance"]])) {
      stop("Cannot transform from counts -> denisity ",
           "if source abundance is NULL")
    } 
    return(TRUE)
    
  }
  stop("programming error, unexpected catalog.type ", t$catalog.type)
}

#' density -> <anything>
#' density.signature -> density.signature, counts.signature
#' 
#' @keywords internal
TCFromDenSigDen <- function(s, t) {

  if (IsSignature(s$catalog.type))
    return(TCFromDenSig(s, t))
  else
    return(TCFromDen(s,t))
}

#' @keywords internal
#' 
#' density.signature -> anything
TCFromDenSig <- function (s, t) {
  if (t$catalog.type == "density") {
    # density.signature -> density
    stop("programming error: density.signature -> density")
  } else if (t$catalog.type == "density.signature") {
    # density.signature -> density.signature
    warning("Transform density.signature -> density signature",
            " is a null operation")
    return(NULL)
  } else if (t$catalog.type == "counts") {
    # density.signature -> counts
    stop("programming error: density.signature -> counts")
  } else if (t$catalog.type == "counts.signature") {
    # density.signature -> counts.signature
    if (is.null(t[["abundance"]])) {
      stop("Cannot transform density.signature -> counts.signature",
           " if target abundance is NULL")
    } else {  return(TRUE)   }
  }
  stop("programming error, unexpected catalog.type ", t$catalog.type)
}

#' @keywords internal
#' density -> anything
TCFromDen <- function (s, t) {  
  if (t$catalog.type == "density") {
    # density -> density
    warning("Transformation from density to density ",
            "is a null operation")
    return(NULL)
  } else if (t$catalog.type == "density.signature") { 
    # density -> density.signature
    return(TRUE)
  } else if (t$catalog.type == "counts") {
    # density -> counts
    if (is.null(t[["abundance"]])) {
      stop("Cannot tranform density -> counts; target abundance is NULL")
    } else {
      return(TRUE)
    }
  } else if (t$catalog.type == "counts.signature") {
    # density -> counts.signature
    if (is.null(t[["abundance"]])) {
      stop("Cannot tranform density -> counts.signature; target abundance is NULL")
    } else {
      return(TRUE)
    }    
  }
  stop("programming error, unexpected catalog.type ", t$catalog.type)
}


#' Transform between counts and density spectrum catalogs
#' and counts and density signature catalogs.
#'
#' @details Only the following transformations are legal:
#'
#' \enumerate{
#'
#' \item \code{counts -> counts} (used to transform
#'    between \code{target.ref.genome} and/or \code{target.region})
#'
#' \item \code{counts -> density}
#'
#' \item \code{counts -> (counts.signature, density.signature)}
#'
#' \item \code{density -> counts} (the semantics are to
#' infer the genome-wide or exome-wide counts based on the
#' densities)
#'
#' \item \code{density -> (counts.signature, density.signature)}
#'
#' \item \preformatted{counts.signature -> (counts.signature, density.signature)}
#'
#' \item \preformatted{density.signature -> counts.signature}
#'
#' \item \code{density.signature -> density.signature} (a null operation)
#'
#' \item \code{density -> density} (a null operation)
#' }
#'
#'
#' @param catalog An SBS or DBS catalog as described in \code{\link{ICAMS}};
#'  must \strong{not} be an ID (indel) catalog.
#'
#' @param target.ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param target.region A \code{region} argument; see \code{\link{as.catalog}}
#' and \code{\link{ICAMS}}.
#'
#' @param target.catalog.type A character string acting as a catalog type
#'   identifier, one of "counts", "density", "counts.signature",
#'   "density.signature"; see \code{\link{as.catalog}}.
#'   
#' @param target.abundance Optional, only needed when
#'   \code{target.ref.genome} is not one of the reference genomes
#'   known to ICAMS (see \code{\link{ICAMS}}. 
#'   The abundance should contain the counts of different
#'   source sequences for mutations. See
#'   \code{\link{all.abundance}}.
#'   
#' @return A catalog as defined in \code{\link{ICAMS}}.
#'
#' @export
TransformCatalog <-
  function(catalog, 
           target.ref.genome,
           target.region, 
           target.catalog.type,
           target.abundance = NULL) {
  
    # Check and normalize the arguments  
    args <-
      CheckAndNormalizeTranCatArgs(
        catalog             = catalog,
        target.ref.genome   = target.ref.genome,
        target.catalog.type = target.catalog.type,
        target.region       = target.region,
        target.abundance    = target.abundance)
    
    s <- args$s
    t <- args$t
    test <- IsTransformationLegal(s, t)
    # cat(test, "\n")

    # t$ref.genome <- NormalizeGenomeArg(t$ref.genome)
    
   if (t$catalog.type != "counts.signature") {

    } else {
      if (s$catalog.type == "counts.signature" && 
          attr(catalog, "region", exact = TRUE) == target.region) {
        return(catalog)
      } else if (s$catalog.type == "counts" && 
                 attr(catalog, "region", exact = TRUE) == target.region) {
        out <- apply(catalog, MARGIN = 2, function (x) x / sum(x))
        return(as.catalog(out, t$ref.genome,
                          target.region, target.catalog.type)) 
        
      }
    }
    
    if (s$catalog.type == "density.signature" &&
        target.catalog.type == "density.signature") {
      return(catalog)
    }

    if (attributes(catalog)$catalog.type == "density" &&
        target.catalog.type == "density") {
      return(catalog)
    }

    source.abundance <- attributes(catalog)$abundance
    
    if (is.null(target.abundance)) {
      inferred.abundance <-
        InferAbundance(catalog, 
                       t$ref.genome, 
                       target.region,
                       target.catalog.type)
      if (is.null(inferred.abundance)) stop("Cannot infer abundance")
      target.abundance <- inferred.abundance
    } else {
      target.abundance <- target.abundance
    }
    
    stopifnot(names(source.abundance) == names(target.abundance))

    factor <- target.abundance / source.abundance
    if (any(is.infinite(factor)) || any(is.na(factor))) {
      factor[is.infinite(factor)] <- 0
      factor[is.na(factor)] <- 0
    }

    names(factor) <- names(target.abundance)
    out.catalog <- catalog

    for(source.n.mer in names(source.abundance)) {
      # CAUTION: this loop depends on how mutations are encoded in
      # the row names of catalogs!
      
      # For 96 and 192 SBS, source.n.mer is e.g. "ACT" (for the encoding of ACT >
      # AGT as "ACTG"); for SBS1536 the n-mer for AACAG > AATAG is AACAG, in the
      # encoding AACAGT. For DBS78 and DBS144 TGGA represents TG >GA, and the
      # source n-mer is TG. For DBS136, TTGA represents TTGA > TNNA, and the
      # source n-mer is TTGA.
      
      # First, get the rows with the given source.n.mer
      rows <- grep(paste("^", source.n.mer, sep=''), rownames(out.catalog))
      # Then update those rows using the factor for that source.n.mer
      
      out.catalog[rows, ] <- out.catalog[rows, ] * factor[source.n.mer]
    }

    if (target.catalog.type %in% c("counts.signature", "density.signature")) {
      out2 <- apply(out.catalog, MARGIN = 2, function (x) x / sum(x))
      return(as.catalog(out2, t$ref.genome,
                        target.region, target.catalog.type))
    } else {
      return(as.catalog(out.catalog, t$ref.genome,
                        target.region, target.catalog.type))
    }
  }

#' Standardize the chromosome name annotations for a data frame.
#'
#' @param df A data frame whose first column contains the Chromosome name
#'
#' @return A data frame whose Chromosome names are only in the form of 1:22, "X"
#'   and "Y".
#'
#' @keywords internal
StandardChromName <- function(df) {
  # Is there any row in df whose Chromosome names have "GL"?
  if (sum(grepl("GL", df[[1]])) > 0) {
    df <- df[-grep("GL", df[[1]]), ]
  }

  # Is there any row in df whose Chromosome names have "KI"?
  if (sum(grepl("KI", df[[1]])) > 0) {
    df <- df[-grep("KI", df[[1]]), ]
  }

  # Is there any row in df whose Chromosome names have "random"?
  if (sum(grepl("random", df[[1]])) > 0) {
    df <- df[-grep("random", df[[1]]), ]
  }

  # Is there any row in df whose Chromosome names are "Hs37D5"?
  if (sum(grepl("^Hs", df[[1]])) > 0) {
    df <- df[-grep("^Hs", df[[1]]), ]
  }

  # Is there any row in df whose Chromosome names contain "M"?
  if (sum(grepl("M", df[[1]])) > 0) {
    df <- df[-grep("M", df[[1]]), ]
  }

  # Remove the "chr" character in the Chromosome's name
  df[, 1] <- sub(pattern = "chr", replacement = "", df[[1]])

  return(df)
}

#' Create a transcript range file from the raw GFF3 File
#'
#' @param file The name/path of the raw GFF3 File, or a complete URL.
#'
#' @importFrom  stringi stri_split_fixed
#'
#' @return A data table which contains chromosome name, start, end position,
#'   strand information and gene name. It is keyed by chrom, start, and
#'   end. Only genes that are associated with a CCDS ID are kept for
#'   transcriptional strand bias analysis.
#'
#' @keywords internal
CreateTransRanges <- function(file) {
  df <- read.csv(file, header = FALSE, fill = TRUE, nrows = 20)
  # Count the number of comment lines
  n <- sum(grepl("#", df[, 1]))

  # Read in the raw GFF3 File while skipping the comment lines
  dt <- data.table::fread(file, header = FALSE, sep = "\t", fill = TRUE, skip = n)

  # Extract the gene ID associated with CCDS ID
  dt1 <- dt[grep("CCDS", dt$V9), ]
  info <- stri_split_fixed(dt1$V9, ";")
  gene.id.CCDS.idx <- lapply(info, grep, pattern = "gene_id")
  ExtractInfo <- function(idx, list1, list2) {
    return(list1[[idx]][list2[[idx]]])
  }
  gene.id.CCDS <-
    unique(sapply(1:length(info), ExtractInfo,
                  list1 = info, list2 = gene.id.CCDS.idx))
  gene.id.CCDS <- sapply(stri_split_fixed(gene.id.CCDS, "="), "[", 2)

  # Extract the gene ID and gene name from entries in dt which belong to the
  # type "gene"
  dt2 <- dt[dt$V3 == "gene", ]
  attributes.info <- stri_split_fixed(dt2$V9, ";")
  gene.id.idx <- lapply(attributes.info, grep, pattern = "gene_id")
  gene.name.idx <- lapply(attributes.info, grep, pattern = "gene_name")
  len <- length(attributes.info)
  gene.id <-
    sapply(1:len, ExtractInfo, list1 = attributes.info, list2 = gene.id.idx)
  gene.id <- sapply(stri_split_fixed(gene.id, "="), "[", 2)
  gene.name <-
    sapply(1:len, ExtractInfo, list1 = attributes.info, list2 = gene.name.idx)
  gene.name <- sapply(stri_split_fixed(gene.name, "="), "[", 2)
  dt3 <- dt2[, c("gene.id", "gene.name") := .(gene.id, gene.name)]

  # Only keep the entries in dt3 with gene ID that is associated with CCDS ID.
  # Select the necessary columns and standardize the chromosome names.
  dt4 <- dt3[gene.id %in% gene.id.CCDS, c(1, 4, 5, 7, 11)]
  colnames(dt4) <- c("chrom", "start", "end", "strand", "gene.name")
  dt5 <- StandardChromName(dt4)

  # Reorder dt5 according to chromosome name, start and end position
  chrOrder <- c((1:22), "X", "Y")
  dt5$chrom <- factor(dt5$chrom, chrOrder, ordered = TRUE)
  setkeyv(dt5, c("chrom", "start", "end"))

  # Combine the overlapping ranges of the same gene if there is any
  dt6 <- dt5[, count := .N, by = .(chrom, gene.name)]
  dt7 <- dt6[count == 1, ]
  dt8 <- dt6[count != 1, ]
  if (nrow(dt8) == 0) {
    return(dt7[, c(1:5)])
  } else {
    gr <- GenomicRanges::makeGRangesFromDataFrame(dt8, keep.extra.columns = TRUE)
    gr1 <- GenomicRanges::reduce(gr, with.revmap = TRUE)
    GetGeneName <- function(idx, names){
      return(names[idx])
    }
    mat <- t(sapply(gr1$revmap, FUN = GetGeneName, names = gr$gene.name))
    unique.gene.name <- apply(mat, MARGIN = 1, FUN = unique)
    GenomicRanges::mcols(gr1) <- unique.gene.name
    dt9 <- as.data.table(gr1)[, c(1:3, 5:6)]
    dt10 <- rbind(dt7[, c(1:5)], dt9, use.names = FALSE)
    return(setkeyv(dt10, c("chrom", "start", "end")))
  }
}

#' Create exome transcriptionally stranded regions
#'
#' @param file Path to a SureSelect BED file which contains unstranded exome
#'   ranges.
#'
#' @param trans.ranges A data.table which contains transcript range and strand
#'   information. Please refer to \code{\link{TranscriptRanges}} for more
#'   details.
#'
#' @import data.table
#'
#' @return A data table which contains chromosome name, start, end position,
#'   strand information. It is keyed by chrom, start, and
#'   end.
#'
#' @keywords internal
CreateExomeStrandedRanges <- function(file, trans.ranges) {
  exome.ranges <- ReadBedRanges(file)
  colnames(exome.ranges) <- c("chrom", "exome.start", "exome.end")

  # Remove ranges which fall on transcripts on both strands and get
  # transcriptionally stranded regions
  trans.ranges <- RemoveRangesOnBothStrand(trans.ranges)

  # Find range overlaps between the exome.ranges and trans.ranges
  dt <- foverlaps(exome.ranges, trans.ranges, type = "any", mult = "all")
  dt1 <- dt[!is.na(strand)]

  # Find out and remove exome ranges that fall on transcripts on both strands
  dt2 <- dt1[, bothstrand := "+" %in% strand && "-" %in% strand,
             by = .(chrom, exome.start, exome.end)]
  dt3 <- dt2[bothstrand == FALSE]

  # One unstranded exome range can be represented by more than 1 row in dt3 if it
  # overlaps with multiple transcriptionally stranded regions. When creating the
  # exome transcriptionally stranded regions, we only need to count those ranges once.
  dt4 <- dt3[, count := .N, by = .(chrom, exome.start, exome.end)]
  dt4 <- dt3[, .(strand = strand[1]), by = .(chrom, exome.start, exome.end)]

  # Intersect dt4 ranges with the transcriptionally stranded regions
  gr1 <- GenomicRanges::makeGRangesFromDataFrame(dt4)
  gr2 <- GenomicRanges::makeGRangesFromDataFrame(trans.ranges)
  gr3 <- GenomicRanges::intersect(gr1, gr2)

  dt5 <- as.data.table(gr3)[, c(1:3, 5)]
  colnames(dt5)[1] <- "chrom"
  chrOrder <- c((1:22), "X", "Y")
  dt5$chrom <- factor(dt5$chrom, chrOrder, ordered = TRUE)
  return(setkeyv(dt5, c("chrom", "start", "end")))
}

#' @keywords internal
PyrTri <- function(mutstring) {
  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  output <-
    ifelse(substr(mutstring, 2, 2) %in% c("A", "G"),
           paste0(revc(substr(mutstring, 1, 3)),
                  revc(substr(mutstring, 4, 4))),
           mutstring)
  return(output)
}

#' @keywords internal
PyrPenta <- function(mutstring) {
  stopifnot(nchar(mutstring) == rep(6, length(mutstring)))
  output <-
    ifelse(substr(mutstring, 3, 3) %in% c("A", "G"),
           paste0(revc(substr(mutstring, 1, 5)),
                  revc(substr(mutstring, 6, 6))),
           mutstring)
  return(output)
}

#' Reverse complement every string in \code{string.vec}.
#' 
#' Based on \code{\link{reverseComplement}}.
#' Handles IUPAC ambiguity codes but not "u" (uracil).
#' (see <https://en.wikipedia.org/wiki/Nucleic_acid_notation>).
#'
#' @param string.vec A character vector.
#'
#' @importFrom Biostrings reverseComplement DNAStringSet
#'
#' @return A character vector with the reverse complement of every
#'   string in \code{string.vec}.
#'
#' @export
#' 
#' @examples 
#' # Suggested by Steve:
#' 
#' revc("aTgc") # GCAT
#' 
#' # A vector and strings with ambiguity codes
#' revc(c("ATGC", "aTGc", "wnTCb")) # GCAT GCAT VGANW
#' 
#' \dontrun{
#' revc("ACGU") # An error}
#' 
#' # Reverse complement a character of DNA sequence
#' string1 <- "ATGC"
#' string2 <- "atgc"
#' revc(string1)        # "GCAT"
#' revc(string2)        # "GCAT"
#' 
#' # Reverse complement a vector of character of DNA sequences
#' strings <- c("ATGC", "TCGT")
#' revc(strings)        # "GCAT" "ACGA"
revc <- function(string.vec) {
  return(
    as.character(reverseComplement(DNAStringSet(string.vec)))
  )
}

#' @title Reverse complement strings that represent stranded SBSs
#'
#' @param mutstring A vector of 4-character strings representing
#' stranded SBSs in trinucleotide context,
#' for example "AATC" represents AAT > ACT mutations.
#'
#' @return Return the vector of
#' reverse complements of the first 3 characters
#' concatenated with the reverse complement of the
#' last character, e.g. "AATC" returns "ATTG".
#'
#' @keywords internal
RevcSBS96 <- function(mutstring) {
  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  context <- revc(substr(mutstring, 1, 3))
  target  <- revc(substr(mutstring, 4, 4))
  return(paste0(context, target))
}

#' @title Reverse complement strings that represent stranded DBSs
#'
#' @param mutstring A vector of 4-character strings representing
#' stranded DBSs, for example "AATC" represents AA > TC mutations.
#'
#' @return Return the vector of
#' reverse complements of the first 2 characters
#' concatenated with the reverse complement of the second
#' 2 characters, e.g. "AATC" returns "TTGA".
#'
#' @keywords internal
RevcDBS144 <- function(mutstring) {
  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  context <- revc(substr(mutstring, 1, 2))
  target  <- revc(substr(mutstring, 3, 4))
  return(paste0(context, target))
}

#' Read transcript ranges and strand information from a gff3 format file.
#' Use this one for the new, cut down gff3 file (2018 11 24)
#'
#' @param file Path to the file with the transcript information with 1-based
#'   start end positions of genomic ranges.
#'
#' @return A data.table keyed by chrom, start, and end.
#'
#' @keywords internal
ReadTranscriptRanges <- function(file) {
  dt <- data.table::fread(file)
  colnames(dt) <- c("chrom", "start", "end", "strand", "gene.name")
  chrOrder <- c((1:22), "X", "Y")
  dt$chrom <- factor(dt$chrom, chrOrder, ordered = TRUE)
  data.table::setkeyv(dt, c("chrom", "start", "end"))
  return(dt)
}

#' Read chromosome and position information from a bed format file.
#'
#' @param file Path to the file in bed format.
#'
#' @return A data.table keyed by chrom, start, and end.
#'
#' @keywords internal
ReadBedRanges <- function(file) {
  dt <- data.table::fread(file)
  dt1 <- StandardChromName(dt[, 1:3])
  colnames(dt1) <- c("chrom", "start", "end")

  # Delete duplicate entries in the BED file
  dt2 <- dplyr::distinct(dt1, chrom, start, end, .keep_all = TRUE)

  # Bed file are 0 based start and 1 based end (an oversimplification).
  # We need to add 1L and not 1, otherwise the column turns to a double
  # we get a warning from data.table.
  dt2$start <- dt2$start + 1L

  chrOrder <- c((1:22), "X", "Y")
  dt2$chrom <- factor(dt2$chrom, chrOrder, ordered = TRUE)
  return(data.table::setkeyv(dt2, c("chrom", "start", "end")))
}

#' Create trinucleotide abundance
#'
#' @param file Path to the file with the nucleotide abundance information with 3
#'   base pairs.
#'
#' @return A numeric vector whose names indicate 32 different types of 3 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreateTrinucAbundance <- function(file) {
  dt <- fread(file)
  colnames(dt) <- c("3bp", "occurrences")
  dt$type <-
    ifelse(substr(dt[[1]], 2, 2) %in% c("A", "G"), revc(dt[[1]]), dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create stranded trinucleotide abundance
#'
#' @param file Path to the file with the nucleotide abundance information with 3
#'   base pairs.
#'
#' @return A numeric vector whose names indicate 64 different types of 3 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreateStrandedTrinucAbundance <- function(file) {
  dt <- fread(file)
  colnames(dt) <- c("3bp", "occurrences")
  dt$type <- dt[[1]]
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create dinucleotide abundance
#'
#' @param file Path to the file with the nucleotide abundance information with 2
#'   base pairs.
#'
#' @import data.table
#'
#' @return A numeric vector whose names indicate 10 different types of 2 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreateDinucAbundance <- function(file) {
  dt <- fread(file)
  colnames(dt) <- c("2bp", "occurrences")
  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  dt$type <-
    ifelse((dt[[1]]) %in% canonical.ref, dt[[1]], revc(dt[[1]]))
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create stranded dinucleotide abundance
#'
#' @param file Path to the file with the nucleotide abundance information with 2
#'   base pairs.
#'
#' @import data.table
#'
#' @return A numeric vector whose names indicate 16 different types of 2 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreateStrandedDinucAbundance <- function(file) {
  dt <- fread(file)
  colnames(dt) <- c("2bp", "occurrences")
  dt$type <- dt[[1]]
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create tetranucleotide abundance
#'
#' @param file Path to the file with the nucleotide abundance information with 4
#'   base pairs.
#'
#' @import data.table
#'
#' @return A numeric vector whose names indicate 136 different types of 4 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreateTetranucAbundance <- function(file) {
  dt <- fread(file)
  colnames(dt) <- c("4bp", "occurrences")
  dt$type <- CanonicalizeQUAD(dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create pentanucleotide abundance
#'
#' @param file Path to the file with the nucleotide abundance information
#'   with 5 base pairs.
#'
#' @import data.table
#'
#' @return A numeric vector whose names indicate 512 different types of 5 base
#'   pairs combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreatePentanucAbundance <- function(file) {
  dt <- fread(file)
  colnames(dt) <- c("5bp", "occurrences")
  dt$type <-
    ifelse(substr(dt[[1]], 3, 3) %in% c("A", "G"), revc(dt[[1]]), dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Take strings representing a genome and return the \code{\link{BSgenome}} object.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @return If \code{ref.genome} is 
#' a \code{\link{BSgenome}} object, return it.
#' Otherwise return the \code{\link{BSgenome}} object identified by the
#' string \code{ref.genome}.
#'
#' @keywords internal
NormalizeGenomeArg <- function(ref.genome) {
  
  if (is.null(ref.genome)) stop("Need a non-NULL ref.genome")
  if (class(ref.genome) == "BSgenome") return(ref.genome)
  
  stopifnot(class(ref.genome) == "character")
  
  if (ref.genome %in%
      c("GRCh38", "hg38", "BSgenome.Hsapiens.UCSC.hg38")) {
    ref.genome <- BSgenome.Hsapiens.UCSC.hg38
  } else if (ref.genome %in%
             c("GRCh37", "hg19", "BSgenome.Hsapiens.1000genomes.hs37d5")) {
    ref.genome <- BSgenome.Hsapiens.1000genomes.hs37d5
  } else {
    stop("Unrecoginzed ref.genome:\n", ref.genome,
         "\nNeed NULL or a BSgenome reference genome\n",
         "or one of the character strings GRCh38, hg38, GRCh37, hg19")
  }
  
  return(ref.genome)
}

#' Stop if \code{region} is illegal for an in-transcript catalogs
#'
#' @param region The region to test (a character string)
#' 
#' @keywords internal

StopIfTranscribedRegionIllegal <- function(region) {
  if (!region %in% c("transcript", "exome")) {
    stop("Require region to be one of transcript, exome, or unknown for ",
         "SBS192Catalog and DBS144Catalog\n",
         "Got ", region)
  }
}

#' Stop if the number of rows in \code{object} is illegal
#' 
#' @param object A \code{catalog}, numeric \code{matrix}, or numeric \code{data.fram}
#' 
#' @keywords internal
StopIfNrowIllegal <- function(object) {
  if(!nrow(object) %in% c(96, 192, 1536, 78, 144, 136, 83)) {
    stop("\nThe number of rows in the input object must be one of\n",
         "96 192 1536 78 144 136 83\ngot ", nrow(object))
  }

}

#' Stop if \code{region} is illegal.
#'
#' @param region Character string to check.
#' 
#' @keywords internal
StopIfRegionIllegal <- function(region) {
  if (!region %in% c("genome", "exome", "transcript", "unknown")) {
    stop("Unrecoginzed region identifier: ", region,
         "\nNeed one of genome, exome, transcript, unknown")
  }
 return(NULL)
}

#' Stop if \code{catalog.type} is illegal.
#'
#' @param catalog.type Character string to check.
#'
#' @keywords internal
StopIfCatalogTypeIllegal <- function(catalog.type) {
  if (is.null(catalog.type)) return (NULL)
  if (!catalog.type %in% c("counts", "density",
                           "counts.signature", "density.signature")) {
    stop("Unrecoginzed catalog type identifier: ", catalog.type,
         "\nNeed one of counts, density, counts.signature, density.signature")
  }
  return(NULL)
}

#' Infer the class of catalog in a file.
#'
#' @param file Path to a catalog on disk in the standardized format.
#'
#' @return A character string with class appropriate for the catalog on disk.
#'
#' @keywords internal
InferClassOfCatalogForRead <- function(file) {
  cos <- data.table::fread(file)
  StopIfNrowIllegal(cos)
  class.string <- InferCatalogClassString(cos) 
  return(structure("ClassofCatalog", class = class.string))
}

#' @keywords internal
InferCatalogClassString <- function(object) {
  
  StopIfNrowIllegal(object)
  nrow <- nrow(object)

  if(nrow == 96)   return("SBS96Catalog")
  if(nrow == 192)  return("SBS192Catalog")
  if(nrow == 1536) return("SBS1536Catalog")
  if(nrow == 78)   return("DBS78Catalog")
  if(nrow == 144)  return("DBS144Catalog")
  if(nrow == 136)  return("DBS136Catalog")
  if(nrow == 83)   return("IndelCatalog")

}

#' Test if object is \code{BSgenome.Hsapiens.1000genome.hs37d5}.
#'
#' @param x Object to test.
#' 
#' @return TRUE if \code{x} is \code{BSgenome.Hsapiens.1000genome.hs37d5}.
#' 
#' @keywords internal
IsGRCh37 <- function(x) {
  if (is.null(x)) return(FALSE)
  return(NormalizeGenomeArg(x)@pkgname == 
           "BSgenome.Hsapiens.1000genomes.hs37d5")
}

#' Test if object is \code{BSgenome.Hsapiens.1000genome.hs37d5}.
#'
#' @param x Object to test.
#' 
#' @return TRUE if \code{x} is \code{BSgenome.Hsapiens.1000genome.hs37d5}.
#' 
#' @keywords internal
IsGRCh38 <- function(x) {
  if (is.null(x)) return(FALSE)
  return(NormalizeGenomeArg(x)@pkgname == 
           "BSgenome.Hsapiens.UCSC.hg38")
}


#' Infer \code{abundance} given a matrix-like \code{object} and additional information.
#'
#' @param object A numeric matrix, numeric data frame, or \code{catalog}.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @param catalog.type A character string for \code{catalog.type}
#' as described in \code{\link{ICAMS}}.
#'
#' @return A value that can be set as the abundance attribute of
#' a \code{catalog} (which may be \code{NULL} if no abundance
#' can be inferred).
#'
#' @keywords internal
InferAbundance <- function(object, ref.genome, region, catalog.type) {
    
    StopIfNrowIllegal(object)
    StopIfRegionIllegal(region)
    StopIfCatalogTypeIllegal(catalog.type)

    if (IsDensity(catalog.type)) {
      ab <- flat.abundance[[as.character(nrow(object))]]
      stopifnot(!is.null(ab))
      return(ab)
    }
    
    if (is.null(ref.genome)) return(NULL)
    ref.genome <- NormalizeGenomeArg(ref.genome)
    ab <- all.abundance[[ref.genome@pkgname]]
    if (is.null(ab)) return(NULL)

    ab2 <- ab[[region]]
    if (is.null(ab2)) return(NULL)

    ab3 <- ab2[[as.character(nrow(object))]]

    return(ab3)

  }
  

#' Create a catalog from a numeric \code{matrix} or numeric \code{data.frame}.
#'
#' @param object A numeric \code{matrix} or numeric \code{data.frame}. 
#'  This object must have
#'   rownames to denote the mutation types. See \code{\link{CatalogRowOrder}}
#'   for more details.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region A character string designating a region, one of
#' \code{genome}, \code{transcript}, \code{exome}, \code{unknown};
#' see \code{\link{ICAMS}}.
#'
#' @param catalog.type One of "counts", "density", "counts.signature",
#'   "density.signature".
#'
#' @param abundance Optional, only needed when \code{ref.genome}
#'  is not one of
#'  the reference genomes known to ICAMS (see \code{\link{ICAMS}}.
#'  The abundance should
#'  contain the counts of different source sequences for mutations. \cr
#'  See \code{\link{all.abundance}}.
#'
#' @return A catalog as described in \code{\link{ICAMS}}.
#'
#' @export
as.catalog <- function(object, 
                       ref.genome = NULL, 
                       region = "unknown", 
                       catalog.type = "counts", 
                       abundance = NULL) {
  if (!is.matrix(object)) {
    if (is.data.frame(object)) {
      object <- as.matrix(object)
    } else {
      stop("object must be numeric matrix or data frame")
    }
  }
  stopifnot(mode(object) == "numeric")

  stopifnot(!is.null(rownames(object)))
  
  StopIfRegionIllegal(region)
  
  StopIfNrowIllegal(object)
  
  StopIfCatalogTypeIllegal(catalog.type)
  
  if (!is.null(ref.genome)) {
    ref.genome <- NormalizeGenomeArg(ref.genome)
  }
  attr(object, "ref.genome") <- ref.genome

  attr(object, "catalog.type") <- catalog.type
  
  if (is.null(abundance)) {
    abundance <- InferAbundance(object, ref.genome, region, catalog.type)
  } 
  attr(object, "abundance") <- abundance

  class.string  <- InferCatalogClassString(object)
  class(object) <- append(class(object), class.string, after = 0)
  class(object) <- unique(attributes(object)$class)
  
  # TODO(Steve): check that his has been moved 
  # to the VCF -> catalog functions
  if (attributes(object)$class[1] %in% c("SBS192Catalog", "DBS144Catalog")) {
    StopIfTranscribedRegionIllegal(region)
  }
  attr(object, "region") <- region
  
  return(object)
}

#' Generate all possible k-mers of length k.
#'
#' @param k Length of k-mers (k>=2)
#'
#' @return Character vector containing all possible k-mers.
#'
#' @keywords internal
GenerateKmer <- function(k) {
  base <- c("A", "C", "G", "T")
  list.of.base <- list()
  for (i in 1:k){
    list.of.base[[i]] <- base
  }
  permutation <- expand.grid(list.of.base, stringsAsFactors = FALSE)
  all.kmer.list <- character(nrow(permutation))
  for (i in 1:k) {
    all.kmer.list <- stringi::stri_c(all.kmer.list, permutation[[i]])
  }
  all.kmer.list <- stringi::stri_sort(all.kmer.list)
}

#' Generate an empty matrix of k-mer abundance
#'
#' @param k Length of k-mers (k>=2)
#'
#' @return An empty matrix of k-mer abundance
#'
#' @keywords internal
GenerateEmptyKmerCounts <- function(k) {
  all.kmer.list <- GenerateKmer(k)
  kmer.counts <- matrix(0, nrow = length(all.kmer.list))
  rownames(kmer.counts) <- all.kmer.list
  colnames(kmer.counts) <- "Freq"
  return(kmer.counts)
}

#' Generate k-mer abundance from given nucleotide sequences
#'
#' @param sequences A vector of nucleotide sequences
#'
#' @param k Length of k-mers (k>=2)
#'
#' @return Matrix of the counts of each k-mer inside \code{sequences}
#'
#' @keywords  internal
GetSequenceKmerCounts <- function(sequences, k) {
  kmer.counts <- GenerateEmptyKmerCounts(k)

  for(start_idx in 1:k){
    temp.seqs <- substring(sequences, start_idx, nchar(sequences))
    temp.kmers <-
      stringi::stri_extract_all_regex(
        temp.seqs, pattern = paste(rep(".", each = k), collapse = ""))
    temp.kmers <- unlist(temp.kmers)
    temp.kmer.counts <- data.frame(table(temp.kmers))
    row.names(temp.kmer.counts) <- temp.kmer.counts[, 1]
    if (any(grepl("N", temp.kmer.counts[, 1]))) {
      temp.kmer.counts <- temp.kmer.counts[-grep("N", temp.kmer.counts[, 1]), ]
    }

    kmer.counts[row.names(temp.kmer.counts), ] <-
      kmer.counts[row.names(temp.kmer.counts), ]  +
      temp.kmer.counts$Freq
  }
  return(kmer.counts)
}

#' Generate k-mer abundance from a given genome
#'
#' @param k Length of k-mers (k>=2)
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param filter.path If given, homopolymers will be masked from
#'   genome(sequence). Only simple repeat masking is accepted now.
#'
#' @importFrom GenomicRanges GRanges
#'
#' @importFrom IRanges IRanges
#'
#' @return Matrix of the counts of each k-mer across the \code{ref.genome}
#'
#' @keywords internal
GetGenomeKmerCounts <- function(k, ref.genome, filter.path) {
  kmer.counts <- GenerateEmptyKmerCounts(k)
  genome <- NormalizeGenomeArg(ref.genome)

  # Remove decoyed chromosomes and mitochondrial DNA
  chr.list <- seqnames(genome)[which(nchar(seqnames(genome)) <= 5)]
  if (any(grepl("M", chr.list))) {
    chr.list <- chr.list[-grep("M", chr.list)]
  }

  if (!missing(filter.path)) {
    filter.df <- fread(filter.path, header = F, stringsAsFactors = F)
    filter.df <- filter.df[filter.df$V6 <= 6]
    filter.df <- StandardChromName(filter.df[, 2:ncol(filter.df)])
    # Check whether chromosome names in filter.df are the same as in ref.genome
    if (!(seqnames(genome)[1] %in% filter.df$V2)){
      filter.df$V2 <- paste0("chr", filter.df$V2)
    }
  }

  print("Start counting by chromosomes")

  for (idx in 1:length(chr.list)) {
    print(chr.list[idx])

    if (!missing(filter.path)) {
      chr.filter.df <- filter.df[which(filter.df$V2 == chr.list[idx]), ]
      filter.chr <- with(chr.filter.df, GRanges(V2, IRanges(V3 + 1, V4)))
      genome.ranges.chr <-
        GRanges(chr.list[idx],
                IRanges(1, as.numeric(GenomeInfoDb::seqlengths(genome)[idx])))
      filtered.genome.ranges.chr <- GenomicRanges::setdiff(genome.ranges.chr, filter.chr)
      genome.seq <- BSgenome::getSeq(genome, filtered.genome.ranges.chr,
                                     as.character = TRUE)
      #Filter shorter homopolymer and microsatellites by regex
      genome.seq <- gsub(homopolymer.ms.regex.pattern, "N", genome.seq)

    } else {
      genome.seq <- BSgenome::getSeq(genome, chr.list[idx], as.character = TRUE)
    }

    kmer.counts <- kmer.counts + GetSequenceKmerCounts(genome.seq, k)
  }
  return(kmer.counts)
}

#' Remove ranges that fall on both strands
#'
#' @param stranded.ranges A keyed data table which has stranded ranges information.
#' It has four columns: chrom, start, end and strand.
#'
#' @return A data table which has removed ranges that fall on both strands from
#'   the input \code{stranded.ranges}.
#'
#' @keywords internal
RemoveRangesOnBothStrand <- function(stranded.ranges) {
  dt <- as.data.table(stranded.ranges)
  dt.plus <- dt[strand == "+", ]
  dt.minus <- dt[strand == "-", ]
  gr.plus <-
    GenomicRanges::makeGRangesFromDataFrame(dt.plus, keep.extra.columns = TRUE)
  gr.minus <-
    GenomicRanges::makeGRangesFromDataFrame(dt.minus, keep.extra.columns = TRUE)
  gr.plus.reduced <- GenomicRanges::reduce(gr.plus, with.revmap = TRUE)
  gr.minus.reduced <- GenomicRanges::reduce(gr.minus, with.revmap = TRUE)

  # Find ranges that have transcripts on both strands and remove these
  # ranges from each strand
  gr <-
    GenomicRanges::intersect(gr.plus.reduced, gr.minus.reduced,
                             ignore.strand = TRUE)

  if (length(gr) != 0) {
    gr1 <- gr
    gr1@strand@values <- "+"
    gr2 <- gr
    gr2@strand@values <- "-"
    gr3 <- GenomicRanges::setdiff(gr.plus.reduced, gr1)
    gr4 <- GenomicRanges::setdiff(gr.minus.reduced, gr2)
  } else {
    gr3 <- gr.plus.reduced
    gr4 <- gr.minus.reduced
  }

  # Get a new data table which does not have ranges on both strands
  dt1 <- as.data.table(gr3)
  dt2 <- as.data.table(gr4)
  dt3 <- rbind(dt1, dt2)
  dt4 <- dt3[, c(1:3, 5)]
  colnames(dt4) <- c("chrom", "start", "end", "strand")
  return(setkeyv(dt4, c("chrom", "start", "end")))
}

#' Generate stranded k-mer abundance from a given genome and gene annotation file
#'
#' @param k Length of k-mers (k>=2)
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param filter.path If given, homopolymers will be masked from
#'   genome(sequence). Only simple repeat masking is accepted now.
#'
#' @param stranded.ranges A keyed data table which has stranded ranges
#'   information. It has four columns: chrom, start, end and strand.
#'
#' @importFrom stats start end
#'
#' @importFrom GenomicRanges GRanges
#'
#' @importFrom IRanges IRanges
#'
#' @return Matrix of the counts of each stranded k-mer across the \code{ref.genome}
#'
#' @keywords internal
GetStrandedKmerCounts <- function(k, ref.genome, stranded.ranges, filter.path) {
  stranded.ranges <- RemoveRangesOnBothStrand(stranded.ranges)
  stranded.ranges <- StandardChromName(stranded.ranges)
  genome <- NormalizeGenomeArg(ref.genome)
  kmer.counts <- GenerateEmptyKmerCounts(k)

  # Check whether chromosome names in stranded.ranges are the same as in ref.genome
  if (!(seqnames(genome)[1] %in% stranded.ranges$chrom)) {
    stranded.ranges$chrom <- paste0("chr", stranded.ranges$chrom)
  }

  if (!missing(filter.path)) {
    filter.df <- fread(filter.path, header = F, stringsAsFactors = F)
    filter.df <- filter.df[filter.df$V6 <= 6]
    filter.df <- StandardChromName(filter.df[, 2:ncol(filter.df)])
    # Check whether chromosome names in filter.df are the same as in ref.genome
    if (!(seqnames(genome)[1] %in% filter.df$V2)){
      filter.df$V2 <- paste0("chr", filter.df$V2)
    }
  }

  print("Start counting by chromosomes")

  for (chr in unique(stranded.ranges$chrom)) {
    print(chr)
    temp.stranded.ranges <- stranded.ranges[stranded.ranges$chrom == chr, ]
    stranded.ranges.chr <-
      with(temp.stranded.ranges,
           GRanges(chrom, IRanges(start, end), strand = strand))

    # Remove the overlapping ranges in stranded.ranges.chr
    stranded.ranges.chr <- IRanges::reduce(stranded.ranges.chr)

    if (!missing(filter.path)) {
      chr.filter.df <- filter.df[which(filter.df$V2 == chr), ]

      # Add strand information to chr.filter.df because
      # the setdiff function is strand specific
      filter.chr <-
        c(with(chr.filter.df, GRanges(V2, IRanges(V3 + 1, V4), strand = "+")),
          with(chr.filter.df, GRanges(V2, IRanges(V3 + 1, V4), strand = "-")))
      filtered.stranded.ranges.chr <-
        GenomicRanges::setdiff(stranded.ranges.chr, filter.chr)

      stranded.seq <- BSgenome::getSeq(genome, filtered.stranded.ranges.chr,
                                       as.character = TRUE)
      # Filter shorter homopolymer and microsatellites by regex
      stranded.seq <- gsub(homopolymer.ms.regex.pattern, "N", stranded.seq)
    } else {
      stranded.seq <- BSgenome::getSeq(genome, stranded.ranges.chr,
                                       as.character = TRUE)
    }
    kmer.counts <- kmer.counts + GetSequenceKmerCounts(stranded.seq, k)
  }
  return(kmer.counts)
}

#' Generate exome k-mer abundance from a given reference genome
#'
#' @param k Length of k-mers (k>=2)
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param exome.range A keyed data table which has exome ranges information.
#' It has three columns: chrom, start and end.
#'
#' @param filter.path If given, homopolymers will be masked from
#'   genome(sequence). Only simple repeat masking is accepted now.
#'
#' @importFrom stats start end
#'
#' @importFrom GenomicRanges GRanges
#'
#' @importFrom IRanges IRanges
#'
#' @return Matrix of the counts of exome k-mer across the \code{ref.genome}
#'
#' @keywords internal
GetExomeKmerCounts <- function(k, ref.genome, exome.ranges, filter.path) {
  exome.ranges <- StandardChromName(exome.ranges)
  genome <- NormalizeGenomeArg(ref.genome)
  kmer.counts <- GenerateEmptyKmerCounts(k)

  # Check whether chromosome names in exome.ranges are the same as in ref.genome
  if (!(seqnames(genome)[1] %in% exome.ranges$chrom)) {
    exome.ranges$chrom <- paste0("chr", exome.ranges$chrom)
  }

  if (!missing(filter.path)) {
    filter.df <- fread(filter.path, header = F, stringsAsFactors = F)
    filter.df <- filter.df[filter.df$V6 <= 6]
    filter.df <- StandardChromName(filter.df[, 2:ncol(filter.df)])
    # Check whether chromosome names in filter.df are the same as in ref.genome
    if (!(seqnames(genome)[1] %in% filter.df$V2)){
      filter.df$V2 <- paste0("chr", filter.df$V2)
    }

  }
  print("Start counting by chromosomes")

  for (chr in unique(exome.ranges$chrom)) {
    print(chr)
    temp.exome.ranges <- exome.ranges[exome.ranges$chrom == chr, ]
    exome.range.chr <-
      with(temp.exome.ranges, GRanges(chrom, IRanges(start, end)))

    # Remove the overlapping ranges in exome.range.chr
    exome.range.chr <- IRanges::reduce(exome.range.chr)

    if (!missing(filter.path)) {
      chr.filter.df <- filter.df[which(filter.df$V2 == chr), ]
      filter.chr <- with(chr.filter.df, GRanges(V2, IRanges(V3 + 1, V4)))
      filtered.exome.range.chr <-
        GenomicRanges::setdiff(exome.range.chr, filter.chr)
      exome.seq <- BSgenome::getSeq(genome, filtered.exome.range.chr,
                                    as.character = TRUE)
      #Filter shorter homopolymer and microsatellites by regex
      exome.seq <- gsub(homopolymer.ms.regex.pattern, "N", exome.seq)

    } else {
      exome.seq <- BSgenome::getSeq(genome, exome.range.chr,
                                    as.character = TRUE)
    }
    kmer.counts <- kmer.counts + GetSequenceKmerCounts(exome.seq, k)
  }
  return(kmer.counts)
}

# Redefine the [ methods for catalogs
#' @export
`[.SBS96Catalog` <- function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==
                                1) {
  y <- NextMethod("[")
  if (class(y) %in% c("integer", "numeric")) {
    return(y)
  } else {
    class(y) <- class(x)
    for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
      attr(y, at) <- attr(x, at, exact = TRUE)
    }
    return(y)
  }
}

#' @export
`[.SBS192Catalog` <- function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==
                                1) {
  y <- NextMethod("[")
  if (class(y) %in% c("integer", "numeric")) {
    return(y)
  } else {
    class(y) <- class(x)
    for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
      attr(y, at) <- attr(x, at, exact = TRUE)
    }
    return(y)
  }
}

#' @export
`[.SBS1536Catalog` <- function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==
                                 1) {
  y <- NextMethod("[")
  if (class(y) %in% c("integer", "numeric")) {
    return(y)
  } else {
    class(y) <- class(x)
    for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
      attr(y, at) <- attr(x, at, exact = TRUE)
    }
    return(y)
  }
}

#' @export
`[.DBS78Catalog` <- function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==
                                  1) {
  y <- NextMethod("[")
  if (class(y) %in% c("integer", "numeric")) {
    return(y)
  } else {
    class(y) <- class(x)
    for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
      attr(y, at) <- attr(x, at, exact = TRUE)
    }
    return(y)
  }
}

#' @export
`[.DBS144Catalog` <- function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==
                                1) {
  y <- NextMethod("[")
  if (class(y) %in% c("integer", "numeric")) {
    return(y)
  } else {
    class(y) <- class(x)
    for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
      attr(y, at) <- attr(x, at, exact = TRUE)
    }
    return(y)
  }
}

#' @export
`[.DBS136Catalog` <- function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==
                                 1) {
  y <- NextMethod("[")
  if (class(y) %in% c("integer", "numeric")) {
    return(y)
  } else {
    class(y) <- class(x)
    for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
      attr(y, at) <- attr(x, at, exact = TRUE)
    }
    return(y)
  }
}

#' @export
`[.IndelCatalog` <- function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==
                                 1) {
  y <- NextMethod("[")
  if (class(y) %in% c("integer", "numeric")) {
    return(y)
  } else {
    class(y) <- class(x)
    for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
      attr(y, at) <- attr(x, at, exact = TRUE)
    }
    return(y)
  }
}

# Redefine the cbind methods for catalogs
#' @export
`cbind.SBS96Catalog` <- function (..., deparse.level = 1) {
  x <- base::cbind.data.frame(..., deparse.level = deparse.level)
  x <- data.matrix(x)
  class(x) <- class(..1)
  for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
    attr(x, at) <- attr(..1, at, exact = TRUE)
  }
  return(x)
}

#' @export
`cbind.SBS192Catalog` <- function (..., deparse.level = 1) {
  x <- base::cbind.data.frame(..., deparse.level = deparse.level)
  x <- data.matrix(x)
  class(x) <- class(..1)
  for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
    attr(x, at) <- attr(..1, at, exact = TRUE)
  }
  return(x)
}

#' @export
`cbind.SBS1536Catalog` <- function (..., deparse.level = 1) {
  x <- base::cbind.data.frame(..., deparse.level = deparse.level)
  x <- data.matrix(x)
  class(x) <- class(..1)
  for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
    attr(x, at) <- attr(..1, at, exact = TRUE)
  }
  return(x)
}

#' @export
`cbind.DBS78Catalog` <- function (..., deparse.level = 1) {
  x <- base::cbind.data.frame(..., deparse.level = deparse.level)
  x <- data.matrix(x)
  class(x) <- class(..1)
  for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
    attr(x, at) <- attr(..1, at, exact = TRUE)
  }
  return(x)
}

#' @export
`cbind.DBS144Catalog` <- function (..., deparse.level = 1) {
  x <- base::cbind.data.frame(..., deparse.level = deparse.level)
  x <- data.matrix(x)
  class(x) <- class(..1)
  for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
    attr(x, at) <- attr(..1, at, exact = TRUE)
  }
  return(x)
}

#' @export
`cbind.DBS136Catalog` <- function (..., deparse.level = 1) {
  x <- base::cbind.data.frame(..., deparse.level = deparse.level)
  x <- data.matrix(x)
  class(x) <- class(..1)
  for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
    attr(x, at) <- attr(..1, at, exact = TRUE)
  }
  return(x)
}

#' @export
`cbind.IndelCatalog` <- function (..., deparse.level = 1) {
  x <- base::cbind.data.frame(..., deparse.level = deparse.level)
  x <- data.matrix(x)
  class(x) <- class(..1)
  for (at in c("ref.genome", "catalog.type", "abundance", "region")) {
    attr(x, at) <- attr(..1, at, exact = TRUE)
  }
  return(x)
}
