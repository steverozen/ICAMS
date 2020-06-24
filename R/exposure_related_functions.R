#' @title Read an exposure matrix from a file
#'
#' @param file CSV file containing an exposure matrix.
#'
#' @param check.names Passed to \code{\link[utils]{read.csv}}.
#' \strong{IMPORTANT}: If \code{TRUE} this will replace the double
#' colon in identifiers of the form <tumor_type>::<sample_id>
#' with two periods (i.e. <tumor_type>..<sample_id>.
#' If \code{check.names} is true, generate a warning
#' if double colons were present.
#'
#' @return Matrix of exposures.
#'
#' @importFrom utils read.csv
#' 
#' @export
ReadExposure <- function(file, check.names = TRUE) {
  if (check.names) {
    headers <- read.csv(file, nrow = 1, header = FALSE, stringsAsFactors = FALSE)
    double.colon <- grep("::", unlist(headers)[-1], fixed = TRUE)
    if (length(double.colon) > 0) {
      warning(":: in sample ID replaced by ..; suggest calling with check.names = FALSE")
    }
  }
  retval <- read.csv(file, row.names = 1, check.names = check.names)
  if (any(duplicated(colnames(retval)))) {
    stop("There is duplicated column name in the input file")
  }
  return(data.matrix(retval))
}

#' @title Write an exposure matrix to a file
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @param file File to which to write the exposure matrix (as a CSV file).
#'
#' @importFrom utils write.csv
#' 
#' @export
WriteExposure <- function(exposure, file) {
  old.digits <- getOption("digits")
  options(digits = 22)
  write.csv(exposure, file, row.names = TRUE)
  on.exit(options(digits = old.digits)) 
}

#' Sort columns of an exposure matrix from largest to smallest (or vice versa)
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'   
#' @param decreasing If \code{TRUE}, sort from largest to smallest.
#' 
#' @return The original \code{exposure} with columns sorted.
#'
#' @export
SortExposure <- function(exposure, decreasing = TRUE) {
  retval <- exposure[, order(colSums(exposure), decreasing = decreasing)]
  return(retval)
}

#' Plot a single exposure plot
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'    signatures in rows and samples in columns. Rownames are taken
#'    as the signature names and column names are taken as the
#'    sample IDs. If you want \code{exp} sorted from largest to smallest
#'    use \code{\link{SortExp}}. Do not use column names that start
#'    with multiple underscores. The exposures will often be mutation
#'    counts, but could also be e.g. mutations per megabase.
#'
#' @param plot.proportion Plot exposure proportions rather than counts.
#'
#' @param plot.legend If \code{TRUE}, plot a legend.
#'
#' @param ... Parameters passed to \code{\link[graphics]{barplot}}.
#'
#' @import graphics 
#' 
#' @return An invisible numeric vector giving the coordinates of all the bar
#'   midpoints drawn, useful for adding to the graph.
#'
#' @keywords internal
PlotExposureInternal <-
  function(exposure, # This is actually the exposure "counts"
           plot.proportion = FALSE,
           plot.legend     = TRUE,
           ...
  ) {
    exposure <- as.matrix(exposure) # in case it is a data frame
    num.sigs  <- nrow(exposure)
    num.samples <- ncol(exposure)
    
    three.dots <- list(...)
    if (is.null(three.dots$col)) {
      if (num.sigs <= 8) {
        three.dots$col <- 
          c("red", "black", "grey", "yellow", "blue", "brown", "green4", "skyblue")
      } else {
        # Lots of signatures; use shaded lines to differentiate
        three.dots$col <- grDevices::rainbow(num.sigs, alpha = 1)
      }
    }
    if (num.sigs <= 12) {
      p.dense = -1 # will repeat as needed, -1 = solid
      p.angle = 0  # ditto
    } else {
      # Lots of signatures; use shaded lines to differentiate
      p.dense = c(-1, 50, 50, 50, 50) # will repeat as needed, -1 = solid
      p.angle = c(0, 45, 135, 45, 135)  # ditto
    }
    # For legend, we need to put in reverse order. make sure this matches
    # order used in barplot
    num.repeats = ceiling(num.sigs/length(p.dense))
    p.dense.rev = rev(rep(p.dense, num.repeats)[1:num.sigs])
    p.angle.rev = rev(rep(p.angle, num.repeats)[1:num.sigs])
    
    
    # l.cex = if (num.sigs > 15) 1.2 else 1 # char expansion for legend (was 0.7)
    # direction = 2 # 1 = always horizontal, 2 = perpendicular to axis
    
    if (plot.proportion) {
      # Matrix divided by vector goes col-wise, not row-wise, so transpose twice
      plot.what <- t(t(exposure)/colSums(exposure))
      
      # If some columns in exposure matrix are "padded", i.e. all the values
      # are 0, then the division above will make NaNs. We try to replace all
      # NaN with 0
      plot.what[is.na(plot.what)] <- 0
    } else {
      plot.what <- exposure
    }
    
    # plot.what <-
    #  cbind(plot.what,
    #        matrix(0,
    #               nrow = nrow(plot.what),
    #               ncol = round(ncol(plot.what) * 0.15)))
    
    # Ignore column names; we'll plot them separately to make them fit
    bp = do.call(
      barplot,
      args = c(list(height   = plot.what,
                    las      = 1,
                    #yaxt     = "s",
                    xaxt     = "n", # Do not plot the X axis
                    density  = p.dense,
                    angle    = p.angle,
                    border   = "white", # ifelse(num.samples>200,NA,1),
                    #cex.main = 1.2,
                    xlim     = c(0, round(ncol(plot.what) * 1.25)),
                    ylim     = c(0, max(colSums(plot.what)) * 1.1)),
               three.dots))
    
    # Get max y values for plot region, put legend at top right
    #dims = par("usr") # c(x.min, x.max, y.min, y.max)
    #y.max = dims[4]
    
    if (plot.legend) {
      # Less space between rows (y.intersp), and between box & label (x.intersp)
      # reverse the order, so sig 1 is at bottom (to match the stacked bar graph)
      #legend.x <- dims[2] * 0.95 #ncol(exposure) * 0.7   # Nanhai, we could pass in legend.x and legend.y as optional arguments from the caller
      #legend.y <- y.max
      legend(#x         = legend.x,
        #y         = legend.y,
        "right",
        inset     = -0.08,
        legend    = rev(row.names(exposure)),
        density   = p.dense.rev,
        angle     = p.angle.rev,
        #bg        = NA,
        xpd       = NA,
        fill      = three.dots$col[num.sigs:1],
        x.intersp = 0.3, #0.4,
        y.intersp = 1, #0.8,
        bty       = "n",
        border    = "white",
        cex       = 0.9,
        title     = "Signature")
      #text(x= legend.x, y = legend.y, labels = "Signature", 
      #xpd = NA, adj = -0.09)#-0.09)
    }
    
    # Now add sample names, rotated to hopefully fit
    # don't even try to show all if there are too many
    if (num.samples <= 200) {
      if (length(bp) < 50) {
        size.adj <- 0.75 
      } else if (length(bp) < 80) {
        size.adj <- 0.65  
      } else if (length(bp) < 100) {
        size.adj <- 0.4 # .5
      } else if (length(bp) < 120) {
        size.adj <- 0.4
      } else if (length(bp) < 150) {
        size.adj <- 0.3
      } else {
        size.adj <- 0.3
      }
      cnames <- colnames(exposure)
      cnames <- sub("_____.*", "", cnames)
      mtext(cnames, side = 1, at = bp, las = 2, cex = size.adj)
    }
    
    invisible(bp)
  }

#' Plot exposures in multiple plots each with a manageable number of samples
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs. If you want
#'   \code{exposure} sorted from largest to smallest use \code{\link{SortExp}}.
#'   Do not use column names that start with multiple underscores. The exposures
#'   will often be mutation counts, but could also be e.g. mutations per
#'   megabase.
#'
#' @param plot.proportion Plot exposure proportions rather than counts.
#'
#' @param samples.per.line Number of samples to show in each plot.
#'
#' @param ... Other arguments passed to \code{\link{PlotExposure}}. If
#'   \code{ylab} is not included, it defaults to a value depending on
#'   \code{plot.proportion}. If \code{col} is not supplied the function tries to
#'   do something reasonable.
#'
#' @export
PlotExposure <- function(exposure,
                         samples.per.line    = 30,
                         plot.proportion     = FALSE,
                         ...
) {
  new.xlim = c(0, samples.per.line * 1.25)
  args <- list(...)
  ylab <- args$ylab
  if (is.null(ylab)) {
    ylab <- ifelse(plot.proportion,
                   "Proportion of mutations",
                   "Number of mutations")
  }
  
  n.sample <- ncol(exposure)
  num.ranges <- n.sample %/% samples.per.line
  size.of.last.range <- n.sample %% samples.per.line
  if (size.of.last.range > 0) {
    padding.len <- samples.per.line - size.of.last.range
    padding <- matrix(0,nrow = nrow(exposure), ncol = padding.len)
    # The column names starting with lots of underscore
    # will not be plotted in the final output.
    colnames(padding) <- paste("_____", 1:ncol(padding), sep = "_")
    exposure <- cbind(exposure, padding)
    starts <- 0:num.ranges * samples.per.line + 1
  } else {
    starts <- 0:(num.ranges - 1) *samples.per.line + 1
  }
  ends   <- starts + samples.per.line - 1
  
  plot.legend <- TRUE
  for (i in 1:length(starts)) {
    PlotExposureInternal(exposure[ , starts[i]:ends[i]],
                         plot.proportion = plot.proportion,
                         plot.legend    = plot.legend,
                         ...)
    #plot.legend <- FALSE
  }
}