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
#' 
#' @examples 
#' file <- system.file("extdata",
#'                     "synthetic.exposure.csv",
#'                     package = "ICAMS")
#' exposure <- ReadExposure(file)
ReadExposure <- function(file, check.names = FALSE) {
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
#' 
#' @examples 
#' file <- system.file("extdata",
#'                     "synthetic.exposure.csv",
#'                     package = "ICAMS")
#' exposure <- ReadExposure(file)
#' WriteExposure(exposure, file = file.path(tempdir(), "synthetic.exposure.csv"))
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
#' 
#' @examples 
#' file <- system.file("extdata",
#'                     "synthetic.exposure.csv",
#'                     package = "ICAMS")
#' exposure <- ReadExposure(file)
#' exposure.sorted <- SortExposure(exposure)
SortExposure <- function(exposure, decreasing = TRUE) {
  retval <- exposure[, order(colSums(exposure), decreasing = decreasing),
                     drop = FALSE]
  return(retval)
}

#' Plot a single exposure plot
#'
#' @param exposure Exposures as a numerical \code{matrix} (or \code{data.frame})
#'   with signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs. If you want
#'   \code{exposure} sorted from largest to smallest, use
#'   \code{\link{SortExposure}}. Do not use column names that start with
#'   multiple underscores. The exposures will often be mutation counts, but
#'   could also be e.g. mutations per megabase.
#'
#' @param plot.proportion Plot exposure proportions rather than counts.
#'
#' @param xlim,ylim Limits for the x and y axis. If \code{NULL}(default), the
#'   function tries to do something reasonable.
#'
#' @param legend.x,legend.y The x and y co-ordinates to be used to position the
#'   legend.
#'   
#' @param cex.legend A numerical value giving the amount by which legend
#'   plotting text and symbols should be magnified relative to the default.
#'   
#' @param ... Other arguments passed to \code{\link[graphics]{barplot}}. If
#'   \code{ylab} is not included, it defaults to a value depending on
#'   \code{plot.proportion}. If \code{col} is not supplied the function tries to
#'   do something reasonable.
#'   
#' @import graphics 
#'
#' @return An \strong{invisible} list whose first element is a logic value
#'   indicating whether the plot is successful. The second element is a numeric
#'   vector giving the coordinates of all the bar midpoints drawn, useful for
#'   adding to the graph.
#'
#' @keywords internal
PlotExposureInternal <-
  function(exposure, # This is actually the exposure "counts"
           plot.proportion = FALSE,
           xlim            = NULL,
           ylim            = NULL,
           legend.x        = NULL,
           legend.y        = NULL,
           cex.legend      = 0.9,
           ...
  ) {
    exposure <- as.matrix(exposure) # In case it is a data frame
    num.sigs  <- nrow(exposure)
    num.samples <- ncol(exposure)
    args <- list(...)
    
    ylab <- args$ylab
    if (is.null(ylab)) {
      ylab <- ifelse(plot.proportion,
                     "Proportion of mutations",
                     "Number of mutations")
    }
    
    if (is.null(args$col)) {
      if (num.sigs <= 8) {
        args$col <- 
          c("red", "black", "grey", "yellow", "blue", "brown", "green4", "skyblue")
      } else {
        # Lots of signatures; use shading lines to differentiate
        args$col <- grDevices::rainbow(num.sigs, alpha = 1)
      }
    }
    
    if (num.sigs <= 12) {
      # Specify the density of shading lines, in lines per inch, for the bars or
      # bar components. Will repeat as needed, non-positive values of density
      # inhibit the drawing of shading lines.
      density <- -1 
      # Sepcify the slope of shading lines, given as an angle in degrees
      # (counter-clockwise), for the bars or bar components. Will repeat as
      # needed.
      angle <- 0  
    } else {
      # Lots of signatures; use shading lines to differentiate
      density <- c(-1, 50, 50, 50, 50) # Will repeat as needed
      angle <- c(0, 45, 135, 45, 135)  # Ditto
    }
    # For legend, we need to put density and angle in reversed order to make
    # sure this matches order used in barplot.
    num.repeats <- ceiling(num.sigs / length(density))
    density.rev <- rev(rep(density, num.repeats)[1:num.sigs])
    angle.rev <- rev(rep(angle, num.repeats)[1:num.sigs])
    
    if (plot.proportion) {
      # Matrix divided by vector goes column-wise, not row-wise, so transpose twice
      plot.what <- t(t(exposure)/colSums(exposure))
      
      # If some columns in exposure matrix are "padded", i.e. all the values
      # are 0, then the division above will make NaNs. We try to replace all
      # NaN with 0
      plot.what[is.na(plot.what)] <- 0
    } else {
      plot.what <- exposure
    }
    
    # Set the limits for x axis, leave extra space for plotting legend
    if (is.null(xlim)) {
      xmax <- round(ncol(plot.what) * 1.7)
      xlim <- c(0, xmax)
    } else {
      xmax <- xlim[2]
    }
    
    if (is.null(ylim)) {
      ymax <- round(max(colSums(plot.what))) 
      ylim <- c(0, ymax)
    } else {
      ymax <- ylim[2]
    }
    
    # Ignore column names, we'll plot them separately to make them fit.
    mp <- do.call(
      barplot,
      args = c(list(height   = plot.what,
                    ylab     = ylab,
                    # The amount of space left before each bar
                    space    = xmax * 0.01, 
                    xaxs     = "i", # No extra spacing at each end of x axis
                    xaxt     = "n", # Do not plot the X axis
                    yaxt     = "n", # Do not plot the Y axis
                    density  = density,
                    angle    = angle,
                    border   = "white", 
                    xlim     = xlim,
                    ylim     = ylim),
               args))
    
    # Get locations for y axis annotations
    y.axis.values <- seq(0, ymax, ymax/4)
    if (plot.proportion) {
      y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
    } else {
      y.axis.labels <- round(y.axis.values)
    }
    
    # Draw y axis and labels
    Axis(side = 2, at = y.axis.values, las = 1, labels = FALSE)
    text(x = -xmax* 0.03, y = y.axis.values, labels = y.axis.labels, 
         las = 1, adj = 1, xpd = NA)
    
    # Setting the default parameters for legend plotting
    if (is.null(legend.x)) {
      legend.x <- xmax * 0.9
    }
    if (is.null(legend.y)) {
      legend.y <- ymax 
    }
    
    legend(x         = legend.x,
           y         = legend.y,
           legend    = rev(row.names(exposure)),
           density   = density.rev,
           angle     = angle.rev,
           xpd       = NA,
           fill      = rev(args$col),
           x.intersp = 0.3, 
           y.intersp = 1, 
           bty       = "n",
           border    = "white",
           cex       = cex.legend,
           title     = "Signature")
    
    # Now add sample names, rotated to hopefully fit,
    # don't even try to show all if there are too many
    if (num.samples <= 200) {
      if (length(mp) < 50) {
        size.adj <- 0.75 
      } else if (length(mp) < 80) {
        size.adj <- 0.65  
      } else if (length(mp) < 100) {
        size.adj <- 0.4 
      } else if (length(mp) < 120) {
        size.adj <- 0.4
      } else if (length(mp) < 150) {
        size.adj <- 0.3
      } else {
        size.adj <- 0.3
      }
      cnames <- colnames(exposure)
      cnames <- sub("_____.*", "", cnames)
      mtext(cnames, side = 1, line = 0.38, at = mp, las = 2, cex = size.adj)
    }
    
    invisible(list(plot.success = TRUE, mp.coordinates = mp))
  }

#' Plot exposures in multiple plots each with a manageable number of samples
#'
#' @inheritParams PlotExposureInternal
#'
#' @param samples.per.line Number of samples to show in each plot.
#' 
#' @param xlim,ylim Limits for the x and y axis. If \code{NULL}(default), the
#'   function tries to do something reasonable.
#'
#' @param legend.x,legend.y The x and y co-ordinates to be used to position the
#'   legend.
#'   
#' @return An \strong{invisible} list whose first element is a logic value
#'   indicating whether the plot is successful. The second element is a numeric
#'   vector giving the coordinates of all the bar midpoints drawn, useful for
#'   adding to the graph.
#'
#' @export
#' 
#' @examples 
#' file <- system.file("extdata",
#'                     "synthetic.exposure.csv",
#'                     package = "ICAMS")
#' exposure <- ReadExposure(file)
#' PlotExposure(exposure[, 1:30])
PlotExposure <- function(exposure,
                         samples.per.line = 30,
                         plot.proportion  = FALSE,
                         xlim             = NULL,
                         ylim             = NULL,
                         legend.x         = NULL,
                         legend.y         = NULL,
                         cex.legend       = 0.9,
                         ...
) {
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
  
  for (i in 1:length(starts)) {
    list <- PlotExposureInternal(exposure[ , starts[i]:ends[i]],
                                 plot.proportion = plot.proportion,
                                 xlim            = xlim,
                                 ylim            = ylim,
                                 legend.x        = legend.x,
                                 legend.y        = legend.y,
                                 cex.legend      = cex.legend,
                                 ...             = ...)
  }
  invisible(list(plot.success = TRUE, mp.coordinates = list$mp.coordinates))
}

#' Plot exposures in multiple plots each with a manageable number of samples to PDF
#' 
#' @inheritParams PlotExposureInternal
#' 
#' @param file The name of the PDF file to be produced.
#' 
#' @param mfrow A vector of the form \code{c(nr, nc)}.
#' Subsequent figures will be drawn in an \code{nr}-by-\code{nc}
#' array on the device by rows.
#' 
#' @param mar A numerical vector of the form \code{c(bottom,
#' left, top, right)} which gives the number of lines of margin to be
#' specified on the four sides of the plot.
#' 
#' @param oma A vector of the form \code{c(bottom, left, top,
#' right)} giving the size of the outer margins in lines of text.
#' 
#' @param samples.per.line Number of samples to show in each plot.
#' 
#' @param xlim,ylim Limits for the x and y axis. If \code{NULL}(default), the
#'   function tries to do something reasonable.
#'
#' @param legend.x,legend.y The x and y co-ordinates to be used to position the
#'   legend.
#'   
#' @return An \strong{invisible} list whose first element is a logic value
#'   indicating whether the plot is successful. The second element is a numeric
#'   vector giving the coordinates of all the bar midpoints drawn, useful for
#'   adding to the graph.
#' 
#' @export
#' 
#' @examples 
#' file <- system.file("extdata",
#'                     "synthetic.exposure.csv",
#'                     package = "ICAMS")
#' exposure <- ReadExposure(file)
#' PlotExposureToPdf(exposure, file = file.path(tempdir(), "exposure.pdf"))
PlotExposureToPdf <- function(exposure,
                              file,  
                              mfrow            = c(2, 1),
                              mar              = c(2, 4, 3, 2),
                              oma              = c(3, 2, 0, 2),
                              samples.per.line = 30,
                              plot.proportion  = FALSE,
                              xlim             = NULL,
                              ylim             = NULL,
                              legend.x         = NULL,
                              legend.y         = NULL,
                              cex.legend       = 0.9,
                              ...
) {
  # Setting the width and length for A4 size plotting
  grDevices::cairo_pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
  
  opar <- par(mfrow = mfrow, mar = mar, oma = oma)
  on.exit(par(opar))
  
  list <- PlotExposure(exposure         = exposure, 
                       samples.per.line = samples.per.line,
                       plot.proportion  = plot.proportion,
                       xlim             = xlim,
                       ylim             = ylim,
                       legend.x         = legend.x,
                       legend.y         = legend.y,
                       cex.legend       = cex.legend,
                       ...              = ...)
  
  grDevices::dev.off()
  invisible(list(plot.success = TRUE, mp.coordinates = list$mp.coordinates))
}