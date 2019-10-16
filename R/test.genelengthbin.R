genelengthbin <- function (df, start='trans.start.pos', end='trans.end.pos') {
  df <- data.frame(df)
  df$length <- abs(df[, start] - df[, end])
  df$lengthBin <- NULL
  range <- max(df$length) - min (df$length)
  df$lengthBin[df$length <= 0.25*range] <- 1
  df$lengthBin[df$length <= 0.5*range & df$length > 0.25*range] <- 2
  df$lengthBin[df$length <= 0.75*range & df$length > 0.5*range] <- 3
  df$lengthBin[df$length > 0.75*range] <- 4
}