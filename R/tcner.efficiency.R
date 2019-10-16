tcnerEfficiency <- function(annotated.SBS.vcf) {
  
  df <- annotated.SBS.vcf
  df <- subset(df, df$trans.start.pos != "NA")

  df <- data.frame(df)
  df$length <- abs(df[, 'trans.start.pos'] - df[, 'trans.end.pos'])
  df$lengthBin <- NULL
  
  df$lengthBin[df$length <= quantile(df$length,c(0.25,0.5,0.75))[[1]]] <- 1
  df$lengthBin[df$length <= quantile(df$length,c(0.25,0.5,0.75))[[2]] & df$length > quantile(df$length,c(0.25,0.5,0.75))[[1]]] <- 2
  df$lengthBin[df$length <= quantile(df$length,c(0.25,0.5,0.75))[[3]] & df$length > quantile(df$length,c(0.25,0.5,0.75))[[2]]] <- 3
  df$lengthBin[df$length > quantile(df$length,c(0.25,0.5,0.75))[[3]]] <- 4
  
  dfsmall <- df[df$length <= 1e+06, ]
  dflarge <- df[df$lengthBin == 4, ]
  dflargest <- df[df$length > quantile(df$length, 0.9), ]
  
  Plotdist2TSS(dist2TSS(annotated.SBS.vcf = dfsmall), 'Small Genes')
  Plotdist2TSS(dist2TSS(dflarge), 'Large Genes')
  Plotdist2TSS(dist2TSS(dflargest), 'Largest Genes')

}
