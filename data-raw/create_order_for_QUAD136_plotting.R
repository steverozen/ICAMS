CreateOrderForQUAD136Plotting <- function() {
  ref.order <- c("AC", "AT", "GC", "CC", "CG", "CT", "TA", "TC", "TG", "TT")
  base <- c("A", "C", "G", "T")
  rev.base <- c("T", "G", "C", "A")
  output <- character(0)
  for (i1 in ref.order) {
    for (i2 in rev.base) {
      for (i3 in base) {
        output <- c(output, paste0(i2, i1, i3))
      }
    }
  }
  return(output)
}

order.for.QUAD136.plotting <- CreateOrderForQUAD136Plotting()
