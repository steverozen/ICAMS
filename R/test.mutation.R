library(data.table)
foo <- output
dt <- data.table(foo)
class <- colnames(dt)[c(4:6, 10:12)]
revc1 <- function(string) {
  return(paste0(ICAMS:::revc(substr(string, 1, 1)), '>', 
                             ICAMS:::revc(substr(string, 3, 3))))
}


PlotStrandBiasTss <- function(matrix, class) {
  
  stopifnot(class %in% c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
  dt <- data.table(matrix)
  revc1 <- function(string) {
    return(paste0(ICAMS:::revc(substr(string, 1, 1)), '>', 
                  ICAMS:::revc(substr(string, 3, 3))))
  }
  output <- dt[, c(revc1(class), class), with = FALSE]
  colnames(output) <- c("antiSense", "Sense")
  barplot(t(output), beside=TRUE, ylim = c(0, max(output) + 10), 
          col = c("grey80", "grey32"), cex.axis = 0.8, cex.names = 0.4, ylab="# of mutations",
          names.arg=paste0('group', 1:10))
  legend("topright", c("antiSense", "Sense"), bty="n", 
         fill=c("grey80", "grey32"))
  title(c(class, "mutations"))
  return(output)
}

PlotStrandBiasTssAll <- function(matrix) {
  
  dt <- data.table(matrix)
  revc1 <- function(string) {
    return(paste0(ICAMS:::revc(substr(string, 1, 1)), '>', 
                  ICAMS:::revc(substr(string, 3, 3))))
  }

  for (class in c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")) {
    output <- dt[, c(revc1(class), class), with = FALSE]
    barplot(t(output), beside=TRUE, ylim = c(0, max(output) *1.2), 
            col = c("grey80", "grey32"), cex.axis = 0.8, cex.names = 0.4, ylab="# of mutations",
            names.arg=paste0('group', 1:10))
    legend("topright", c("antiSense", "Sense"), bty="n", 
           fill=c("grey80", "grey32"))
    title(c(class, "mutations"))
  }
}
