narrowPeaktoGRanges <- function(file) {
  c <- read.table(file,sep="\t",header=FALSE)
  peak <- GRanges(
    seqnames=as.character(c[[1]]),
    ranges=IRanges(start=c[[2]],end=c[[3]])
  )
  mcols(peak) <- c
  return(peak)
}