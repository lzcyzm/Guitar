

makeGuitarCoordsFromGRangesList <- function(comp, 
                                            noBins=100,
                                            collapseGene=FALSE,
                                            width = 51) {
  
  # get all the check points
  tx_length <- as.numeric(sum(width(comp)))
  checkpoints_interval <- tx_length/noBins
  
  # get transcript name and strand
  tx_name <- names(comp)
  granges <- unlist(comp)
  tx <- granges[tx_name]
  strand <- as.character(as.data.frame(strand(tx))[[1]])
  chr <- as.character(as.data.frame(seqnames(tx))[[1]])
  
  # get coordinates
  t <- lapply(X=1:noBins, 
              FUN=.makeGuitarCoordsForSingleIndex, 
              comp=comp,
              noBins=noBins,
              checkpoints_interval=checkpoints_interval,
              strand=strand,
              chr=chr,
              tx_name=tx_name)  
  
  GuitarCoords <- .combineListOfGRanges(t)
  
  if (collapseGene) {
    temp <- split(GuitarCoords, mcols(GuitarCoords)$pos, drop=FALSE)
    mcols(temp) <- data.frame(pos=unique(mcols(GuitarCoords)$pos))
    GuitarCoords <- temp
  }
  
  # return the result
  return(GuitarCoords)
}

.combineListOfGRanges <- function(t){
  txt <- "c(t[[1]]"
  for (i in 2:length(t)) {
    txt <- paste(txt,",t[[",i,"]]",sep="")
  }
  txt <- paste(txt,")",sep="")
  c <- parse(text=txt)
  
  # suppressWarnings
  # GuitarCoords <- eval(c)
  GuitarCoords <- suppressWarnings(eval(c))
  
  return(GuitarCoords)
}
.makeGuitarCoordsForSingleIndex <- function(
  index, comp, noBins, checkpoints_interval,
  strand, chr, tx_name,
  binWidth=51) {
  
  # get checkpoints
  checkpoints <- round(checkpoints_interval*(index-0.5))
  checkpoints_transcript <- GRanges(seqnames=tx_name,
                                    IRanges(start=checkpoints, end=checkpoints, names=tx_name),
                                    strand=strand)
  
  # convert to genomic coordinates
  checkPoints_genomic <- mapFromTranscripts(checkpoints_transcript, comp)
  
  # resize
  # binWidth <- 4*round(checkpoints_interval)+1
  checkRegion_genomic <- resize(x=checkPoints_genomic, 
                                width=binWidth, 
                                fix="center")
  

  # add annotation information
  mcols(checkRegion_genomic) <- data.frame(txid=tx_name, pos=(index-0.5)/noBins, interval=checkpoints_interval)
  
  return(checkRegion_genomic)
}