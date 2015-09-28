

getNeighborhood <- function(comp, side=5, Width=1000) {
  # transcript info
  tx_name <- names(comp)
  granges <- unlist(comp)
  tx <- granges[tx_name]
  strand <- as.character(as.data.frame(strand(tx))[[1]])
  chr <- as.character(as.data.frame(seqnames(tx))[[1]])
  
  if (side == 5) {
    checkpoints <- rep(1,length(tx))
    checkpoints_transcript <- GRanges(seqnames=tx_name,
                                      IRanges(start=checkpoints, end=checkpoints, names=tx_name),
                                      strand=strand)
    # convert to genomic coordinates
    checkPoints_genomic <- mapFromTranscripts(checkpoints_transcript, comp)
    # resize
    checkRegion_genomic <- resize(x=checkPoints_genomic, width=Width, fix="end")
  } else if (side == 3) {
    checkpoints <- sum(width(comp))
    checkpoints_transcript <- GRanges(seqnames=tx_name,
                                      IRanges(start=checkpoints, end=checkpoints, names=tx_name),
                                      strand=strand)
    # convert to genomic coordinates
    checkPoints_genomic <- mapFromTranscripts(checkpoints_transcript, comp)
    # resize
    checkRegion_genomic <- resize(x=checkPoints_genomic, width=Width, fix="start")
  }
  
  # convert to list
  names(checkRegion_genomic) <- rep("",length(tx))
  sidelist <- split(checkRegion_genomic, tx_name, drop=TRUE)
  sidelist <- sidelist[tx_name]
  mapped_chr <- as.character(as.data.frame(seqnames(checkRegion_genomic))[[1]])
  mcols(sidelist) <- data.frame(mapped_chr)
  
  return(sidelist)
}