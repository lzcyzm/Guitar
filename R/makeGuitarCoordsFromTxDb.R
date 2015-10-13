

# make Guitar Coordinates from TranscriptDb object
makeGuitarCoordsFromTxDb <- function(txdb, 
                                     maximalAmbiguity = 3, 
                                     minimalComponentLength = 100,
                                     minimalNcRNALength = 300,
                                     noBins=100){
  
  parameter = list()
  parameter$txdb <- txdb
  parameter$maximalAmbiguity <- maximalAmbiguity # whether overlap with another transcript
  parameter$minimalComponentLength <- minimalComponentLength # minimal length required for each component
  parameter$minimalNcRNALength <- minimalNcRNALength
  parameter$noBins <- noBins
  
  # prepare the bins
  component <- .extractComponent(parameter)
  side <- .get2sides(component)
  
  # print
  print("Building Guitar Coordinates. It may take a few minutes ...")
  utr3_bin <- makeGuitarCoordsFromGRangesList(component[["utr3"]],parameter$noBins)
  utr5_bin <- makeGuitarCoordsFromGRangesList(component[["utr5"]],parameter$noBins)
  cds_bin <- makeGuitarCoordsFromGRangesList(component[["cds"]],parameter$noBins)
  ncRNA_bin <- makeGuitarCoordsFromGRangesList(component[["ncRNA"]],parameter$noBins)
  mRNA_front_bin <- makeGuitarCoordsFromGRangesList(side[["mrna_front"]],parameter$noBins)
  mRNA_back_bin <- makeGuitarCoordsFromGRangesList(side[["mrna_back"]],parameter$noBins)
  ncRNA_front_bin <- makeGuitarCoordsFromGRangesList(side[["ncrna_front"]],parameter$noBins)
  ncRNA_back_bin <- makeGuitarCoordsFromGRangesList(side[["ncrna_back"]],parameter$noBins)
  print("Guitar Coordinates Built ...")
  
  # group together
  mcols(utr3_bin) <- data.frame(mcols(utr3_bin),comp="UTR3",category="mRNA")
  mcols(utr5_bin) <- data.frame(mcols(utr5_bin),comp="UTR5",category="mRNA")
  mcols(cds_bin) <- data.frame(mcols(cds_bin),comp="CDS",category="mRNA")
  mcols(mRNA_front_bin) <- data.frame(mcols(mRNA_front_bin),comp="Front",category="mRNA")
  mcols(mRNA_back_bin) <- data.frame(mcols(mRNA_back_bin),comp="Back",category="mRNA")
  
  mcols(ncRNA_bin) <- data.frame(mcols(ncRNA_bin),comp="lncRNA",category="lncRNA")
  mcols(ncRNA_front_bin) <- data.frame(mcols(ncRNA_front_bin),comp="Front",category="lncRNA")
  mcols(ncRNA_back_bin) <- data.frame(mcols(ncRNA_back_bin),comp="Back",category="lncRNA")
  
  
  GuitarCoords <- suppressWarnings(c(mRNA_front_bin,utr5_bin, cds_bin, utr3_bin, mRNA_back_bin,
                                     ncRNA_front_bin,ncRNA_bin,ncRNA_back_bin))
  
  return(GuitarCoords)}

.extractComponent <- function(parameter){
  
  txdb <- parameter$txdb
  # ambiguity filter
  exons <- exonsBy(txdb, by = "tx",use.names=TRUE)
  noTx <- length(exons)
  print(paste("total",noTx,"transcripts extracted ..."));
  
  temp <- countOverlaps(exons, exons)
  ambiguityFilteredTx <- names(exons[temp < (parameter$maximalAmbiguity+2)])
  noTxLeft <- length(ambiguityFilteredTx)
  print(paste("total",noTxLeft,"transcripts left after ambiguity filter ..."))
  exons <- exons[ambiguityFilteredTx]
  
  
  # extract important components
  cds <- cdsBy(txdb, by = "tx",use.names=TRUE)
  utr5 <- fiveUTRsByTranscript(txdb, use.names=TRUE)
  utr3 <- threeUTRsByTranscript(txdb, use.names=TRUE)
  
  # extract mRNAs
  flag_utr5 <- (sum(width(utr5)) > parameter$minimalComponentLength)
  name_utr5 <- names(utr5)[flag_utr5]
  flag_utr3 <- (sum(width(utr3)) > parameter$minimalComponentLength)
  name_utr3 <- names(utr3)[flag_utr3]
  flag_cds <- (sum(width(cds)) > parameter$minimalComponentLength)
  name_cds <- names(cds)[flag_cds]
  name_mRNA <- intersect(intersect(name_utr5,name_utr3),name_cds)
  name_filtered_mRNA <- intersect(name_mRNA,names(exons))
  cds_filtered <- cds[name_filtered_mRNA]
  utr5_filtered <- utr5[name_filtered_mRNA]
  utr3_filtered <- utr3[name_filtered_mRNA]
  print(paste("total",length(cds_filtered),"mRNAs left after component length filter ..."))
  
  # extract mRNAs
  all_mRNA <- unique(c(names(utr5),names(utr3),names(cds))) 
  name_ncRNA <- setdiff(names(exons),all_mRNA) 
  ncRNA <- exons[name_ncRNA] 
  flag_ncRNA <- 
    (sum(width(ncRNA)) > parameter$minimalComponentLength) & 
    (sum(width(ncRNA)) > parameter$minimalNcRNALength)
  name_ncRNA <- names(ncRNA)[flag_ncRNA]
  ncRNA_filtered <- ncRNA[name_ncRNA]
  print(paste("total",length(ncRNA_filtered),"ncRNAs left after ncRNA length filter ..."))
  
  # return the result
  comp <- list(cds=cds_filtered,utr3=utr3_filtered,utr5=utr5_filtered,ncRNA=ncRNA_filtered)
  return(comp)}

.get2sides <- function(component) {
  utr3 <- component[["utr3"]]
  utr5 <- component[["utr5"]]
  ncrna <- component[["ncRNA"]]
  
  mrna_front <- getNeighborhood(comp=utr5, side=5)
  mrna_back <- getNeighborhood(comp=utr3, side=3)
  ncrna_front <- getNeighborhood(comp=ncrna, side=5)
  ncrna_back <- getNeighborhood(comp=ncrna, side=3)  
  
  result <- list(mrna_front=mrna_front,
                 mrna_back=mrna_back,
                 ncrna_front=ncrna_front,
                 ncrna_back=ncrna_back)
  return(result)
}




