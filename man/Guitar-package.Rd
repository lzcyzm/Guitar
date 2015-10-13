\name{Guitar-package}
\alias{Guitar-package}
\alias{Guitar}
\docType{package}
\title{
Guitar
}
\description{
RNA Landmarks Guided Transcriptomic View of RNA-related Genomic Features. 
}
\details{
The package is designed for transcriptomic visualization of RNA-related genomic features represented with genome-based coordinates with respect to the landmarks of RNA transcripts, i.e.,	transcription starting site, start codon, stop codon and transcription ending site. 
}
\author{
Jia Meng <jia.meng@hotmail.com>
}
\references{
~~ Literature or other references for background information ~~
}

\examples{

# read genomic features
narrowPeak <- system.file(
  "extdata", "m6A_hg19_1000peaks_macs2.narrowPeak", 
  package="Guitar")

# bam imported as GAlignments
m6A_Bcell <- narrowPeaktoGRanges(narrowPeak) 

# generate a list of genomic features
m6A_Bcell_1 <- m6A_Bcell[1:300]
m6A_Bcell_2 <- m6A_Bcell[301:600]
m6A_Bcell_3 <- m6A_Bcell[601:900]
feature_hg19 <- list(m6A_Bcell_1, m6A_Bcell_2, m6A_Bcell_3) 
names(feature_hg19) <- c("m6A_1","m6A_2","m6A_3")

# Make Guitar coordiantes
txdb_file <- system.file("extdata", "hg19_toy.sqlite", 
                         package="Guitar")
txdb <- loadDb(txdb_file)
gc_txdb <- makeGuitarCoordsFromTxDb(txdb,noBins =10)


# Plot
GuitarPlot(feature_hg19, 
           GuitarCoordsFromTxDb = gc_txdb)

}
