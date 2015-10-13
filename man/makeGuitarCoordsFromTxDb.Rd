\name{makeGuitarCoordsFromTxDb}
\alias{makeGuitarCoordsFromTxDb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
makeGuitarCoordsFromTxDb
}
\description{
Make a Guitar Coordinates from TranscriptDb object, i.e., making Guitar coordinates for 8 different components, including, 5'UTR, CDS, 3'UTR, lncRNA, Pomoter and Tail of mRNA and lncRNA. Additional filters will discard transcripts that are too short or has too much ambigous on Genome to increase the sensitivity of the analysis.
}
\usage{
makeGuitarCoordsFromTxDb(txdb, maximalAmbiguity = 3, 
    minimalComponentLength = 100, minimalNcRNALength = 300, 
    noBins = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{txdb}{
A transcriptDb object, which can be generated from makeTxDbFromUCSC or other functions.
}
  \item{maximalAmbiguity}{
If a transcript overlap with more number of transcripts than this number, this transcript will be used in the analysis. By filtering out a number of transcripts, this filter also decrease memory usage and 
computation time. Default: 3. 
}
  \item{minimalComponentLength}{
The mimimal length of the components (5'UTR, CDS, 3'UTR) of a mRNA. Unfortunately, some mRNAs do not all 3 components or some components can be too short and cannot provide effective resolution for the analysis. These mRNAs will be filtered out from the analysis. Default: 100
}
  \item{minimalNcRNALength}{
non-coding RNAs with length smaller than this value will not be used in the analysis.
}
  \item{noBins}{
The number of sections the "transcript" of the GRangesList object will be equally divided into. This is the resolution of the Guitar coordinates. The larger noBins is, the clearer the visualization will be; however, more computation time and memory resource will be required. You may want to set a smaller number when using a slow computer. Default: 100.
}
}

\value{
A Guitar coordiantes (GRanges object) will be returned, with Transcript ID, the relative position of each GRanges on the RNA transcript, the interval (bp) between different coordiantes on a transcript component: Front (DNA), Back (DNA), 5'UTR, CDS, 3'UTR, lncRNA. 
}

\author{
Jia Meng <jia.meng@hotmail.com>
}

\examples{
txdb_file <- system.file("extdata", "hg19_knownGene_sample.sqlite",
                         package="GenomicFeatures")
txdb <- loadDb(txdb_file)
gc_txdb <- makeGuitarCoordsFromTxDb(txdb,noBins =3)
gc_txdb
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
