\name{makeGuitarCoordsFromGRangesList}
\alias{makeGuitarCoordsFromGRangesList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
makeGuitarCoordsFromGRangesList
}
\description{
make Guitar Coordinates From a GRangesList object
}
\usage{
makeGuitarCoordsFromGRangesList(comp, noBins = 100, collapseGene = FALSE, width = 51)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{comp}{
A GRangesList object
}
  \item{noBins}{
The number of sections the "transcript" of the GRangesList object will be equally divided into. This is the resolution of the Guitar coordinates. The larger noBins is, the clearer the visualization will be; however, more computation time and memory resource will be required. You may want to set a smaller number when using a slow computer. Default: 100. 
}
  \item{collapseGene}{
Whether merge the Guitar coordinates from different transcripts. Default: FALSE. 
}
  \item{width}{
The width of each check points. Default: 51
}
}

\value{
Guitar coordinates will be generated, which provides a fast reference between genomic coordinates and the transcriptomic coordinates. By default, the returned value should be in GRanges format; if collapseGene = TRUE, GRangesList object will be returned.
}

\author{
Jia Meng <jia.meng@hotmail.com>
}

\examples{
gr1 <-  GRanges(seqnames = "chr2", ranges = IRanges(3000, 6000),
          strand = "+", score = 5L, GC = 0.45)
gr2 <-  GRanges(seqnames = c("chr1", "chr1"),
          ranges = IRanges(c(7000,13000), width = 3),
          strand = c("+"), score = 3:4, GC = c(0.3, 0.5))
grl <- GRangesList("gr1" = gr1, "gr2" = gr2)
tc <- makeGuitarCoordsFromGRangesList(grl, noBins = 5)
tc
mcols(tc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Guitar }
\keyword{ GRanges }
\keyword{ GRangesList }% __ONLY ONE__ keyword per line
