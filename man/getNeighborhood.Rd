\name{getNeighborhood}
\alias{getNeighborhood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getNeighborhood
}
\description{
This functon extract the neighborhood DNA regions of a GRangesList object. 
The neighborhood DNA regions represent the promoter region (5'end) and complementary 
DNA region on the 3'end, which we call "Tail" here. The 5'end neighborhood DNA can often be useful 
as a kind of negative control region that is not associated with a clear biological function.

This function is different from the "promoters" function available from GenomicRanges R package, which is
designed to extract only the promoter region (5' side DNA), but not the 3' side DNA region.   
}
\usage{
getNeighborhood(comp, side = 5, Width = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{comp}{
a GRangesList object, whose 5' or 3' side DNA neighborhood will be extracted.
}
  \item{side}{
either 5 or 3. When it is 5, the promoter region of the GRangesList object will be extracted; if it is 3, the region on the 3' side complementary to the promoter region will be extracted.  
}
  \item{Width}{
The width of the neighborhood region will be extracted. default: 1000
}
}
\details{
All inner list elements of the input GRangesList object  must have the same strand
}
\value{
A GRangesList, of the same length as the input, will be returned. You may use unlist to convert it into a GRanges object.
}

\author{
Jia Meng <jia.meng@hotmail.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
gr1 <-
  GRanges(seqnames = "chr2", ranges = IRanges(3000, 6000),
          strand = "+", score = 5L, GC = 0.45)
gr2 <-
  GRanges(seqnames = c("chr1", "chr1"),
          ranges = IRanges(c(7000,13000), width = 3),
          strand = c("+"), score = 3:4, GC = c(0.3, 0.5))
grl <- GRangesList("gr1" = gr1, "gr2" = gr2)
promoter <- getNeighborhood(grl,5)
promoter
promoter_GRanges <- unlist(promoter)
promoter_GRanges
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
