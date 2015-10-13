\name{BED12toGRangesList}
\alias{BED12toGRangesList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BED12toGRangesList
}
\description{
read bed12 format into R as GRangesList object
}
\usage{
BED12toGRangesList(filepath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepath}{
the path where the bed12 file is located
}
}
\details{
The function read bed12 into R as GRangesList object, with the introns spliced out, making it different from other functions such as import.bed.
}
\value{
The returned GRangsList object has the same number of GRangesLists as the number of rows in the bed12 file.
}
\references{
For more infomration about the bed format, please refer to: 
https://genome.ucsc.edu/FAQ/FAQformat.html#format1
}
\author{
Jia Meng <jia.meng@hotmail.com> 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
bed12=system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed", package="Guitar")
m6A_HepG2 <- BED12toGRangesList(bed12) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bed12 }
\keyword{ GRangesList }% __ONLY ONE__ keyword per line
