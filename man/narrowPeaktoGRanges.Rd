\name{narrowPeaktoGRanges}
\alias{narrowPeaktoGRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
narrowPeaktoGRanges
}
\description{
read the narrowpeak format from MACS software into GRangesList object.

NAME_peaks.narrowPeak is BED6+4 format file which contains the peak locations together with peak summit, 
pvalue and qvalue. You can load it to UCSC genome browser. Definition of some specific columns are:
5th: integer score for display
7th: fold-change
8th: -log10pvalue
9th: -log10qvalue
10th: relative summit position to peak start.
}
\usage{
narrowPeaktoGRanges(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ a string specifies where the narrow peak file is located.
}
}

\value{
A GRanges object will be returned.
}
\references{
MACS2 software: https://github.com/taoliu/MACS/tree/master/MACS2
}

\examples{
narrowPeak=system.file("extdata", "m6A_hg19_1000peaks_macs2.narrowPeak", package="Guitar")
m6A_Bcell <- narrowPeaktoGRanges(narrowPeak) # bam imported as GAlignments
m6A_Bcell
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

