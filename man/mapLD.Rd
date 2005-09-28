\name{mapLD}
\alias{mapLD}
\title{Construct Haplotype Blocks}
\description{
  mapLD calculates confidence intervals for Lewontin's D' (1964) and
  constructs haplotype blocks using the approach described by Gabriel, et
  al (2002).
}
\usage{
mapLD(SNPdata, locusID.col, subjectID.col, allele.cols, WhichGene = NA, outgraph = NA)
}
\arguments{
  \item{SNPdata}{Data frame containing SNP data. At least 4 fields are
    required, see arguments \code{locusID.col}, \code{subjectID.col} and
  \code{allele.cols} for details.}
  \item{locusID.col}{Column name or index for a marker ID.}
  \item{subjectID.col}{Column name or index for a sample ID.}
  \item{allele.cols}{A vector of 2 for the column names or indices for
    the values of two alleles.}
  \item{WhichGene}{Gene name. Default is NA.}
  \item{outgraph}{Name of an EPS file to which the heatmap of LD is
    printed. If not specified, the heatmap is printed to a file called
    LDmap.eps under the working directory where the current R session is
  running.}
}
\details{
  EM algorithm is not used for estimating two-locus haplotype
  frequencies. Instead, A faster one-dimensional golden section search
  combined with parabolic interpolation is applied to finding the MLE for
  two-locus haplotype frequencies (Weir BS (1996)). In addition, a
  recursive function is defined to compute the haplotype blocks
  according to Gabriel et al's (2002) definition later on
  modified by Wall and Pritchard (2003).
}
\value{
  A list of the following components:
  \item{LDinfo}{A dataframe containing information about the point
    estimate and the corresponding 90\% coverage confidence interval
    for the pair-wise D', as well as the four 2-locus haplotype freqiencies.}
  \item{LDblock}{A dataframe containing the head and tail SNPs for each
    of the haplotype blocks.}
  \item{LocusIndex}{A numeric index for each of the SNP markers under
    investigation.}

  
  In addition, an pair-wise LD heatmap with haplotype block boundaries
  labelled is printed to a postscript file. See argument \code{outgraph}
  for details.
}
\references{1. Gabriel SB, et al (2002). The Structure of
  Haplotypes in the Human Genome. Science, 296(5576):2225-9

  
  2.Wall JD and Pritchard JK (2003). Assessing the Performance
  of the Haplotype Block Model of Linkage Disequilibrium. Am J Hum
  Genet., 73(3):502-15


  3. Lewontin RC (1964). The interaction of selection and
  linkage. I. General considerations: heterotic models. Genetics,
  49:49-67.


  4. Weir BS (1996). Genetic Data Analysis II. Sinauer, Sunderland, MA}
  
\author{Peter Hu and Jared Lunceford}
\note{A LD map is most meaningful when the markers are ordered by
  physical locations on the chromosome. It is therefore recommended markers
  be sorted before data are passed into the \code{mapLD} function.}

 \section{Warning}{If the value of sample ID is an integer exceeding
   the precision setup in R for integers, then merging two loci marker data by
   sample ID is error-prone. Such a problem could be avoided by specifying
   option colClasses = 'character' when data are read using function
   \code{read.table}.} 

 \examples{
data(SNPdata)
getLD <- mapLD(SNPdata = SNPdata,
locusID.col = 'markerID',
subjectID.col = 'subjectID',
allele.cols = 1:2,
WhichGene = NA,
outgraph = NA)
}
\keyword{ hplot }
\keyword{ methods }
