\name{SNPdata}
\alias{SNPdata}
\docType{data}
\title{An Example of Genotyping Data}
\description{
   The example SNP data are a subset of real SNP data. Due to the
   concerns on confidentiality, Real subject IDs and marker IDs are
   replaced with their corresponding indices.
}
\usage{data(SNPdata)}
\format{
  A data frame with 100 observations and 13 SNP markers on the following 5 variables:
  \describe{
    \item{Allele1:}{Allele variant 1}
    \item{Allele2:}{Allele variant 2}
    \item{subjectID:}{Sample ID}
    \item{markerID:}{SNP Marker ID}
    \item{position:}{Relative position of the SNP marker in the chromosome}
  }
}
\examples{
data(SNPdata)
}
\keyword{datasets}
