\name{oligo.freq}
\alias{oligo.freq}
\title{
Oligonucleotide Frequency
}
\description{
This function will calculate oligonucleotide frequency of each sequence or contig from a FASTA file.
}
\usage{
oligo.freq(fasta_file, f)
}
\arguments{
  \item{fasta_file}{
Name of the fasta or multifasta file
}
  \item{f}{
Length of the oligonucleotide
}
}
\value{
Frequency value of each oligonucleotide of length specified by the user
}
\author{
Dipro Sinha <\email{diprosinha@gmail.com}>,Sayanti Guha Majumdar, Anu Sharma, Dwijesh Chandra Mishra
}
\examples{
library(metaCluster)
freq <- oligo.freq(fasta_file = system.file("extdata/sample1.fasta", package = "metaCluster"),4)
}
\keyword{ oligonucleotide frequency }
