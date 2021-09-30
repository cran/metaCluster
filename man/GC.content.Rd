\name{GC.content}
\alias{GC.content}
\title{
Calculation of GC content
}
\description{
This function will calculate GC content from each sequence or contigs of a FASTA file.
}
\usage{
GC.content(fasta_file)
}
\arguments{
  \item{fasta_file}{
Name of the fasta or multifasta file
}
}
\value{
Value of the GC content of each sequence or contig.
}
\author{
Dipro Sinha <\email{diprosinha@gmail.com}>,Sayanti Guha Majumdar, Anu Sharma, Dwijesh Chandra Mishra
}
\examples{
library(metaCluster)
library(seqinr)
sample_data <- read.fasta(file = system.file("extdata/sample1.fasta", package = "metaCluster"),
seqtype = "DNA")
gc <- GC.content(sample_data)
}
\keyword{ GC content }
