\name{opt.clust.num}
\alias{opt.clust.num}

\title{
Finding Optimum Number of Cluster for Metagenomics Data
}
\description{
This function will give optimum number of clusters based on Within Sum of Squares (wss) plot.
}
\usage{
opt.clust.num(data, nc, seed = 1234)
}
\arguments{
  \item{data}{
Feature matrix consisting of different genomic features.Each row represents features corresponding to a particular individual or contig and each column represents different genomic features.
}
  \item{nc}{
Probable number of clusters
}
  \item{seed}{
Seed value for iteration
}
}
\value{
WSS plot
}
\author{
Dipro Sinha <\email{diprosinha@gmail.com}>,Sayanti Guha Majumdar, Anu Sharma, Dwijesh Chandra Mishra
}
\examples{
library(metaCluster)
data(metafeatures)
wss_plot <- opt.clust.num(metafeatures[1:200,], nc=10, seed = 1234)
}
\keyword{ WSS }
