\name{clust.suite}
\alias{clust.suite}
\title{
Determination of Suitable Clustering Algorithm for Metagenomics Data
}
\description{
This function will give the best clustering algorithm for a given metagenomics data based on silhouette index for kmeans clustering, kmedoids clustering, fuzzy kmeans clsutering, DBSCAN clustering and hierarchical clsutering. }
\usage{
clust.suite(data, k, eps, minpts)
}
\arguments{
  \item{data}{
Feature matrix consisting of different genomic features.Each row represents features corresponding to a particular individual or contig and each column represents different genomic features.
}
  \item{k}{
Optimum number of clusters
}
  \item{eps}{
Radius value for DBSCAN clustering
}
  \item{minpts}{
Minimum point value of DBSCAN clustering
}
}
\value{
  \item{kmeans }{Output of kmeans clustering}
  \item{kmedoids }{Output of kmedoids clustering}
  \item{fkmeans }{Output of fuzzy kmeans clustering}
  \item{dbscan }{Output of dbscan clustering}
  \item{hierarchical }{Output of hierarchical clustering}
  \item{silhouette.kmeans }{Silhouette plot of kmeans clustering}
  \item{silhouette.kmedoids }{Silhouette plot of kmedoids clustering}
  \item{silhouette.fkmeans }{Silhouette plot of fuzzy kmeans clustering}
  \item{silhouette.dbscan }{Silhouette plot of dbscan clustering}
  \item{silhouette.hierarchical }{Silhouette plot of hierarchical clustering}
  \item{best.clustering.method}{Best clustering algorithm based on silhouette index}
  \item{silhouette.summary }{Average silhouette width of each clustering algorithm}
}
\author{
Dipro Sinha <\email{diprosinha@gmail.com}>,Sayanti Guha Majumdar, Anu Sharma, Dwijesh Chandra Mishra
}
\examples{

library(metaCluster)
data(metafeatures)
result <- clust.suite(metafeatures[1:200,],8,0.5,10)
}

\keyword{ kmeans }
\keyword{ kmedoids }
\keyword{ Fuzzy-kmeans }
\keyword{ DBSCAN }
\keyword{ Hierarchical }
\keyword{ Metagenomics }
\keyword{ Binning }
