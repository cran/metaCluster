#' @title Finding Optimum Number of Cluster for Metagenomics Data
#'
#' @description This function will give optimum number of clusters based on Within Sum of Squares (wss) plot.
#'
#' @param data
#'
#' @param nc
#'
#' @param seed
#'
#' @return WSS plot
#'
#' @examples
#'
#' @export
  opt.clust.num <- function(data, nc, seed=1234){
  cl_data<-data
  cl_data.features=cl_data
  cl_data.features$class<-NULL
  wss <- (nrow(cl_data.features)-1)*sum(apply(cl_data.features,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(cl_data.features, centers=i)$withinss)}
  wss.plot <- plot(1:nc, wss, type="b", xlab="Number of Clusters",
                   ylab="Within groups sum of squares")
return(wss.plot)
  }
