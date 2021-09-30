#' @title Determination of Suitable Clustering Algorithm for Metagenomics Data
#'
#' @description This function will give the best clustering algorithm for a given metagenomics data based on silhouette index for kmeans clustering, kmedoids clustering, fuzzy kmeans clsutering, DBSCAN clustering and hierarchical clsutering.
#'
#' @param data
#'
#' @param k
#'
#' @param eps
#'
#' @param minpts
#'
#' @return Best clustering algorithm
#'
#' @examples
#'
#' @export
clust.suite <- function(data,k,eps, minpts){
##########################################

#load_library
  requireNamespace("factoextra")
  requireNamespace("cluster")
  requireNamespace("dbscan")
  requireNamespace("dplyr")
  requireNamespace("grDevices")
#########################################
#pre-processing
cl_data<- data
cl_data.features=cl_data[, -1]
cl.scaled <- scale(cl_data[, -1])
#########################################
#k-means clustering
res_km<-kmeans(cl_data.features, k, iter.max = 10, nstart = 1)
sil_km <- silhouette(res_km$cluster, dist(cl.scaled))
plot_kmeans <- fviz_silhouette(sil_km, xlab = "Silhouette Plot of kmeans Clustering")
#########################################

#k-medoids clustering
res_kmo<-pam(cl_data.features, k)
sil_kmo <- silhouette(res_kmo$cluster, dist(cl.scaled))
plot_kmedoids <- fviz_silhouette(sil_kmo, xlab = "Silhouette Plot of kmedoids Clustering")
#########################################
#FKM clustering
res_fkm<-fanny(cl_data.features, k, metric = "euclidean", stand = FALSE)
sil_fkm <- silhouette(res_fkm$cluster, dist(cl.scaled))
plot_fkmean <- fviz_silhouette(sil_fkm, xlab = "Silhouette Plot of Fuzzy kmeans Clustering")
#########################################
#DBSCAN Clustering
res_db <- dbscan(cl_data.features, eps, minpts)
sil_db <- silhouette(res_db$cluster, dist(cl.scaled))
plot_dbscan <- fviz_silhouette(sil_db, xlab = "Silhouette Plot of DBSCAN Clustering")
########################################
#Hierarchical Clustering
distance_mat <- dist(cl_data.features, method = 'euclidean')
res_hc <- hclust(distance_mat, method = "average")
cut_avg <- cutree(res_hc, k)
res_hierarchical <- mutate(cl_data, cluster = cut_avg)
sil_hc <- silhouette(res_hierarchical$cluster, dist(cl.scaled))
plot_hc <- fviz_silhouette(sil_hc, xlab = "Silhouette Plot of Hierarchical Clustering")
#########################################

kmeans <- mean(plot_kmeans$data$sil_width)
kmedoids <- mean(plot_kmedoids$data$sil_width)
fuzzy_kmeans <- mean(plot_fkmean$data$sil_width)
db_scan <- mean(plot_dbscan$data$sil_width)
hierarchical <- mean(plot_hc$data$sil_width)
method_names <- c("kmeans","kmedoids","fuzzy_kmeans","db_scan","hierarchical")
avg.silhouette.width <- c(kmeans,kmedoids,fuzzy_kmeans,db_scan,hierarchical)
final_meth <- cbind(method_names,avg.silhouette.width)
best_method = final_meth[which(final_meth[,2] == max(final_meth[,2]))]

result <- list(kmeans=res_km,kmedoids = res_kmo,fkmeans = res_fkm,dbscan = res_db,hierarchical = res_hierarchical, silhouette.kmeans = plot_kmeans,silhouette.kmedoids = plot_kmedoids,silhouette.fkmeans = plot_fkmean,silhouette.dbscan = plot_dbscan,silhouette.hierarchical = plot_hc,best.clustering.method = best_method,silhouette.summary=final_meth)
return(result)

}
