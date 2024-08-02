
label<-read.table("AT2_label.txt",header = T)
data<-readRDS("AT2.rds")
label<-label$Sub
library(dplyr)
library(scCAN)
library(umap)
data<-t(data)
set.seed(123)
library(Rtsne)
rtsne_out <- Rtsne(as.matrix(data),
                   dims = 3,
                   initial_dims=50,
                   perplexity=10,
                   pca = T
)
X1 <- rtsne_out$Y[,1:2]
plot(X1)
pca_out <- prcomp(data)
X2 <- pca_out$x[,1:2]
umap_result <- umap(data)
X3<-umap_result$layout
latent<-X1
data1 <- list()

data1$all.latent <- list(X1,X2,X3)

library(SNFtool)
all.sim <- list()
for(i in 1 : length(data1$all.latent)){
  dat <- data1$all.latent[[i]]
  tmp <- dist2(dat,dat)^(1/2)
  all.sim[[i]] <- tmp
}
all.aff <- list()
for(i in 1 : length(all.sim)){
  dat <- all.sim[[i]]
  tmp <- affinityMatrix(dat)
  all.aff[[i]] <- tmp
}

K<-20
T<-10
W = SNF(all.aff, K, T)
res = estimateNumberOfClustersGivenGraph(W, NUMC=3:15)
k1 <- res$`Eigen-gap best`
k2 =res$`Eigen-gap 2nd best`
k <- min(k1,k2)
groups = spectralClustering(W,k)
allCluster<-list(groups)
idx = which(label == "E14")
start.idx<-idx
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
data2<-list()
data2$expr<-data
cluster_dist <- function(latent, data, start.idx) {
  cl <- data$cl
  start.clus = getmode(cl[start.idx])
  
  latent <- cbind(latent, cl)
  latent <- as.data.frame(latent)
  cell_stages <- factor(cl)
  
  x <- latent
  total_col_num <- ncol(x)
  # Euclidian center
  clus_cent <- NULL
  # Select the last column, cell type/clustering result, which has to be represent by numeric id.
  for (i in 1:length(unique(x[, total_col_num]))) {
    temp <- x[which(x[, total_col_num] == i), ]
    temp_cent <- apply(temp, 2, mean)
    clus_cent <- rbind(clus_cent, temp_cent)
  }
  
  clus_cent_row_names <- paste("clus", seq(1:length(unique(x[, total_col_num]))))
  clus_cent <- rbind(clus_cent[start.clus, ], clus_cent[-start.clus, ])
  clus_cent_row_names <- c(clus_cent_row_names[start.clus], clus_cent_row_names[-start.clus])
  rownames(clus_cent) <- clus_cent_row_names
  
  clus_dist <- as.matrix(stats::dist(clus_cent[, 1:(total_col_num - 1)], method = "euclidian"))[, 1]
  clus_cent <- as.data.frame(clus_cent)
  
  dist_order <- NULL
  for (i in 1:length(latent[, total_col_num])) {
    dist_order <- c(dist_order, clus_dist[which(clus_cent$cl == latent[, total_col_num][i])])
  }
  
  out <- list(center = clus_cent, clus_dist = clus_dist, dist_order = dist_order)
}

pseudotime_list <- lapply(allCluster, function(x) rep(0, length(x)))

for (i in seq_along(allCluster)) {
  data2$cl <- allCluster[[i]]
  out <- cluster_dist(latent, data2, start.idx)
  tmp <- out$dist_order %>% as.numeric()
  tmp <- tmp / max(tmp)
  
  for (j in seq_along(allCluster[[1]])) {
    pseudotime_list[[1]][j] <- tmp[j]
  }
}
pseudotime <- lapply(seq_along(pseudotime_list[[1]]), function(j) {
  sum(vapply(seq_along(pseudotime_list), function(i) pseudotime_list[[i]][j], numeric(1)))
})
pseudotime <-unlist(pseudotime)
