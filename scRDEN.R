library(Hmisc)
library(ESEA)
library(ggplot2)
# Load your gene expression matrix 
data<-read.csv("AT2.csv",check.names = FALSE,row.names = 1)
data<-as.matrix(data)
data<-getFilteredData(data,min.cells = 0.1*ncol(data))
max(data)
mu <- apply(data, 1, mean)
sigma <- apply(data, 1, sd)
lower <- mu - 4*sigma

upper <- mu + 4*sigma

data <- ifelse(data > upper, upper, data)
data <- ifelse(data < lower, lower, data)

library(IEntropy)
library(rsvd)
max(data)
data<-log(data+1,base = 2)
int_res <- Get_entropy(data,5,1)
ie<-int_res$Gene
data<-data[which(rownames(data) %in% ie),]

gene_expression_matrix<-data
getFilteredData <- function(data, min.cells = 0.05*ncol(data)) {
  message("Running ...")
  filt.data <- data[Matrix::rowSums(data > 0) > min.cells, ]
  message("Expression Filtering Done.")

  filt.data <- filt.data[!grepl(pattern = "^MT-|^ERCC-|^RPS|^RPL", x = rownames(filt.data), ignore.case = TRUE), ]
  
  message("Mitochondrial, Ribosomal and Pseudo Genes Filtering Done.")
  
  return(filt.data)
}
# 1.compute GGI network
log_gene_expression_matrix<-t(gene_expression_matrix)
# Calculate Pearson correlation coefficients
correlation_matrix <- rcorr(log_gene_expression_matrix)
library(WGCNA)
th<-pickHardThreshold.fromSimilarity(
  abs(correlation_matrix$r),
  RsquaredCut = 0.8, 
  cutVector = seq(0.1, 0.9, by = 0.05),
  moreNetworkConcepts=FALSE, 
  removeFirst = FALSE, nBreaks = 10)
significant_correlations <- abs(correlation_matrix$r) * (correlation_matrix$P <= 0.01)

threshold <- th$cutEstimate

GGI_network <- abs(significant_correlations) >= threshold

removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
GGI_network2<-removeRowsAllNa(GGI_network)
indices <- which(GGI_network2==TRUE, arr.ind = TRUE)
result=data.frame(V1=rownames(GGI_network2)[indices[,1]],
                  V2=colnames(GGI_network2)[indices[,2]])
result<-as.matrix(result)


library(tidyverse)
tmp<-ifelse(result[,1]>result[,2],paste0(result[,1],result[,2]),paste0(result[,2],result[,1]))

out <- result[!(duplicated(tmp) ), ]
new_data<-gene_expression_matrix
rank.matrix <- function(x){
  rankmatrix <- apply(x, 2, function(column) rank(column, ties.method = "average"))
  colnames(rankmatrix) = colnames(x)
  row.names(rankmatrix) = row.names(x)
  return(rankmatrix)
}

delta.rank <- function(rank.data, net) {
  index1 <- match(net[, 1], rownames(rank.data))
  index2 <- match(net[, 2], rownames(rank.data))
  deltarank <- rank.data[index1, ] - rank.data[index2, ]
  rownames <- paste(net[, 1], net[, 2], sep = "|")
  row.names(deltarank) <- rownames
  colnames(deltarank) <- colnames(rank.data)
  return(deltarank)}
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}

rank_data=rank.matrix(new_data)
net<-out
delta_data<-delta.rank(rank_data,net)
norm_data <-dynutils::scale_minmax(delta_data)

find_hv_edges = function(count, I, J){
  count_nzero = lapply(1:I, function(i) setdiff(count[i, ], log10(1.01)))
  mu = sapply(count_nzero, mean)
  mu[is.na(mu)] = 0
  sd = sapply(count_nzero, sd)
  sd[is.na(sd)] = 0
  cv = sd/mu
  cv[is.na(cv)] = 0
  # sum(mu >= 1 & cv >= quantile(cv, 0.25), na.rm = TRUE)
  high_var_genes = which( cv >= quantile(cv, 0.99))
  if(length(high_var_genes) < 500){ 
    high_var_genes = 1:I}
  count_hv = count[high_var_genes, ]
  return(count_hv)
}
I = nrow(norm_data)
J = ncol(norm_data)
count_hv<-find_hv_edges(norm_data,I,J)
count_hv[1:5,1:5]
saveRDS(count_hv,"AT2.rds")




