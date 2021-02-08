if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("zellkonverter")
BiocManager::install("sva")
library(sva)
library(MASS)
library(data.table)
library(dplyr)
library(zellkonverter)

#############
# Obtaining an H5AD file.
sce <- readH5AD("adata.h5ad")
sce
total_counts <- sce@colData@listData[["total_counts"]]
spatial <- sce@int_colData@listData[["reducedDims"]]@listData[["spatial"]]
countmatrix <- sce@assays@data@listData[["X"]] %>% as.matrix()
t_countmatrix <- countmatrix %>% as.data.frame() %>% transpose()
colnames(t_countmatrix) <- rownames(countmatrix)
rownames(t_countmatrix) <- c(1:dim(countmatrix)[2])
#############
cells <- data.frame(X = paste0(spatial[, 1], 'x', spatial[, 2]), x = spatial[, 1], y = spatial[, 2], total_counts = total_counts)
counts <- data.frame(X = paste0(spatial[, 1], 'x', spatial[, 2]), t_countmatrix)
#############
#counts <- read.table(file='/Users/phuong/GPcounts/data/MouseOB/Rep11_MOB_0.csv', header=TRUE, sep=',')
#cells <- read.table(file='/Users/phuong/GPcounts/data/MouseOB/locations.csv', header=TRUE, sep=',')
rownames(counts) = counts[,1]
rownames(cells) = cells[,1]
counts = counts[,-1]
cells = cells[,-1]

counts = cbind(counts,cells$total_counts)
colnames(counts)[dim(counts)[2]] = 'total_counts'
coeffs_nb <-list()
scales_nb <- list()
total_counts = c(cells$total_counts)

for (x in c(1:100)){
  results<-glm.nb(formula=counts[,x]~0+counts$total_counts,link='identity', data=counts)
  coeffs = as.data.frame(results$coefficients)
  scales <- results$coefficients*total_counts
  scales = as.data.frame(scales)
  scales_nb <- append(scales_nb,scales)
}

for (x in c(1:100)){
  results<-glm.nb(formula=counts[,x]~0+counts$total_counts,link='identity', data=counts, start=c(2,0.1,0.1,0.1))
  coeffs = as.data.frame(results$coefficients)
  scales <- results$coefficients*total_counts
  scales = as.data.frame(scales)
  scales_nb <- append(scales_nb,scales)
}

for (x in c(1:100)){
  results<-glm.nb(formula=counts[,x]~0+counts$total_counts,link='identity', data=counts, start=c(14484))
  coeffs = as.data.frame(results$coefficients)
  scales <- results$coefficients*total_counts
  scales = as.data.frame(scales)
  scales_nb <- append(scales_nb,scales)
}

for (x in c(1:100)){
  results<-glm.nb(formula=counts[,x]~0+counts$total_counts,link='identity', data=counts, start=c(36602))
  coeffs = as.data.frame(results$coefficients)
  scales <- results$coefficients*total_counts
  scales = as.data.frame(scales)
  scales_nb <- append(scales_nb,scales)
}

scales_nb=write.table(scales_nb,sep = "\t","scales_nb.txt",col.names = TRUE, row.names = F)
