### Plotting the Biotagging 2017 clusters transcription by strand, binned across the span of the clusters

# load libraries
library(PCAtools)
library(data.table)
library(dplyr)
options(scipen=999)

# loading data and formatting matrix
bins <- c(paste(seq(-8,-1,by=1),"plus",sep=""),
          paste(seq(1,8,by=1),"plus",sep=""),
          paste(seq(-8,-1,by=1),"minus",sep=""),
          paste(seq(1,8,by=1),"minus",sep="")) 
            
filelist <- c("ATACCrossRef/Galaxy351-[Biot_cluster1_1.6kb_100bp__Heatmap_values].tabular",
            "ATACCrossRef/Galaxy353-[Biot_cluster2_1.6kb_100bp__Heatmap_values].tabular",
            "ATACCrossRef/Galaxy355-[Biot_cluster3_1.6kb_100bp__Heatmap_values].tabular",
            "ATACCrossRef/Galaxy359-[Biot_cluster4_1.6kb_100bp__Heatmap_values].tabular",
            "ATACCrossRef/Galaxy363-[Biot_cluster5_1.6kb_100bp__Heatmap_values].tabular",
            "ATACCrossRef/Galaxy365-[Biot_cluster6_1.6kb_100bp__Heatmap_values].tabular",
            "ATACCrossRef/Galaxy367-[Biot_cluster7_1.6kb_100bp__Heatmap_values].tabular",
            "ATACCrossRef/Galaxy349-[Biot_cluster8_1.6kb_100bp__Heatmap_values].tabular",
            "ATACCrossRef/Galaxy369-[Biot_cluster9_1.6kb_100bp__Heatmap_values].tabular",
            "ATACCrossRef/Galaxy361-[Biot_cluster10_1.6kb_100bp__Heatmap_values].tabular")

one <- fread(filelist[1],header=FALSE)
df <- one[,1:32]
col <- paste("cluster1_",1:nrow(one),sep="")
for(i in 2:10){
     tab <- fread(filelist[i],header=FALSE)
     df <- rbind(df,tab[,1:32],use.names=FALSE)
     col <- c(col,paste("cluster",i,"_",1:nrow(tab),sep=""))
}

matrix <- as.matrix(df)
matrix <- t(matrix)
colnames(matrix) <- col
row.names(matrix) <- bins    

# metadata for the PCA
colData <- data.frame(cluster=gsub("_\\d+",'', col)) 
row.names(colData) <- col

# get the PCA object
p <- pca(matrix, metadata = colData)
pscale <- pca(matrix, metadata = colData,scale=TRUE)

#plots
screeplot(p, axisLabSize = 18, titleLabSize = 22)
screeplot(pscale, axisLabSize = 18, titleLabSize = 22)

biplot(pscale, showLoadings = TRUE,
    labSize = 1, pointSize = 5, sizeLoadingsNames=5, 
    colby='cluster',legendPosition='right') 

pairsplot(pscale)  

plotloadings(pscale,
    rangeRetain = 0.01,
    labSize = 4.0,
    title = 'Loadings plot',
    subtitle = 'PC1, PC2, PC3, PC4, PC5',
    caption = 'Top 1% variables',
    shape = 24,
    col = c('limegreen', 'black', 'red3'),
    drawConnectors = TRUE)

# convert to prcomp to add new data
pscale.prcomp <- list(sdev = pscale$sdev,
    rotation = data.matrix(pscale$loadings),
x = data.matrix(pscale$rotated),
    center = TRUE, scale = TRUE)

class(pscale.prcomp) <- 'prcomp'

