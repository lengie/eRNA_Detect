### tSNEBidirReg.r
### t-SNE and subsequent clustering on raw bidirectional regions
###
### Last updated: July 2021
### Input: tabular data from createMatrix in deepTools, binning transcription across plus and minus strands of bidirectional regions identified through custom code
### (Current) Output: t-SNE plot, heirarchical and k-menas clustering of the t-SNE in plots

library(gridExtra)
library(ggplot2)
library(dplyr) 
library(data.table)
library(parallel)
library(stringr)
library(factoextra)
library(Rtsne)
library(tictoc)
library(fastcluster)
library(cluster)
options(scipen=999)

load <- fread("sox10_Zv9BidirReprod100Rem500bpNearestSpannedUnder1500bpClusterLabeledBinnedTranscription.csv")
smallbins <- fread("Galaxy415-[100Rem500bpNearestSpanned_under1500bp_1.5kb_50bp__Heatmap_values].tabular")
newunder15 <- as.matrix(smallbins[,1:60])
row.names(newunder15) <- colnames(load)[2:length(load)]
colnames(newunder15) <- c(paste(rep(-7:-1,each=2),c("50bp_plus","00bp_plus"),sep=""),"-50bp_plus",
                          paste(rep(1:7,each=2),c("00bp_plus","50bp_plus"),sep=""),
                          paste(rep(-7:-1,each=2),c("50bp_minus","00bp_minus"),sep=""),"-50bp_minus",
                          paste(rep(0:7,each=2),c("00bp_minus","50bp_minus"),sep=""))   #bins are a bit off at the 0 point

noZ <- newunder15[rowSums(newunder15 != 0) > 0,]
data <- as.data.frame(noZ)
scaledata <- scale(noZ)
scaledata <- as.data.frame(scaledata)

clusterno <- row.names(noZ)
clusterno <- replace(clusterno,grepl("Bidir",clusterno),11)
clusterno <- gsub('cluster([0-9]+)_.*','\\1',clusterno)
geom <- row.names(noZ)
geom <- gsub('cluster([0-9]+)_.*','\\1',geom) 
geom <- stringr::str_replace_all(geom,setNames(c(rep("3",5),rep("2",5)),as.character(c(10:1))))
geom <- replace(geom,grepl("Bidir",geom),1)
geom <- as.factor(as.numeric(geom))
data$clusterno <- as.factor(clusterno)
scaledata$clusterno <- as.factor(clusterno)

unique <- unique(scaledata)
set.seed(333)
ind <- sample(dim(unique)[1], 600000)
samp <- unique[ind,]
        

tic()
tsnescaled_100 <- Rtsne(samp[,1:60], dims = 2, perplexity=100, verbose=TRUE, max_iter = 100,check_duplicates=FALSE,num_threads=0)
toc()
d_tsne_100 = as.data.frame(tsnescaled_100$Y)
write.table(d_tsne_100,"sox10_Zv9BidirReprod100Rem500bpNearestSpannedUnder1500bp100bp_tsne600krand333_100perp.csv",quote=FALSE,sep='\t')
summary(tsnescaled_100)
summary(d_tsne_100)

png("tsne_under1500bp100bin_seed333rand600k_100perp.png",width=1080,height=1080)
ggplot(d_tsne_100, aes(x=V1, y=V2,color=as.factor(samp[,61]))) +
  geom_point(size=0.25) +
  scale_color_manual(name="Orig. clusters",values=c("black","purple","gray","gold","orange","red","magenta","green","green4","cyan","cornflowerblue")) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_light(base_size=20) +    
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) 
dev.off()

#6 clusters? 9? Looking back at the k-means elbows, looks like 6 is best for scaled data (though that was pre-tSNE remap)
## Creating k-means clustering model, and assigning the result to the data used to create the tsne
fit_cluster_kmeans=kmeans(scale(d_tsne_100), 6)
d_tsne_100$cl_kmeans = factor(fit_cluster_kmeans$cluster)

set.seed(111)
postind <- sample(dim(d_tsne_100)[1], 150000)
forplotting <- d_tsne_100[postind,]
postsamp <- data.matrix(d_tsne_100[postind, c(1,2)])

fit_cluster <- clara(scale(postsamp),k=6,metric="euclidean") #vs manhattan or jaccard distance? Maybe jaccard
#fit_cluster_jaccard <- clara(scale(postsamp),k=6,metric="jaccard")
fit_cluster_hierarchical <- fastcluster::hclust(dist(scale(postsamp)),method="centroid")
  
forplotting$cl_kmeans = factor(fit_cluster$cluster)
#forplotting$cl_kmeansjaccard = factor(fit_cluster_jaccard$cluster)
forplotting$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=6))

#plotting by new kmeans clusters, ignoring old groupings 
plot_cluster=function(data, var_cluster, palette)
{
  ggplot(data, aes_string(x="V2", y="V3", color=var_cluster)) +
  geom_point(size=0.25) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "horizontal") + 
    scale_colour_brewer(palette = palette) 
}

plot_k=plot_cluster(forplotting, "cl_kmeans", "Accent")
plot_h=plot_cluster(forplotting, "cl_hierarchical", "RdYlGn")
## and finally: putting the plots side by side with gridExtra lib...
png("tsne_under1500bp100bin_seed222rand600k_100perp_kmeans_seed111rand150kpoints.png",width=1080,height=1080)
    grid.arrange(plot_k, plot_h,  ncol=2)
dev.off()

#all the points into a k-means
d_tsne_100 <- fread("sox10_Zv9BidirReprod100Rem500bpNearestSpannedUnder1500bp100bp_tsne600krand222_100perp.csv",header=FALSE)
d_tsne_100 <- d_tsne_100[,c(2,3)]
all <- data.matrix(d_tsne_100)
fit_cluster <- clara(scale(all),k=6,metric="euclidean")
d_tsne_100$orig_clusters <- factor(as.numeric(unique[ind,61]))
d_tsne_100$euc_kmeansclusters <- factor(fit_cluster$cluster)

plot_o <- ggplot(d_tsne_100, aes_string(x="V2", y="V3", color="orig_clusters")) +
  geom_point(size=0.25) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "horizontal") +
  scale_color_manual(name="Orig. clusters",values=c("black","purple","gray","gold","orange","red","magenta","green","green4","cyan","cornflowerblue"))
plot_k <- plot_cluster(d_tsne_100,"euc_kmeansclusters","RdYlGn") 

png("tsne_under1500bp100bin_seed222rand600k_100perp_kmeansclara_allpoints.png",width=1920,height=1080)
    grid.arrange(plot_o, plot_k,  ncol=2)
dev.off()
