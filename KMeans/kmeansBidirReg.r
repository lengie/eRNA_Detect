### kmeansBidirReg.r
###
### Purpose: Take a bed file of regions of interest, as well as the createMatrix tabular data of transcription across bins spanning the ROIs, 
###          and output k-means clustering of the transcriptional profile 
###
###
### Written by Liana Engie
### Last updated: June 2021
###
### 
### Input:
### (Current) Output: 

library(GenomicFeatures)
library(GenomicAlignments) 
library(ggplot2)
library(dplyr) 
library(data.table)
library(pvclust)
library(parallel)
library(dendextend)
library(stringr)
library(factoextra)
library(foreach)
library(doMC)
library(Rtsne)
options(scipen=999)

##regions under 1.5kb
# createMatrix data of transcription in bins spanning length of genomic location
under <- fread("Galaxy381-[100Rem500bpNearestSpanned_under1500bp_1.4kb_100bp__Heatmap_values].tabular",header=FALSE)
under15 <- t(under) %>% replace(is.na(.),0)

#1400bp centered on middle of region
binsun <- c(paste(seq(-7,-1,by=1),"00bp_plus",sep=""),
          paste(seq(1,7,by=1),"00bp_plus",sep=""),
          paste(seq(-7,-1,by=1),"00bp_minus",sep=""),
          paste(seq(1,7,by=1),"00bp_minus",sep="")) 
row.names(under15) <- binsun

# labeling by number to track
colnames(under15) <- paste("BidirReg",1:dim(under15)[2],sep="")
                      
# bed file of actual ROI genomic locations, minus one that Galaxy/deepTools had issue with
onefive <- fread("sox10_Zv9FlankedWiderNoncodingReprod100Rem500bpNearestSpanned<1.5kb.bed")
colnames(onefive) <- c("chr","start","end","ID","score","strand")
onefiveless <- onefive[-1356503,]
underGR <- GRanges(onefiveless)

# cluster list from sox10 paper
sox10_1 <- fread('ATACCrossRef/Sox10nuclear_cluster1.counts.txt')
sox10_2 <- fread('ATACCrossRef/Sox10nuclear_cluster2.counts.txt')
sox10_3 <- fread('ATACCrossRef/Sox10nuclear_cluster3.counts.txt')
sox10_4 <- fread('ATACCrossRef/Sox10nuclear_cluster4.counts.txt')
sox10_5 <- fread('ATACCrossRef/Sox10nuclear_cluster5.counts.txt')
sox10_6 <- fread('ATACCrossRef/Sox10nuclear_cluster6.counts.txt')
sox10_7 <- fread('ATACCrossRef/Sox10nuclear_cluster7.counts.txt')
sox10_8 <- fread('ATACCrossRef/Sox10nuclear_cluster8.counts.txt')
sox10_9 <- fread('ATACCrossRef/Sox10nuclear_cluster9.counts.txt')
sox10_10 <- fread('ATACCrossRef/Sox10nuclear_cluster10.counts.txt')
cluster1 <- GRanges(sox10_1)
cluster2 <- GRanges(sox10_2)
cluster3 <- GRanges(sox10_3)
cluster4 <- GRanges(sox10_4)
cluster5 <- GRanges(sox10_5)
cluster6 <- GRanges(sox10_6)
cluster7 <- GRanges(sox10_7)
cluster8 <- GRanges(sox10_8)
cluster9 <- GRanges(sox10_9)
cluster10 <- GRanges(sox10_10)

clist <- GRangesList(cluster10,cluster9,cluster8,cluster7,cluster6,cluster5,cluster4,cluster3,cluster2,cluster1)
cnames <- c("cluster10","cluster9","cluster8","cluster7","cluster6","cluster5","cluster4","cluster3","cluster2","cluster1") 
for(i in 1:10){
    # overlap between the bed file and each cluster
    overlaps <- findOverlaps(underGR,clist[[i]],ignore.strand=TRUE)
    if(length(overlaps)>0){
        # change the column name of the samples that overlap, so they're not 'bidirregion#' but "clusteri_#' and we can see how cluster overlaps are distributed in the new k-means 
        colnames(under15)[queryHits(overlaps)] <- paste(cnames[i],"_",1:length(overlaps),sep="")
        
        #save the list of regions that overlap with each cluster, as one regino may overlap with several clusters
        x <- paste(cnames[i],"list",sep="")
        assign(x,underGR[queryHits(overlaps),])
    }
}

# remove all columns with no reads at all
noZ <- under15[, colSums(under15 != 0) > 0]

# set colors for leaves of dendrogram
clusterno <- colnames(noZ)
clusterno <- replace(clusterno,grepl("Bidir",clusterno),11)
clusterno <- gsub('cluster([0-9]+)_.*','\\1',clusterno)
clustercol <- stringr::str_replace_all(clusterno, 
                                        setNames(c("gray","purple","black","gold","orange","red","magenta","green","green4","cyan","cornflowerblue"),
                                                 as.character(c(11,10,1:9))))


#k means
#set new factors for sox10 cluster1-5 and 6-10 as fviz_cluster() cannot plot more than 6 geom_point types simultaneously
lab_col <- colnames(noZ)
lab_col <- gsub('cluster([0-9]+)_.*','\\1',lab_col) 
lab_col <- stringr::str_replace_all(lab_col,setNames(c(rep("3",5),rep("2",5)),as.character(c(10:1))))
lab_col <- replace(lab_col,grepl("Bidir",lab_col),1)
lab_col <- as.numeric(lab_col)

data <- as.data.frame(t(noZ[2:29,]))
scaledata <- scale(data)
data$clusterno <- as.factor(lab_col)
scaledata$clusterno <- as.factor(lab_col)

registerDoMC(15)
ans.kmeans.par <- foreach(i=1:15) %dopar% {
    return(kmeans(scaledata[,1:60],centers=9,nstart=100,iter.max=10,algorithm="MacQueen"))
}

fviz_cluster(ans.kmeans.par, data=scaledata[,1:60], geom = "point", repel=TRUE) + ggtitle("sox10 Reprod Bidir Reg Scaled, k = 9") + geom_point(aes(shape=lab_col))

# getting median and max of each cluster

sumstatsCluster <- function(dt,clusterlist){
    for(i in 1:length(clusterlist)){
        dd <- cbind(as.data.frame(dt), cluster = clusterlist[[i]]$cluster)
        c <- length(clusterlist[[i]]$size)
        for(j in 1:c){
            temp <- dplyr::filter(dd,cluster==j)
            print(paste("Run ",i,", cluster ",j,", sox10 median and max are",sep=""))
            print(median(c(temp$sox10_871_HisatDefault.primary.plus,temp$sox10_873_HisatDefault.primary.plus)))
            print(max(c(temp$sox10_871_HisatDefault.primary.plus,temp$sox10_873_HisatDefault.primary.plus)))
            print("sox10 minus")
            print(median(c(temp$sox10_871_HisatDefault.primary.minus,temp$sox10_873_HisatDefault.primary.minus)))
            print(max(c(temp$sox10_871_HisatDefault.primary.minus,temp$sox10_873_HisatDefault.primary.minus)))
            print("b-actin plus")
            print(median(c(temp$bactin_869_HisatDefault.primary.plus,temp$bactin_870_HisatDefault.primary.plus)))
            print(max(c(temp$bactin_869_HisatDefault.primary.plus,temp$bactin_870_HisatDefault.primary.plus)))
            print("b-actin minus")
            print(median(c(temp$bactin_869_HisatDefault.primary.minus,temp$bactin_870_HisatDefault.primary.minus)))
            print(max(c(temp$bactin_869_HisatDefault.primary.minus,temp$bactin_870_HisatDefault.primary.minus)))
            print("ribo-seq plus")
            print(median(c(temp$sox10_riboseq867_HisatDefault.primary.plus,temp$sox10_riboseq873_HisatDefault.primary.plus)))
            print(max(c(temp$sox10_riboseq867_HisatDefault.primary.plus,temp$sox10_riboseq873_HisatDefault.primary.plus)))
            print("ribo-seq minus")
            print(median(c(temp$sox10_riboseq867_HisatDefault.primary.minus,temp$sox10_riboseq873_HisatDefault.primary.minus)))
            print(max(c(temp$sox10_riboseq867_HisatDefault.primary.minus,temp$sox10_riboseq873_HisatDefault.primary.minus)))
    }}
}

#tSNE
unique <- unique(scaledata)
tsnescaled <- Rtsne(scaledata[,1:60], dims = 2, perplexity=30, verbose=TRUE, max_iter = 100,num_threads=0)


#functions
run.pvclust.iterations <- function(matrix,colorvector,iterations){
    for(i in 1:iterations){
        seed <- runif(1, min = 100, max = 999)
        set.seed(seed)
        print(seed) #maybe make a not-verbose option

        ind <- sample(dim(noZ)[2], 1000)
        samp <- noZ[,ind]
        colors_to_use <- clustercol[ind]
        clust <- pvclust(samp,method.hclust="centroid",method.dist="cor")
        plotd <- as.dendrogram(test)
        colors_to_use <- colors_to_use[order.dendrogram(plotd)]
        dendnew <- assign_values_to_leaves_edgePar(dend=plotd, value = colors_to_use, edgePar = "col") %>% set("labels_cex", 0)  
    
        pngfile <- paste("pvclust_ReprodUnd1500bpRandSamp",seed,".png",sep="")
            plot(dendnew,main="pvclust with leaves clustered by color") 
        dev.off()
    }
}


run.kmeans.iterations <- function(matrix,iterations){
    for(i in 1:iterations){
        seed <- runif(1, min = 100, max = 999)
        set.seed(seed)
        print(seed) 
        ind <- sample(dim(matrix)[2], 1000)
        samp <- noZ[,ind]
        clustercat <- clustercol[,ind]

        kClust <- kmeans(t(samp), centers=9, nstart = 100, iter.max = 10)
        pngfile <- paste("kmeans_ReprodUnd1500bpIter10N10k_",i,".png",sep="")
        png(pngfile,)
            fviz_cluster(kClusters, data=samp, geom="point", repel=TRUE, shape=clustercat) + ggtitle("k = 2") + scale_shape_manual('sox10 cluster designation', values=c(22,23,24))
        dev.off()
    }
}
