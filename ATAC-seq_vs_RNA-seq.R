### After sox10 ATAC processing and sox10 nuc-seq processing, comparing them to pull out bidirectional regions that overlap with ATAC seq using only RNA-seq data
###
###
### Last Updated: July 2021

## low count filtering RNA-seq
exonfile <- fread("Zv9exons.bed")
colnames(exonfile) <- c("chr","start","end","ID","score","strand")
exonfile <- mutate(exonfile,width=end-start)
reprod <- fread("sox10_Zv9FromFASTQNoncodingUnder10kReprodOnly.bed")
colnames(reprod) <- c("chr","start","end")
reprod <- mutate(reprod,width=end-start)


# just checking that in the process of combining reproducible bidirectional regions, I didn't add coding regions back in
remcoding <- function(gr){
    hits <- GenomicRanges::setdiff(gr,coding)
    hits <- hits[!width(hits)==1]
    return(hits)
}

coding <- fread("Zv9ExonsUTRs500bpFlanking.bed")
colnames(coding) <- c("chr","start","end","width","strand","ID","score")
coding <- GRanges(coding)
seqlevelsStyle(coding) <- "ensembl"
newr <- remcoding(GRanges(reprod))

# create a data frame of features from both exons, bidirectional regions, in both strands
feat <- data.frame(chr=c(exonfile$chr,exonfile$chr,newdf$seqnames,newdf$seqnames),
                   start=c(exonfile$start,exonfile$start,newdf$start,newdf$start),
                   end=c(exonfile$end,exonfile$end,newdf$end,newdf$end),
                   strand=c(rep(c("+","-"),each=nrow(exonfile)),rep(c("+","-"),each=nrow(newdf))),
                   width=c(exonfile$width,exonfile$width,newdf$width,newdf$width)
)
featgr <- GRanges(feat)
# the bam files are in Ensembl style not UCSC
seqlevelsStyle(featgr) <- "ensembl"

# get the read counts for both strands, exons and bidirectional regions
counts <- summarizeOverlaps(features=featgr,reads=bamlist,ignore.strand=FALSE)      
genecounts <- assay(counts)

# set up the read counts per base pair data for plotting
toplot <- data.frame(nuc1=genecounts[,1]/feat$width,
                    nuc2=genecounts[,2]/feat$width,
                    strand=feat$strand)

# create a biological replicate plot
png("biorepplot_exonsbidir.png",width=1200,height=1080)
    ggplot(toplot, aes(x=nuc1, y=nuc2,color=strand)) +
      geom_point(size=0.25) + xlab("Nuc1") + ylab("Nuc2") +
      ggtitle("sox10 biological replicates, exons' and bidir regions' read counts by length") +
      scale_color_manual(name="Strand",
                         labels=c("Plus","Minus"),
                         values=c("tan3","cadetblue4"))
dev.off() 

# low read count filtering
filter <- toplot[which(toplot[,1]==0),]
max(filter$nuc2)
filter <- toplot[which(toplot[,2]==0),]
max(filter$nuc1)
# chose the smaller. Turn this into a function later



## ATAC-seq overlapping with bidirectional regions: what regions are we targetting as putative enhancers?
comp <- function(gr1,gr2){
    overlaps <- findOverlaps(gr1,gr2,ignore.strand=TRUE)
    hits <- gr1[unique(queryHits(overlaps)),]
    nohits <- gr1[-queryHits(overlaps),]
    print(length(hits))
    return(GRangesList(hits,nohits))
}


genrich <- fread("sox10ATAC_genrich.narrowPeaks")
dualMACs <- fread("sox10ATAC_MACS2_narrowPeaks.bed")

colnames(dualMACs) <- c("chr","start","end","name","score","strand","fold_enrichment","pval","qval","blockCount")
colnames(genrich) <- c("chr","start","end","name","score","strand","signalValue","pval","qval","peak")


genp <- dplyr::filter(genrich,pval<=5)
macsp <- dplyr::filter(dualMACs,pval<=5) 

reprod <- fread("sox10_Zv9FromFASTQNoncodingUnder10kReprodOnly.bed")
colnames(reprod) <- c("chr","start","end")

greprod <- GRanges(reprod)
ggen <- GRanges(genrich)
gmacs <- GRanges(dualMACs)

genoverlaps <- comp(greprod,ggen)
macsov <- comp(greprod,gmacs)

# used createMatrix in deeptools to bin transcription across the span of the region
genovfile <- fread("sox10_BidirOverlapGenrich3kb50bpbins.tabular")
# the way these were generated, it's nuc1plus, nuc2minus, nuc2plus, nuc2minus

#macsov <- as.matrix(macsovfile[,61:180])
genov <- as.matrix(genovfile[,61:180])
bins <- c(paste(seq(-30,-1,by=1),"minus",sep=""),
            paste(seq(1,30,by=1),"minus",sep=""),
            paste(seq(-30,-1,by=1),"plus",sep=""),
            paste(seq(1,30,by=1),"plus",sep=""))
#colnames(macsov) <- bins
colnames(genov) <- bins

pv_g <- pvclust(t(genov),method.hclust="centroid",method.dist="cor", parallel=TRUE)

# bidirectional regions not overlapping with ATAC
genNofile <- fread("sox10_BidirNotOverlapGenrich3kb50bpbins.tabular")
genNo <- as.matrix(genNofile[,61:180])
bins <- c(paste(seq(-30,-1,by=1),"minus",sep=""),
            paste(seq(1,30,by=1),"minus",sep=""),
            paste(seq(-30,-1,by=1),"plus",sep=""),
            paste(seq(1,30,by=1),"plus",sep=""))
colnames(genNo) <- bins
rownames(genNo) <- paste("ATACOv",1:nrow(genNo),sep="")

# trying this with 20 64GB nodes
tic()
pv_g <- pvclust(t(genNo),method.hclust="centroid",method.dist="cor", parallel=TRUE) 
toc()

cluster <- cutree(pv_g,k=3)

png("testPCA.png",width=1920,height=1080)
    fviz_dend(as.dendrogram(pv_g), k = 3, # Cut in four groups
              cex = 0.5, # label size
              k_colors = c("#00AFBB", "#E7B800", "#FC4E07"),
              color_labels_by_k = TRUE #, # color labels by groups
              #rect = TRUE # Add rectangle around groups
    )
dev.off()
