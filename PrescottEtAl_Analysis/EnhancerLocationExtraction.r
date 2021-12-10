###
###
###
###
###


chimpenh <- read_excel("Prescott_S2.xls", sheet = 3)
humanenh <- read_excel("Prescott_S2.xls", sheet = 2)
topenh <- read_excel("Prescott_S2.xls", sheet = 1)

test <- as.data.frame(chimpenh[,c(2,3,4,5)])
colnames(test) <- c("chimp","padj","logfold","basemean")

humandf <- as.data.frame(humanenh[,c(1,3,4,5)])
colnames(humandf) <- c("human","padj","logfold","basemean")

chimp1000 <- test %>% separate(chimp,c("chr","start","end"))
human1000 <- humandf %>% separate(human,c("chr","start","end"))
write.table(chimp1000,"Prescott_Chimp1000TopEnh.bed",quote=FALSE, row.names=FALSE, col.names=FALSE,sep="\t")
write.table(human1000,"Prescott_Human1000TopEnh.bed",quote=FALSE, row.names=FALSE, col.names=FALSE,sep="\t")

chimp1000 <- GRanges(chimp1000)
human1000 <- GRanges(human1000)

loadbd3gr <- function(bd){
    bed <- fread(bd)
    colnames(bed) <- c("chr","start","end")
    return(GRanges(bed))
}

saveOv <- function(bidirfileext,enh){
    bidir <- loadbd3gr(paste(bidirfileext,"_Under10kOverlapsMerged.bed",sep=""))
    ov <- findOverlaps(bidir,enh,ignore.strand=TRUE) #bidir is query, enhancer is subject
    ov_enh <- enh[unique(subjectHits),]
    no_enh <- enh[-unique(subjectHits),]
    print(paste("Number of unique enhancers:",length(ov_enh),sep=" "))
    gr_save(ov_enh,paste(bidirfileext,"_BidirEnhOv",sep=""))
    gr_save(no_enh,paste(bidirfileext,"_BidirNoEnh",sep=""))
    return(ov_enh)
}

h1_1 <- loadbd3gr("SRR2096446.1",human1000)
h1_2 <- loadbd3gr("SRR2096447.1",human1000)
h2_1 <- loadbd3gr("SRR2096448.1",human1000)
h2_2 <- loadbd3gr("SRR2096449.1",human1000)
h3_1 <- loadbd3gr("SRR2096450.1",human1000)
h3_2 <- loadbd3gr("SRR2096451.1",human1000)
c1_1 <- loadbd3gr("SRR2096452.1",chimp1000)
c1_2 <- loadbd3gr("SRR2096453.1",chimp1000)
c2_1 <- loadbd3gr("SRR2096454.1",chimp1000)
c2_2 <- loadbd3gr("SRR2096455.1",chimp1000)
