## Written by Liana Engie, Dr. Scott E Fraser lab
## Last updated 1/25/2021

# Investigating the statistical distribution of nuclear RNA reads
# Uses sox10 nuclear RNA pulled from Trinh et al 2017 Biotagging paper: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89670
# Loads the data, combines both strands into single GRanges object, loads annotations, makes count tables

library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicAlignments)
library(data.table)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(fitdistrplus)

# if using bedgraphs
sox10_nuc1minusfile <- "sox10nuc_minus1.bedGraph"
sox10_nuc1plusfile <- "sox10nuc_plus1.bedGraph"
sox10_nuc2minusfile <- "sox10nuc_minus2.bedGraph"
sox10_nuc2plusfile <- "sox10nuc_plus2.bedGraph"

sox10_nuc1minus <- fread(sox10_nuc1minusfile) 
sox10_nuc1plus <- fread(sox10_nuc1plusfile) 
sox10_nuc2minus <- fread(sox10_nuc2minusfile) 
sox10_nuc2plus <- fread(sox10_nuc2plusfile) 

colnames(sox10_nuc1minus) <- c("seqnames","start","end","score")
colnames(sox10_nuc2minus) <- c("seqnames","start","end","score")
colnames(sox10_nuc2plus) <- c("seqnames","start","end","score")
colnames(sox10_nuc1plus) <- c("seqnames","start","end","score")
sox10_nuc1minus <- mutate(sox10_nuc1minus,strand="-")
sox10_nuc2minus <- mutate(sox10_nuc2minus,strand="-")
sox10_nuc2plus <- mutate(sox10_nuc2plus,strand="+")
sox10_nuc1plus <- mutate(sox10_nuc1plus,strand="+")

# OR if using bigwigs
sox10_nuc1minusfile <- "GSM2386488_Sox10nuclear_minus.sort.bam.bg.bw"
sox10_nuc1plusfile <- "GSM2386488_Sox10nuclear_plus.sort.bam.bg.bw"
sox10_nuc2minusfile <- "GSM2386489_Sox10_nuclear1_minus.bw"
sox10_nuc2plusfile <- "GSM2386489_Sox10_nuclear1_plus.bw"

sox10_nuc1minus <- import(sox10_nuc1minusfile,format="BigWig") 
sox10_nuc1plus <- import(sox10_nuc1plusfile,format="BigWig")
sox10_nuc2minus <- import(sox10_nuc2minusfile,format="BigWig")
sox10_nuc2plus <- import(sox10_nuc2plusfile,format="BigWig")

strand(sox10_nuc1minus) <- "-"
strand(sox10_nuc2minus) <- "-"
strand(sox10_nuc2plus) <- "+"
strand(sox10_nuc1plus) <- "+"

polyAplus <- import("GSM2386493_Sox10_Ribo_plus.sort.bam.bg.bw",format="BigWig")
polyAminus <- import("GSM2386493_Sox10_Ribo_minus.sort.bam.bg.bw",format="BigWig")  

# combining into GRangesLists
grl <- GRangesList(sox10_nuc1minus,sox10_nuc1plus)
grl2 <- GRangesList(sox10_nuc2minus,sox10_nuc2plus)
grA <- GRangesList(polyAplus,polyAminus)

# turn the GRangesLists into GRanges
sox10_nuc1 <- unlist(grl)
sox10_nuc2 <- unlist(grl2)
polyA <- unlist(grA)

# retrieve the annotations to generate count tables
gtffile <- "/panfs/qcb-panasas/engie/GRCz11EnhDet/Danio_rerio.GRCz11.99.gtf"
txdb <- makeTxDbFromGFF(gtffile,
                        format="gtf",
                        circ_seqs = character()
                        )
seqlevelsStyle(txdb) <- "UCSC"
transcripts <- transcripts(txdb)
exons <- exons(txdb)

# trying to use the zebrafish ncRNA database
lncfile <- "/panfs/qcb-panasas/engie/NCReadDist/ZFLNC_lncRNA.gtf"	
lncRNA <- rtracklayer::import(lncfile) # will not load as a TxDb file

# generating read count tables
exoncounts1 <- summarizeOverlaps(features=exons,reads=sox10_nuc1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
exoncounts2 <- summarizeOverlaps(features=exons,reads=sox10_nuc2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
exoncounts1tb <- assay(exoncounts1)
exoncounts2tb <- assay(exoncounts2)

nccounts1 <- summarizeOverlaps(features=lncRNA,reads=sox10_nuc1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
nccounts2 <- summarizeOverlaps(features=lncRNA,reads=sox10_nuc2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE) 
nccounts1tb <- assay(nccounts1)
nccounts2tb <- assay(nccounts2)

txcounts1 <- summarizeOverlaps(features=transcripts,reads=sox10_nuc1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
txcounts2 <- summarizeOverlaps(features=transcripts,reads=sox10_nuc2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE) 
txcounts1tb <- assay(txcounts1)
txcounts2tb <- assay(txcounts2)

exonA <- summarizeOverlaps(features=exons,reads=polyA,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
exonAtb <- assay(exonA)

ncA <- summarizeOverlaps(features=lncRNA,reads=polyA,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
ncAtb <- assay(ncA)

txA <- summarizeOverlaps(features=transcripts,reads=polyA,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
txAtb <- assay(txA)

## DESeq Analysis for each type of data
cond <- c("nuc","nuc","polyA")
colData <- data.frame(cond=cond) 
row.names(colData) <- c("sox10nuc1","sox10nuc2","sox10polyA")
colData$cond <- factor(colData$cond)

# exons
exoncounts <- data.frame(one = exoncounts1tb,
                         two = exoncounts2tb,
                         poly = exonAtb)

dds <- DESeqDataSetFromMatrix(countData = exoncounts,
                              colData = colData,
                              design = ~ cond)
dds <- DESeq(dds)
norm <- assays(dds)
normtb <- norm[[1]]

#transcripts
txcounts <- data.frame(one = txcounts1tb,
                       two = txcounts2tb,
                       poly = txAtb)

ddstx <- DESeqDataSetFromMatrix(countData = txcounts,
                              colData = colData,
                              design = ~ cond)
ddstx <- DESeq(ddstx)
normtx <- assays(ddstx)
normtxtb <- normtx[[1]]

#lncRNAs
lnccounts <- data.frame(one = nccounts1tb,
                        two = nccounts2tb,
                        poly = ncAtb)

ddsnc <- DESeqDataSetFromMatrix(countData = lnccounts,
                              colData = colData,
                              design = ~ cond)
ddsnc <- DESeq(ddsnc)
normnc <- assays(ddsnc)
normnctb <- normnc[[1]]

exonNorm <- as.data.frame(normtb)
txNorm <- as.data.frame(normtxtb)
lncNorm <- as.data.frame(normnctb)

#comparing statistical distributions
exon.norm.params <- fitdistr(exonNorm$V1,"normal")$estimate
exon.poisson.params <- fitdistr(exonNorm$V1,"poisson")$estimate
exon.negbinom.params <- fitdistr(exonNorm$V1,"negative binomial", lower=c(0,0), method = "SANN")$estimate
exon.dist.params <- map(list(Normal = exon.norm.params,Poisson = exon.poisson.params,`Negative Binomial` = exon.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

exon2.norm.params <- fitdistr(exonNorm$V2,"normal")$estimate
exon2.poisson.params <- fitdistr(exonNorm$V2,"poisson")$estimate
exon2.negbinom.params <- fitdistr(exonNorm$V2,"negative binomial", lower=c(0,0), method = "SANN")$estimate
exon2.dist.params <- map(list(Normal = exon2.norm.params,Poisson = exon2.poisson.params,`Negative Binomial` = exon2.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

exonA.norm.params <- fitdistr(exonNorm$V3,"normal")$estimate
exonA.poisson.params <- fitdistr(exonNorm$V3,"poisson")$estimate
exonA.negbinom.params <- fitdistr(exonNorm$V3,"negative binomial", lower=c(0,0), method = "SANN")$estimate
exonA.dist.params <- map(list(Normal = exonA.norm.params,Poisson = exonA.poisson.params,`Negative Binomial` = exonA.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

tx.norm.params <- fitdistr(txNorm$V1,"normal")$estimate
tx.poisson.params <- fitdistr(txNorm$V1,"poisson")$estimate
tx.negbinom.params <- fitdistr(txNorm$V1,"negative binomial", lower=c(0,0), method = "SANN")$estimate
tx.dist.params <- map(list(Normal = tx.norm.params,Poisson = tx.poisson.params,`Negative Binomial` = tx.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

tx2.norm.params <- fitdistr(txNorm$V2,"normal")$estimate
tx2.poisson.params <- fitdistr(txNorm$V2,"poisson")$estimate
tx2.negbinom.params <- fitdistr(txNorm$V2,"negative binomial", lower=c(0,0), method = "SANN")$estimate
tx2.dist.params <- map(list(Normal = tx2.norm.params,Poisson = tx2.poisson.params,`Negative Binomial` = tx2.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

txA.norm.params <- fitdistr(txNorm$V3,"normal")$estimate
txA.poisson.params <- fitdistr(txNorm$V3,"poisson")$estimate
txA.negbinom.params <- fitdistr(txNorm$V3,"negative binomial", lower=c(0,0), method = "SANN")$estimate
txA.dist.params <- map(list(Normal = txA.norm.params,Poisson = txA.poisson.params,`Negative Binomial` = txA.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

nc.norm.params <- fitdistr(normnctb$V1,"normal")$estimate
nc.poisson.params <- fitdistr(normnctb$V1,"poisson")$estimate
nc.negbinom.params <- fitdistr(normnctb$V1,"negative binomial", lower=c(0,0), method = "SANN")$estimate
nc.dist.params <- map(list(Normal = nc.norm.params,Poisson = nc.poisson.params,`Negative Binomial` = nc.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

nc2.norm.params <- fitdistr(normnctb$V2,"normal")$estimate
nc2.poisson.params <- fitdistr(normnctb$V2,"poisson")$estimate
nc2.negbinom.params <- fitdistr(normnctb$V2,"negative binomial", lower=c(0,0), method = "SANN")$estimate
nc2.dist.params <- map(list(Normal = nc2.norm.params,Poisson = nc2.poisson.params,`Negative Binomial` = nc2.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

ncA.norm.params <- fitdistr(normnctb$V3,"normal")$estimate
ncA.poisson.params <- fitdistr(normnctb$V3,"poisson")$estimate
ncA.negbinom.params <- fitdistr(normnctb$V3,"negative binomial", lower=c(0,0), method = "SANN")$estimate
ncA.dist.params <- map(list(Normal = ncA.norm.params,Poisson = ncA.poisson.params,`Negative Binomial` = ncA.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

# plots
mybinwidth <- 1

exon1 <- ggplot(data=exonNorm, aes(x=sox10nuc1)) + geom_histogram(binwidth=1) + xlim(0,100)+labs(x="sox10 nuc1 exon read counts") + ylim(0,30000)
exon2 <- ggplot(data=exonNorm, aes(x=sox10nuc2)) + geom_histogram(binwidth=1) + xlim(0,100)+labs(x="sox10 nuc2 exon read counts") + ylim(0,30000)
exonA <- ggplot(data=exonNorm, aes(x=sox10polyA)) + geom_histogram(binwidth=1) + xlim(0,100)+labs(x="sox10 polyA exon read counts") + ylim(0,30000)

tx1 <- ggplot(data=txNorm, aes(x=sox10nuc1)) + geom_histogram(binwidth=1) + xlim(0,100)+labs(x="sox10 nuc1 transcript read counts") + ylim(0,1500)
tx2 <- ggplot(data=txNorm, aes(x=sox10nuc2)) + geom_histogram(binwidth=1) + xlim(0,100)+labs(x="sox10 nuc2 transcript read counts") + ylim(0,1500)
txA <- ggplot(data=txNorm, aes(x=sox10polyA)) + geom_histogram(binwidth=1) + xlim(0,100)+labs(x="sox10 polyA  read counts") + ylim(0,1500)

lnc1 <- ggplot(data=lncNorm, aes(x=sox10nuc1)) + geom_histogram(binwidth=1) + xlim(0,100)+labs(x="sox10 nuc1 exon read counts") + ylim(0,9000)+
            stat_function(aes(color = "black"),
                          fun=function(x,mean,sd) mybinwidth * nrow(lncNorm) * dnorm(x,mean, sd),
                          args=fitdistr(lncNorm$sox10polyA,"normal")$estimate) +
            stat_function(aes(color = "blue"),
                          fun=function(x,lambda) mybinwidth * nrow(lncNorm) * dpois(x,lambda),
                          args=fitdistr(lncNorm$sox10polyA,"poisson")$estimate,
                          xlim=c(1,100), n=100) + 
            stat_function(aes(color = "orange"),
                          fun=function(x,size, mu) mybinwidth * nrow(lncNorm) * dnbinom(x,size = size, mu = mu),
                          args=fitdistr(lncNorm$sox10polyA,"negative binomial", method="SANN")$estimate,
                          lower=c(0,0),n=100) + 
            scale_color_manual("Distribution", 
                                values=c(black="black",blue="blue",orange="orange"), labels=ncA.dist.params)
lnc2 <- ggplot(data=lncNorm, aes(x=sox10nuc2)) + geom_histogram(binwidth=1) + xlim(0,100)+labs(x="sox10 nuc2 exon read counts") + ylim(0,9000)+
            stat_function(aes(color = "black"),
                          fun=function(x,mean,sd) mybinwidth * nrow(lncNorm) * dnorm(x,mean, sd),
                          args=fitdistr(lncNorm$sox10polyA,"normal")$estimate) +
            stat_function(aes(color = "blue"),
                          fun=function(x,lambda) mybinwidth * nrow(lncNorm) * dpois(x,lambda),
                          args=fitdistr(lncNorm$sox10polyA,"poisson")$estimate,
                          xlim=c(1,100), n=100) + 
            stat_function(aes(color = "orange"),
                          fun=function(x,size, mu) mybinwidth * nrow(lncNorm) * dnbinom(x,size = size, mu = mu),
                          args=fitdistr(lncNorm$sox10polyA,"negative binomial", method="SANN")$estimate,
                          lower=c(0,0),n=100) + 
            scale_color_manual("Distribution", 
                                values=c(black="black",blue="blue",orange="orange"), labels=ncA.dist.params)
lncA <- ggplot(data=lncNorm, aes(x=sox10polyA)) + geom_histogram(binwidth=mybinwidth) + xlim(0,100)+labs(x="sox10 polyA lncRNA read counts") + ylim(0,4000) +
            stat_function(aes(color = "black"),
                          fun=function(x,mean,sd) mybinwidth * nrow(lncNorm) * dnorm(x,mean, sd),
                          args=fitdistr(lncNorm$sox10polyA,"normal")$estimate) +
            stat_function(aes(color = "blue"),
                          fun=function(x,lambda) mybinwidth * nrow(lncNorm) * dpois(x,lambda),
                          args=fitdistr(lncNorm$sox10polyA,"poisson")$estimate,
                          xlim=c(1,100), n=100) + 
            stat_function(aes(color = "orange"),
                          fun=function(x,size, mu) mybinwidth * nrow(lncNorm) * dnbinom(x,size = size, mu = mu),
                          args=fitdistr(lncNorm$sox10polyA,"negative binomial", method="SANN")$estimate,
                          lower=c(0,0),n=100) + 
            scale_color_manual("Distribution", 
                                values=c(black="black",blue="blue",orange="orange"), labels=ncA.dist.params)

# perform chi-squared tests for goodness of fit of the distributions
length <- nrow(exonNorm)-1
fit <- dnbinom(0:length,size=0.2569002,mu=8.4660221)
chisq.test(exonNorm$V1, p=fit)
fit <- dnbinom(0:length,size=0.1616598,mu=6.1709287)
chisq.test(exonNorm$V2, p=fit)
fit <- dnbinom(0:length,size=0.0653513,mu=3.1287361)
chisq.test(exonNorm$V3, p=fit)

length <- nrow(txNorm)-1
fit <- dnbinom(0:length,size=0.2568344,mu=8.4563691)
chisq.test(txNorm$V1, p=fit)
fit <- dnbinom(0:length,size=0.1617959,mu=6.1837367)
chisq.test(txNorm$V2, p=fit)
fit <- dnbinom(0:length,size=0.06555382,mu=3.14159576)
chisq.test(txNorm$V3, p=fit)

length <- nrow(lncNorm)-1
fit <- dnbinom(0:length,size=0.1921967,mu=47.5004140)
chisq.test(normnctb$V1, p=fit)
chisq.test(normnctb$V2, p=fit)
fit <- dnbinom(0:length, size=0.06592968, mu=42.50703446)
chisq.test(normnctb$V3, p=fit)

