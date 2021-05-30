### PairwiseGOAnalysis.r
###
### Purpose: Gene ontology analyses between conditions: 48hpf, 72hpf, 96hpf
###
###
### Written by Liana Engie
### Last updated: May 2021
###
### runPairGO(bamfile,gtffile)
### Input: 
### (Current) Output: txt file of GO terms, list of top differentially expressed genes

library(GenomicFeatures)
library(GenomicAlignments) 
library(dplyr) 
library(data.table)
library(DESeq2)
library(ggpubr)
library(org.Dr.eg.db)
options(scipen=999)

gtffile <- "/panfs/qcb-panasas/engie/GRCz11EnhDet/Danio_rerio.GRCz11.99.gtf"
txdb <- makeTxDbFromGFF(gtffile,
                        format="gtf",
                        circ_seqs = character()
                        )
xcripts <- genes(txdb)

4872filelist <- c("GRCz11Star/ct711a_150804_hets_nuc1PrimaryReads.bam",
               "GRCz11Star/ct711a_150804_hets_nuc2PrimaryReads.bam",
               "BiotaggingTimePts/72hpf_19_10_1/72hpf_191001_PrimaryReads.bam",
               "BiotaggingTimePts/72hpf_19_10_29/72hpf_191029_PrimaryReads.bam")

4896filelist <- c("GRCz11Star/ct711a_150804_hets_nuc1PrimaryReads.bam",
               "GRCz11Star/ct711a_150804_hets_nuc2PrimaryReads.bam",
               "BiotaggingTimePts/96hpf_19_8_5/96hpf_190805_PrimaryReads.bam",
               "BiotaggingTimePts/96hpf_19_9_18/96hpf_190918_PrimaryReads.bam")

7292filelist <- c("BiotaggingTimePts/72hpf_19_10_1/72hpf_191001_PrimaryReads.bam",
               "BiotaggingTimePts/72hpf_19_10_29/72hpf_191029_PrimaryReads.bam",
               "BiotaggingTimePts/96hpf_19_8_5/96hpf_190805_PrimaryReads.bam",
               "BiotaggingTimePts/96hpf_19_9_18/96hpf_190918_PrimaryReads.bam")

bamlist <- BamFileList(filelist)

#cond <- c("48","48","72","72")
#cond <- c("48","48","96","96")
#cond <- c("72","72","96","96")

colData <- data.frame(time=cond) 
row.names(colData) <- c("ct711a_150804_hets_nuc1PrimaryReads.bam",
                        "ct711a_150804_hets_nuc2PrimaryReads.bam",
                        "72hpf_191001_PrimaryReads.bam",
                        "72hpf_191029_PrimaryReads.bam")
row.names(colData) <- c("ct711a_150804_hets_nuc1PrimaryReads.bam",
                        "ct711a_150804_hets_nuc2PrimaryReads.bam",
                        "96hpf_190805_PrimaryReads.bam",
                        "96hpf_190918_PrimaryReads.bam")
row.names(colData) <- c("72hpf_191001_PrimaryReads.bam",
                        "72hpf_191029_PrimaryReads.bam",
                        "96hpf_190805_PrimaryReads.bam",
                        "96hpf_190918_PrimaryReads.bam")

# Get read counts for your regions
overlaps <- summarizeOverlaps(features=xcripts,reads=bamlist,singleEnd=FALSE,ignore.strand=FALSE) 
genecountsassay <- assay(overlaps)

deseq <- SummarizedExperiment(assays=genecountsassay, rowRanges=xcripts, colData=colData)
dds <- DESeqDataSet(deseq, design= ~time)

keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
ddsdf <- as.data.frame(assay(dds))
rownames(ddsdf) <- rownames(dds)
colnames(ddsdf) <- colnames(dds)
#write.table(ddsdf,"BiotaggingTimePts/tpm4_devtimepts4872_genecounts.csv",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
#write.table(ddsdf,"BiotaggingTimePts/tpm4_devtimepts4896_genecounts.csv",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
#write.table(ddsdf,"BiotaggingTimePts/tpm4_devtimepts7296_genecounts.csv",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)

fdr.threshold <- 0.1
ddsnew <- DESeq(dds)
res <- results(ddsnew, independentFiltering=FALSE) 
# assayed.genes <- rownames(res) check this
de.genes <- rownames(res)[ which(res$padj < fdr.threshold) ] 
write.table(rownames(dds),"BiotaggingTimePts/tpm4_devtimepts4872_genenames.csv",sep="\t",quote=FALSE)
write.table(de.genes,"BiotaggingTimePts/tpm4_devtimepts4872_DEgenenames.csv",sep="\t",quote=FALSE)

write.table(rownames(dds),"BiotaggingTimePts/tpm4_devtimepts4896_genenames.csv",sep="\t",quote=FALSE)
write.table(de.genes,"BiotaggingTimePts/tpm4_devtimepts4896_DEgenenames.csv",sep="\t",quote=FALSE)

write.table(rownames(dds),"BiotaggingTimePts/tpm4_devtimepts7296_genenames.csv",sep="\t",quote=FALSE)
write.table(de.genes,"BiotaggingTimePts/tpm4_devtimepts7296_DEgenenames.csv",sep="\t",quote=FALSE)


assay.genes <- rownames(dds)
gene.vector=as.integer(assay.genes%in%de.genes)
names(gene.vector)=assay.genes
write.table(gene.vector,"BiotaggingTimePts/tpm4_devtimepts4872_goseqinput.csv",sep="\t",quote=FALSE,col.names=TRUE)
write.table(gene.vector,"BiotaggingTimePts/tpm4_devtimepts4896_goseqinput.csv",sep="\t",quote=FALSE,col.names=TRUE)
write.table(gene.vector,"BiotaggingTimePts/tpm4_devtimepts7296_goseqinput.csv",sep="\t",quote=FALSE,col.names=TRUE)

## get the gene lengths for bias analysis
xcriptsKept <- xcripts[genecounts$V1,]
lengthData <- width(xcriptsKept)
medianLengthData <- median(lengthData)
pwf <- goseq::nullp(gene.vector,bias.data=lengthData)

write.table(pwf,"BiotaggingTimePts/tpm4_devtimepts4872_goseqPWF.csv",sep="\t",quote=FALSE,col.names=TRUE)
write.table(pwf,"BiotaggingTimePts/tpm4_devtimepts4896_goseqPWF.csv",sep="\t",quote=FALSE,col.names=TRUE)
write.table(pwf,"BiotaggingTimePts/tpm4_devtimepts7296_goseqPWF.csv",sep="\t",quote=FALSE,col.names=TRUE)

# Using Wallenius approximation
GO.wall <- goseq::goseq(pwf,"danRer7","ensGene")

enriched.GO=GO.wall$category[GO.wall$over_represented_pvalue<.05]
write.table(GO.wall,"BiotaggingTimePts/tpm4_devtimepts4872_goseqGOWall.csv",sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)
write.table(GO.wall,"BiotaggingTimePts/tpm4_devtimepts4896_goseqGOWall.csv",sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)
write.table(GO.wall,"BiotaggingTimePts/tpm4_devtimepts7296_goseqGOWall.csv",sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)

capture.output(for(go in enriched.GO[1:length(enriched.GO)]) { print(GOTERM[[go]])
cat("--------------------------------------\n")
}, file="BiotaggingTimePts/tpm4_devtimepts4872_goseqEnriched.txt")

capture.output(for(go in enriched.GO[1:length(enriched.GO)]) { print(GOTERM[[go]])
cat("--------------------------------------\n")
}, file="BiotaggingTimePts/tpm4_devtimepts4896_goseqEnriched.txt")

capture.output(for(go in enriched.GO[1:length(enriched.GO)]) { print(GOTERM[[go]])
cat("--------------------------------------\n")
}, file="BiotaggingTimePts/tpm4_devtimepts7296_goseqEnriched.txt")



# random sampling to check against the Wallenius approx, if you want it
GO.samp <- goseq::goseq(pwf,"danRer7","ensGene",method="Sampling",repcnt=1000)
plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1],GO.samp[,1]),2]),
    xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
    xlim=c(-3,0))
abline(0,1,col=3,lty=2)

GO.samp <- goseq::goseq(pwf,"danRer7","ensGene",method="Sampling",use_genes_without_cat=TRUE,repcnt=1000)

plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1],GO.samp[,1]),2]),
    xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",main="Random Sampling INCLUDING uncat genes",
    xlim=c(-3,0))
abline(0,1,col=3,lty=2)


