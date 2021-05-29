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
options(scipen=999)

gtffile <- "/panfs/qcb-panasas/engie/GRCz11EnhDet/Danio_rerio.GRCz11.99.gtf"
txdb <- makeTxDbFromGFF(gtffile,
                        format="gtf",
                        circ_seqs = character()
                        )
xcripts <- genes(txdb)

filelist <- c("../GRCz11Star/ct711a_150804_hets_nuc1PrimaryReads.bam",
               "../GRCz11Star/ct711a_150804_hets_nuc2PrimaryReads.bam",
               "../BiotaggingTimePts/72hpf_19_10_1/72hpf_191001_PrimaryReads.bam",
               "../BiotaggingTimePts/72hpf_19_10_29/72hpf_191029_PrimaryReads.bam",
               "../BiotaggingTimePts/96hpf_19_8_5/96hpf_190805_PrimaryReads.bam",
               "../BiotaggingTimePts/96hpf_19_9_18/96hpf_190918_PrimaryReads.bam")
bamlist <- BamFileList(filelist)

cond <- c("48","48","72","72","96","96")
colData <- data.frame(time=cond) 
row.names(colData) <- c("ct711a_150804_hets_nuc1PrimaryReads.bam",
                        "ct711a_150804_hets_nuc2PrimaryReads.bam",
                        "72hpf_191001_PrimaryReads.bam",
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

function(eset, class_id, control, treatment)
 {
   #require(DESeq2)
   control_inds <- which(pData(eset)[, class_id] == control)
   treatment_inds <- which(pData(eset)[, class_id] == treatment)
   eset.compare <- eset[, c(control_inds, treatment_inds)]
 
   # make deseq2 compliant dataset
   colData <- data.frame(condition=as.character(pData(eset.compare)[, class_id]))
   dds <- DESeqDataSetFromMatrix(exprs(eset.compare), colData, formula( ~ condition))
 
   # set reference to control, otherwise default is alphabetical order
   dds$condition <- factor(dds$condition, levels=c(control,treatment))
 
   # run deseq2
    3 steps:
      1. estimate size factors
      2. estimate dispersion
      3. negative binomial GLM fitting and wald test
   dds_res <- DESeq(dds)
   res <- results(dds_res)
   res$dispersion <- dispersions(dds_res)
   return(res)
}


