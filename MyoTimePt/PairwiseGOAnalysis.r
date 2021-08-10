### PairwiseGOAnalysis.r
###
### Purpose: Gene ontology analyses between conditions: 48hpf, 72hpf, 96hpf
### Used package GoSeq after finding differential expression with DESeq2. 
### Corrected for gene length bias and a Wallenius approximation to calculate p-values when fetching GO annotations. Enriched GO terms had pval>.5
###
###
### Written by Liana Engie
### Last updated: June 2021
###
### runPairGO(bamfilelist,features)
### Input: 
### Output: txt file of GO terms, list of top differentially expressed genes. Can be used in PlotGOAnnotations.r

library(GenomicFeatures)
library(GenomicAlignments) 
library(dplyr) 
library(data.table)
library(DESeq2)
library(ggpubr)
library(goseq)
library(GO.db)
library(org.Dr.eg.db)
options(scipen=999)

# load the gene list from a GTF file
gtffile <- "/panfs/qcb-panasas/engie/GRCz11EnhDet/Danio_rerio.GRCz11.99.gtf"
txdb <- makeTxDbFromGFF(gtffile,
                        format="gtf",
                        circ_seqs = character()
                        )
features <- genes(txdb)

# load the processed BAM files (not normalized, as DESeq requires)
filelist <- c("GRCz11Star/ct711a_150804_hets_nuc1PrimaryReads.bam",
               "GRCz11Star/ct711a_150804_hets_nuc2PrimaryReads.bam",
               "BiotaggingTimePts/72hpf_19_10_1/72hpf_191001_PrimaryReads.bam",
               "BiotaggingTimePts/72hpf_19_10_29/72hpf_191029_PrimaryReads.bam")
label <- c("4872")

filelist <- c("GRCz11Star/ct711a_150804_hets_nuc1PrimaryReads.bam",
              "GRCz11Star/ct711a_150804_hets_nuc2PrimaryReads.bam",
               "BiotaggingTimePts/96hpf_19_8_5/96hpf_190805_PrimaryReads.bam",
               "BiotaggingTimePts/96hpf_19_9_18/96hpf_190918_PrimaryReads.bam")
label <- c("4896")

filelist <- c("BiotaggingTimePts/72hpf_19_10_1/72hpf_191001_PrimaryReads.bam",
               "BiotaggingTimePts/72hpf_19_10_29/72hpf_191029_PrimaryReads.bam",
               "BiotaggingTimePts/96hpf_19_8_5/96hpf_190805_PrimaryReads.bam",
               "BiotaggingTimePts/96hpf_19_9_18/96hpf_190918_PrimaryReads.bam")
label <- c("7296")

# Create the experimental condition data frame
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

PairwiseRunDESeq <- function(filelist,features,label,colData){
  bamlist <- BamFileList(filelist)

  # Get read counts for the genes for all conditions
  overlaps <- summarizeOverlaps(features=features,reads=bamlist,singleEnd=FALSE,ignore.strand=FALSE) 
  genecountsassay <- assay(overlaps)

  deseq <- SummarizedExperiment(assays=genecountsassay, rowRanges=xcripts, colData=colData)
  dds <- DESeqDataSet(deseq, design= ~time)

  # Remove genes with no reads in all conditions
  keep <- rowSums(counts(dds)) > 1
  dds <- dds[keep,]
  ddsdf <- as.data.frame(assay(dds))
  rownames(ddsdf) <- rownames(dds)
  colnames(ddsdf) <- colnames(dds)

  # Differentially expressed genes
  fdr.threshold <- 0.1
  ddsnew <- DESeq(dds)
  res <- results(ddsnew, independentFiltering=FALSE) 

  de.genes <- rownames(res)[ which(res$padj < fdr.threshold) ] 
  filen <- paste("tpm4_devtimepts",label,"_genenames.csv",sep="")
  write.table(rownames(dds),filen,sep="\t",quote=FALSE)
  filen2 <- paste("tpm4_devtimepts",label,"_DEgenenames.csv",sep="")

  # Create a binary matrix for the genes expressed in each condition
  assay.genes <- rownames(dds)
  gene.vector=as.integer(assay.genes%in%de.genes)
  names(gene.vector)=assay.genes
  filen3 <- paste("tpm4_devtimepts",label,"_goseqinput.csv",sep="")
  write.table(gene.vector,filen3,sep="\t",quote=FALSE,col.names=TRUE)

  return(res)
}

DEupdown <- function(df,label,log2FC,fdr.threshold){
    DE <- df[ df$padj < fdr.threshold, ]
    up <- subset(DE, log2FoldChange > log2FC)  
    down <- subset(DE, log2FoldChange < -log2FC)

    write.table(up,paste("tpm4_devtimepts_DE",label,"up.csv",sep=""),sep="\t",quote=FALSE)
    write.table(down,paste("tpm4_devtimepts_DE",label,"down.csv",sep=""),sep="\t",quote=FALSE)
    return(list(up,down))
    assay.genes <- names(res)
    xcriptsKept <- features[names(features) %in% assay.genes,]
    lengthData <- width(xcriptsKept)
    
  for(i in 1:2){
        de.genes <- rownames(list(up,down)[[i]])
        gene.vector=as.integer(as.matrix(assay.genes)%in%de.genes)
        names(gene.vector)=as.matrix(assay.genes)
        pwfn <- paste("pwf",l[i],sep="")
        assign(pwfn,goseq::nullp(gene.vector, "danRer11", "ensGene", bias.data=lengthData))   
        l <- c("up","down")
        labelupdown <- paste(label,l[i],sep="")
        filepwf <- paste("tpm4_devtimepts",labelupdown,"_goseqPWF.csv",sep="")
        write.table(pwf,filepwf,sep="\t",quote=FALSE,col.names=TRUE)
    }
    plist <- list(pwfup,pwfdown)
}

for(i in 1:2){
  # Using Wallenius approximation
  GO.wall <- goseq::goseq(plist[[i]],"danRer11","ensGene")

  enriched.GO=GO.wall$category[GO.wall$over_represented_pvalue<.05]
  filegw <- paste("Pairwise/tpm4_devtimepts",label,"_goseqGOWall.csv",sep="")
  write.table(GO.wall,filegw,sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)
  labelupdown <- paste(label,l[i],sep="")
  filen <- paste("Pairwise/tpm4_devtimepts",labelupdown,"_goseqEnriched.txt",sep="")
  capture.output(for(go in enriched.GO[1:length(enriched.GO)]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
  }
  , file=filen)
}

  # random sampling to check against the Wallenius approx, if you want it
  GO.samp <- goseq::goseq(pwf,"danRer11","ensGene",method="Sampling",repcnt=1000)
  plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1],GO.samp[,1]),2]),
      xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
      xlim=c(-3,0))
  abline(0,1,col=3,lty=2)

  GO.samp <- goseq::goseq(pwf,"danRer11","ensGene",method="Sampling",use_genes_without_cat=TRUE,repcnt=1000)

  plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1],GO.samp[,1]),2]),
      xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",main="Random Sampling INCLUDING uncat genes",
      xlim=c(-3,0))
  abline(0,1,col=3,lty=2)



  ## get the gene lengths for bias analysis
  xcriptsKept <- features[names(features) %in% assay.genes, ]
  
  lengthData <- width(xcriptsKept)
  medianLengthData <- median(lengthData)
  pwf <- goseq::nullp(gene.vector,bias.data=lengthData)

  filepwf <- paste("tpm4_devtimepts",label,"_goseqPWF.csv",sep="")
  write.table(pwf,filepwf,sep="\t",quote=FALSE,col.names=TRUE)

  # Using Wallenius approximation
  GO.wall <- goseq::goseq(pwf,"danRer11","ensGene")

  enriched.GO=GO.wall$category[GO.wall$over_represented_pvalue<.05]
  filegw <- paste("tpm4_devtimepts",label,"_goseqGOWall.csv",sep="")
  write.table(GO.wall,filegw,sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)

  filen <- paste("tpm4_devtimepts",label,"_goseqEnriched.txt",sep="")
  capture.output(for(go in enriched.GO[1:length(enriched.GO)]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
  }
  , file=filen)

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
}

### separating up and down regulated genes
fdr.threshold <- 0.1
FC.threshold <- 2

upanddownfile <- function(genenames){
    genecounts <- fread(genenames)
    xcriptsKept <- features[genecounts$V1,]
    lengthData <- width(xcriptsKept)

    pwf <- goseq::nullp(gene.vector,bias.data=lengthData) 
    pwf.down  = nullp(gene.vector, "danRer11", "ensGene", bias.data = lengthData.down)
    pwf.up  = nullp(gene.vector, "danRer11", "ensGene", bias.data = lengthData.up) 
}

GO.wall.up <- goseq(pwf.up,"danRer11","ensGene")
GO.wall.up$padj <- p.adjust(GO.wall.up$over_represented_pvalue, method="BH")
GO.wall.up
GO.wall.up <- subset(GO.wall.up, GO.wall.up$padj<.05)

GO.wall.down <- goseq(pwf.down, "hg38","ensGene")
GO.wall.down$padj <- p.adjust(GO.wall.down$over_represented_pvalue, method="BH")
GO.wall.down <- subset(GO.wall.down, GO.wall.down$padj<.05) }

 ###example

  bamlist <- BamFileList(filelist)

  # Get read counts for the genes for all conditions
  overlaps <- summarizeOverlaps(features=features,reads=bamlist,singleEnd=FALSE,ignore.strand=FALSE) 
  genecountsassay <- assay(overlaps)

  deseq <- SummarizedExperiment(assays=genecountsassay, rowRanges=xcripts, colData=colData)
  dds <- DESeqDataSet(deseq, design= ~time)

  # Remove genes with no reads in all conditions
  keep <- rowSums(counts(dds)) > 1
  dds <- dds[keep,]
  ddsdf <- as.data.frame(assay(dds))
  rownames(ddsdf) <- rownames(dds)
  colnames(ddsdf) <- colnames(dds)
  #write.table(ddsdf,"BiotaggingTimePts/tpm4_devtimepts4872_genecounts.csv",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
  #write.table(ddsdf,"BiotaggingTimePts/tpm4_devtimepts4896_genecounts.csv",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
  #write.table(ddsdf,"BiotaggingTimePts/tpm4_devtimepts7296_genecounts.csv",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)

  # Differentially expressed genes
  fdr.threshold <- 0.1
  ddsnew <- DESeq(dds)
  res <- results(ddsnew, independentFiltering=FALSE) 

  de.genes <- rownames(res)[ which(res$padj < fdr.threshold) ] 
  write.table(rownames(dds),"BiotaggingTimePts/tpm4_devtimepts4872_genenames.csv",sep="\t",quote=FALSE)
  write.table(de.genes,"BiotaggingTimePts/tpm4_devtimepts4872_DEgenenames.csv",sep="\t",quote=FALSE)

  write.table(rownames(dds),"BiotaggingTimePts/tpm4_devtimepts4896_genenames.csv",sep="\t",quote=FALSE)
  write.table(de.genes,"BiotaggingTimePts/tpm4_devtimepts4896_DEgenenames.csv",sep="\t",quote=FALSE)

  write.table(rownames(dds),"BiotaggingTimePts/tpm4_devtimepts7296_genenames.csv",sep="\t",quote=FALSE)
  write.table(de.genes,"BiotaggingTimePts/tpm4_devtimepts7296_DEgenenames.csv",sep="\t",quote=FALSE)

  # Create a binary matrix for the genes expressed in each condition
  assay.genes <- rownames(dds)
  gene.vector=as.integer(assay.genes%in%de.genes)
  names(gene.vector)=assay.genes
  write.table(gene.vector,"BiotaggingTimePts/tpm4_devtimepts4872_goseqinput.csv",sep="\t",quote=FALSE,col.names=TRUE)
  write.table(gene.vector,"BiotaggingTimePts/tpm4_devtimepts4896_goseqinput.csv",sep="\t",quote=FALSE,col.names=TRUE)
  write.table(gene.vector,"BiotaggingTimePts/tpm4_devtimepts7296_goseqinput.csv",sep="\t",quote=FALSE,col.names=TRUE)

  ## get the gene lengths for bias analysis
  xcriptsKept4872 <- xcripts[1:length(dds),]
  xcriptsKept4896 <- xcripts[1:length(dds),]
  xcriptsKept7296 <- xcripts[1:length(dds),]

  lengthData <- width(xcriptsKept4872)
  lengthData <- width(xcriptsKept4896)
  lengthData <- width(xcriptsKept7296)

  medianLengthData <- median(lengthData)
  pwf <- goseq::nullp(gene.vector,bias.data=lengthData)
  #plot+title("48hpf vs 72hpf")
  #plot+title("48hpf vs 96hpf")
  #plot+title("72hpf vs 96hpf")

  write.table(pwf,"BiotaggingTimePts/tpm4_devtimepts4872_goseqPWF.csv",sep="\t",quote=FALSE,col.names=TRUE)
  write.table(pwf,"BiotaggingTimePts/tpm4_devtimepts4896_goseqPWF.csv",sep="\t",quote=FALSE,col.names=TRUE)
  write.table(pwf,"BiotaggingTimePts/tpm4_devtimepts7296_goseqPWF.csv",sep="\t",quote=FALSE,col.names=TRUE)

  # Using Wallenius approximation
  GO.wall <- goseq::goseq(pwf,"danRer11","ensGene")

  enriched.GO=GO.wall$category[GO.wall$over_represented_pvalue<.05]
  write.table(GO.wall,"BiotaggingTimePts/tpm4_devtimepts4872_goseqGOWall.csv",sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)
  write.table(GO.wall,"BiotaggingTimePts/tpm4_devtimepts4896_goseqGOWall.csv",sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)
  write.table(GO.wall,"BiotaggingTimePts/tpm4_devtimepts7296_goseqGOWall.csv",sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)

  capture.output(for(go in enriched.GO[1:length(enriched.GO)]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
  }
  , file="BiotaggingTimePts/tpm4_devtimepts4872_goseqEnriched.txt")

  capture.output(for(go in enriched.GO[1:length(enriched.GO)]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
  }, file="BiotaggingTimePts/tpm4_devtimepts4896_goseqEnriched.txt")

  capture.output(for(go in enriched.GO[1:length(enriched.GO)]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
  }, file="BiotaggingTimePts/tpm4_devtimepts7296_goseqEnriched.txt")


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
