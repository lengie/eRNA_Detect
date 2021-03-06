### detectEnhancers.R
###
### Purpose: Load RNA-seq data, output locations and FPKMs of regions with bidirectional reads within a given region 
###
###
### Written by Liana Engie
### Last updated: February 2019
###
### detectEnhancers(chromosome, input_start, input_end,strand)
### Input: string chromosome number, int input_start, int input_end, string strand (either "+" or "-")
### Output:

library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)
library(ggplot2)
library(dplyr) 
library(data.table)
library(DESeq2)
library(Rsubread)

detectEnhancers{

	clusters <- findReadRegions(chromosome,input_start,input_end,strand,bed){
		if(strand!="-"||strand!="+"){
			print('Strand should be "-" or "+"')
			break 
		}
		
		interval <- GRanges(seqnames=chromosome,range=IRanges(input_start,input_end),strand=strand)
		list <- subsetByOverlaps(bed,interval) #GRanges obj
		if(length(ranges(list))==0){
			print('No clusters of reads in range.')
		}#clusters <- ranges(list) #IRanges obj
		clusters <- list
		
		#now we have the clusters of reads. Get the read counts
		#Remove read counts that are 0
		#Find read counts of the opposite strand
		#Calculate TPM
		#Report TPM, keep if it's over 2
	}

	range <- reverseStrand(range){
		anti <- strand(range)
		for(i=1:length(runValue(anti))){
			if(runValues(anti)[i]=="-"){
				runValues(anti)[i] <- "+"
			}else if (runValues(anti)[i]=="+"){runValues(anti)[i] <- "-"} #don't want to mess with the *
		}
		strand(range) <- anti
	}

	FPKM <- FPKMCalc(clusters,reads,scaling_factor){
		#Assumes it's given only the region with the suspected enhancer, so it'll check the given regions then the antisense strand
		counts <- summarizeOverlaps(features=clusters,reads=reads,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
		sense_counts <- assay(counts) #a matrix
		#just read strand from GRanges obj
		opp_strand <- reverseStrand(clusters)
		
		counts <- summarizeOverlaps(features=opp_strand,reads=reads,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
		anti_counts <- assay(counts)
		
		
		FPKM <- data.frame(seqnames(clusters),start(clusters),width(clusters),sense_counts,anti_counts)  
		colnames(FPKM) <- c("chr","start","width","sense","antisense")
		FPKM %>% mutate(sense = sense/(scaling_factor*width)) %>% mutate(antisense=antisense/(scaling_factor*width))
	}

	scaling <- TPMScaleFac(reads,merged,totalreads){
		#To calculate the scaling factor for TPM for an experiment, looking at unannotated regions rather than a GTF file's annotated genes
		counts <- summarizeOverlaps(features=merged,reads=reads,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
		allcounts <- assay(counts)
		df <- data.frame(seqnames(clusters),start(clusters),width(clusters),allcounts)  
		colnames(df) <- c("chr","start","width","counts")
		df <- mutate(df, RPK = counts/width)
		scaling <- sum(df$RPK)*totalreads
	}
	
	TPM <- TPMCalc(clusters,reads,scaling){
		#Same as with FPKM, assumes only given region on strand with enhancer
		counts <- summarizeOverlaps(features=clusters,reads=reads,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
		sense_counts <- assay(counts)
		opp_strand <- reverseStrand(clusters)
		
		counts <- summarizeOverlaps(features=opp_strand,reads=reads,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
		anti_counts <- assay(counts)
		
		TPM <- data.frame(seqnames(clusters),start(clusters),width(clusters),sense_counts,anti_counts)  
		colnames(TPM) <- c("chr","start","width","sense","antisense")
		TPM %>% mutate(sense = sense/scaling_factor) %>% mutate(antisense=antisense/scaling_factor)		
	}
	
	overlap <- byChrOverlap(mergedbedchrsplit){
		#input a merged bed file that has been factored by chromosome
		overlap <- data.frame()
		for(i=1:length(mergedbedchrsplit)){
			byStrandTemp <- unlist(mergedbedchrsplit[i])
			strandSplit <- split(byStrandTemp,strand(byStrandTemp))
			compare <- mergeByOverlaps(ranges(strandSplit[[1]]),ranges(strandSplit[[2]]))
			compare <- cbind(rep(as.character(names(mergeSplit)[i]),length(compare)),compare)
			overlap <- rbind(overlap,compare)
		}
		colnames(overlap) <- c('Chr','PlusStrand','MinusStrand')
	}

	### Actual start of the program
	hets1 <- "/auto/cmb-00/rr/engie/RNA/hets1.bam" 
	hets2 <- "/auto/cmb-00/rr/engie/RNA/hets2.bam"
	
	#Number of clusters in RNA-seq data, in millions
	hets1frag <- 50.731011 
	hets2frag <- 56.315336 
	flag <- scanBamFlag(isSecondaryAlignment=FALSE, isDuplicate=FALSE)
	hetsread1 <- readGAlignmentPairs(hets1,param=ScanBamParam(flag=flag))
	hetsread2 <- readGAlignmentPairs(hets2,param=ScanBamParam(flag=flag))
	
	#use merged bed, separate by strand
    #bam file merged with bookends, else '/auto/cmb-00/rr/engie/RNA/merged1_nobookend.bed'
	file <- '/auto/cmb-00/rr/engie/RNA/merged1.bed'
	bed1 <- fread(file,fill=TRUE,verbose=TRUE,data.table=FALSE) 
	bed1R <- GRanges(seqnames=bed1$V1,ranges=IRanges(start=bed1$V2, end=bed1$V3),strand=bed1$V4)
	mergeSplit <- split(bed1R, as.factor(seqnames(bed1R)))
	overlap <- ByChrOverlap(mergeSplit)
	
	#Will use gtf file to ID nc regions later for validation, but calc scaling factor with merge  
	#g <- "/auto/cmb-00/rr/engie/RNA/Danio_rerio.GRCz10.87.chr.gtf" 
	#txdb <- makeTxDbFromGFF(g, format="gtf",circ_seqs = character())
	#sequence <- seqlevels(txdb)
	#new <- mapSeqlevels(sequence,"UCSC")
	#new <- new[complete.cases(new)]
	#txdb <- renameSeqlevels(txdb,new)
	#STILL NEED TO PULL OUT NON CODING REGIONS FROM GTF
	#negative rather than a category itself
	
	#FPKM1 <- FPKMCalc(bed1R,hetsread1,hets1frag) 
	

	
}


	
	
