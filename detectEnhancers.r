### detectEnhancers.R
###
### Purpose: Take RNA-seq data and output locations and FPKMs of regions with bidirectional reads within a given region 
###
###
### Written by Liana Engie
### Last updated: April 2017
###
### detectEnhancers(chromosome, input_start, input_end,strand)
### Input: string chromosome number, 
### Output:
library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)
library(ggplot2)
library(dtplyr) 
library(data.table)

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

	hets1 <- "/auto/cmb-00/rr/engie/RNA/Aligned.sortedByCoord.out.bam" 
	hets2 <- "/auto/cmb-00/rr/engie/RNA/hets2.bam"
	hets1frag <- 50.731011 
	hets2frag <- 56.315336 
	flag <- scanBamFlag(isSecondaryAlignment=FALSE, isDuplicate=FALSE)
	hetsread1 <- readGAlignmentPairs(hets1,param=ScanBamParam(flag=flag))
	hetsread2 <- readGAlignmentPairs(hets2,param=ScanBamParam(flag=flag))
	#Will use gtf file later for validation, but calc scaling factor with merge  
	#g <- "/auto/cmb-00/rr/engie/RNA/Danio_rerio.GRCz10.87.chr.gtf" 
	#txdb <- makeTxDbFromGFF(g, format="gtf",circ_seqs = character())
	#sequence <- seqlevels(txdb)
	#new <- mapSeqlevels(sequence,"UCSC")
	#new <- new[complete.cases(new)]
	#txdb <- renameSeqlevels(txdb,new)
	file <- '/auto/cmb-00/rr/engie/RNA/merged1.bed'
	bed1 <- fread(file,fill=TRUE,verbose=TRUE,data.table=FALSE) 
	bed1R <- GRanges(seqnames=bed1$V1,ranges=IRanges(start=bed1$V2, end=bed1$V3),strand=bed1$V4)
	#FPKM1 <- FPKM(bed1R,hetsread1,hets1frag) 
	
	#FPKM for all merged regions:
	counts <- summarizeOverlaps(features=bed1R,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	read_counts <- assay(counts) 
	FPKM <- data.frame(seqnames(bed1R),start(bed1R),width(bed1R),read_counts)
	colnames(FPKM) <- c("chr","start","width","counts")
	FPKM <- mutate(FPKM,fpkm = counts/(hets1frag*width))
	thresh <- filter(FPKM,fpkm>1)
	#the read counts are incorrect somehow...
	zeroes <- which(FPKM$fpkm==0)
	switch <- bed1R[zeroes]
	anti <- strand(switch)
	#length(which(strand(anti)=="*"))
	anti[anti=="+"] <- "*"
	anti[anti=="-"] <- "+"
	anti[anti=="*"] <- "-"
	strand(switch) <- anti
	testcounts <- summarizeOverlaps(features=switch,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	test_counts <- assay(testcounts) 
	#> length(which(test_counts==0))                                 
	#[1] 19934
	#> length(which(test_counts==1))
	#[1] 45970
	switchFPKM <- data.frame(seqnames(switch),start(switch),width(switch),test_counts)
	colnames(switchFPKM) <- c("chr","start","width","counts")
	switchFPKM <- mutate(switchFPKM,fpkm = counts/(hets1frag*width))
	
	file <- '/auto/cmb-00/rr/engie/RNA/merged2.bed'
	bed2 <- fread(file,fill=TRUE,verbose=TRUE,data.table=FALSE) 
	bed2R <- GRanges(seqnames=bed1$V1,ranges=IRanges(start=bed2$V2, end=bed2$V3),strand=bed2$V4)
	counts2 <- summarizeOverlaps(features=bed2R,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	read_counts2 <- assay(counts2) 
	FPKM2 <- data.frame(seqnames(bed2R),start(bed2R),width(bed2R),read_counts2)
	colnames(FPKM2) <- c("chr","start","width","counts")
	FPKM2 <- mutate(FPKM2,fpkm = counts/(hets2frag*width))
	thresh <- filter(FPKM2,fpkm>1)

	file <- '/auto/cmb-00/rr/engie/RNA/merged2_nobookend.bed'
	bed2nb <- fread(file,fill=TRUE,verbose=TRUE,data.table=FALSE) 
	bed2nbR <- GRanges(seqnames=bed2nb$V1,ranges=IRanges(start=bed2nb$V2, end=bed2nb$V3),strand=bed2nb$V4)
	counts2nb <- summarizeOverlaps(features=bed2nbR,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	read_counts2nb <- assay(counts2nb) 
	FPKM2nb <- data.frame(seqnames(bed2nbR),start(bed2nbR),width(bed2nbR),read_counts2nb)
	colnames(FPKM2nb) <- c("chr","start","width","counts")
	FPKM2nb <- mutate(FPKM2nb,fpkm = counts/(hets2frag*width))
	
	
	TPMscale <- TPMScaleFac(hetsread1,bed1R,hets1frag)
	
	clusters <- findReadRegions("chr1",39693870,39800883,"-",bed1R)
	MAML_TPM <- TPMCalc(clusters,hetsread1,TPMscale)
	MAML_TPM2 <- TPMCalc(clusters,hetsread2,TPMscale)

	bach2 <- findReadRegions("chr20",24073750,24223950, "+",bed1R)
	bach2FPKM <- FPKMCalc(bach2,hetsread1,hets1frag)
	bach2FPKMcheck <- FPKMCalc(bach2,hetsread3,hets2frag)

}
