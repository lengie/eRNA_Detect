### detectEnhancers.R
###
### Purpose: Load RNA-seq data, output locations and FPKMs of regions with bidirectional reads within a given region 
###
###
### Written by Liana Engie
### Last updated: September 2018
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
			byStrandTemp <- split(ranges(mergedbedchrsplit$i),by.strand(mergedbedchrsplit$i))
			compare <- mergeByOverlaps(byStrandTemp[[1]],byStrandTemp[[2]])
			overlap <- data.frame(overlap,compare)
		}
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
	
	#separate by strand and merge overlaps
	byStrand <- split(ranges(hetsread1), as.factor(strand(hetsread1)))
	export.bed(byStrand[[1]],name="Hets1PosStrand")
	export.bed(byStrand[[2]],name="Hets1NegStrand")
	pos1 <- as.data.frame(byStrand[[1]])
	neg1 <- as.data.frame(byStrand[[2]])
	write.csv(pos1,file="Hets1PosStrand")
	write.csv(neg1,file="Hets1NegStrand")
	rm(hetsread1)
	compared <- mergeByOverlaps(byStrand[[1]],byStrand[[2]])
	
	
	#Will use gtf file to ID nc regions later for validation, but calc scaling factor with merge  
	#g <- "/auto/cmb-00/rr/engie/RNA/Danio_rerio.GRCz10.87.chr.gtf" 
	#txdb <- makeTxDbFromGFF(g, format="gtf",circ_seqs = character())
	#sequence <- seqlevels(txdb)
	#new <- mapSeqlevels(sequence,"UCSC")
	#new <- new[complete.cases(new)]
	#txdb <- renameSeqlevels(txdb,new)
	#STILL NEED TO PULL OUT NON CODING REGIONS FROM GTF
	#negative rather than a category itself
	
	#bam file merged with bookends, else '/auto/cmb-00/rr/engie/RNA/merged1_nobookend.bed'
	file <- '/auto/cmb-00/rr/engie/RNA/merged1.bed'
	bed1 <- fread(file,fill=TRUE,verbose=TRUE,data.table=FALSE) 
	bed1R <- GRanges(seqnames=bed1$V1,ranges=IRanges(start=bed1$V2, end=bed1$V3),strand=bed1$V4)
	#FPKM1 <- FPKMCalc(bed1R,hetsread1,hets1frag) 
	
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
	
	#TODO 3/15
	file <- '/auto/cmb-00/rr/engie/RNA/merged1_nobookend.bed'
	bed1nb <- fread(file,fill=TRUE,verbose=TRUE,data.table=FALSE) 
	bed1nbR <- GRanges(seqnames=bed1nb$V1,ranges=IRanges(start=bed1nb$V2, end=bed1nb$V3),strand=bed1nb$V4)
	counts <- summarizeOverlaps(features=bed1nbR,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	nbread_counts <- assay(counts) 
	FPKMnb <- data.frame(seqnames(bed1nbR),start(bed1nbR),width(bed1nbR),nbread_counts)
	colnames(FPKMnb) <- c("chr","start","width","counts")
	FPKMnb <- mutate(FPKMnb,fpkm = counts/(hets1frag*width))
	thresh <- filter(FPKMnb,fpkm>1)
	zeroes <- which(FPKM$fpkm==0)
	switch <- bed1nbR[zeroes]
	anti <- strand(switch)
	#length(which(strand(anti)=="*"))
	anti[anti=="+"] <- "*"
	anti[anti=="-"] <- "+"
	anti[anti=="*"] <- "-"
	strand(switch) <- anti
	testcounts <- summarizeOverlaps(features=switch,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	test_counts <- assay(testcounts) 
	switchFPKMnb <- data.frame(seqnames(switch),start(switch),width(switch),test_counts)
	colnames(switchFPKMnb) <- c("chr","start","width","counts")
	switchFPKMnb <- mutate(switchFPKMnb,fpkm = counts/(hets1frag*width))
	
	file <- '/auto/cmb-00/rr/engie/RNA/merged2.bed'
	bed2 <- fread(file,fill=TRUE,verbose=TRUE,data.table=FALSE) 
	bed2R <- GRanges(seqnames=bed2$V1,ranges=IRanges(start=bed2$V2, end=bed2$V3),strand=bed2$V4)
	counts2 <- summarizeOverlaps(features=bed2R,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	read_counts2 <- assay(counts2) 
	FPKM2 <- data.frame(seqnames(bed2R),start(bed2R),width(bed2R),read_counts2)
	colnames(FPKM2) <- c("chr","start","width","counts")
	FPKM2 <- mutate(FPKM2,fpkm = counts/(hets2frag*width))
	#> length(which(FPKM2==0))
	#[1] 35542
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
	switchFPKM <- data.frame(seqnames(switch),start(switch),width(switch),test_counts)
	colnames(switchFPKM) <- c("chr","start","width","counts")
	switchFPKM <- mutate(switchFPKM,fpkm = counts/(hets1frag*width))
	
	file <- '/auto/cmb-00/rr/engie/RNA/merged2_nobookend.bed'
	bed2nb <- fread(file,fill=TRUE,verbose=TRUE,data.table=FALSE) 
	bed2nbR <- GRanges(seqnames=bed2nb$V1,ranges=IRanges(start=bed2nb$V2, end=bed2nb$V3),strand=bed2nb$V4)
	counts2nb <- summarizeOverlaps(features=bed2nbR,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	read_counts2nb <- assay(counts2nb) 
	FPKM2nb <- data.frame(seqnames(bed2nbR),start(bed2nbR),width(bed2nbR),read_counts2nb)
	colnames(FPKM2nb) <- c("chr","start","width","counts")
	FPKM2nb <- mutate(FPKM2nb,fpkm = counts/(hets2frag*width))

	#> length(start(bed2R))
	#[1] 40839
	#> length(start(bed2nbR))
	#[1] 245760
	#> length(start(bed1R))
	#[1] 275101
	#> length(seqnames(bed1nbR))
	#[1] 277035
	#> length(seqnames(hetsread2))
	#[1] 49779418
	
	hets26 <- "/auto/cmb-00/rr/engie/RNA/26hr1.bam" 
	hets26read1 <- readGAlignmentPairs(hets26,param=ScanBamParam(flag=flag))

	#using the featureCounts/Subreads output from terminal instead
	#This is no bookend and multimapping reads
	file <- "/auto/cmb-00/rr/engie/RNA/output_strand_noB.tsv"
	hets1_mm <- fread(file,fill=TRUE,data.table=FALSE,skip=1,header=TRUE) 
	colnames(hets1_mm)[7] <- "Counts"
	hets1_mm <- mutate(hets1_mm,FPKM = Counts/(hets1frag*Length))
	hets1_mm <- mutate(hets1_mm,RPK = Counts/Length)
	scaling <- sum(hets1_mm$RPK)*hets1frag
	[1] 362997.2
	hets1_mm <- mutate(hets1_mm, TPM = Counts/scaling)
	write.table(hets1_mm,file="hets1_mm_FPKMTPM.csv",row.names=FALSE) 
	
	#Try no bookend and NO multimapping reads
	file <- "/auto/cmb-00/rr/engie/RNA/hets1_noB_noMultim.tsv"
	hets1 <- fread(file,fill=TRUE,data.table=FALSE,skip=1,header=TRUE) 
	colnames(hets1)[7] <- "Counts"
	hets1 <- mutate(hets1,FPKM = Counts/(hets1frag*Length))
	hets1 <- mutate(hets1,RPK = Counts/Length)
	scaling <- sum(hets1$RPK)*hets1frag
	#[1] 241875.2
	hets1 <- mutate(hets1, TPM = Counts/scaling)
	write.table(hets1,file="hets1_FPKMTPM.csv",row.names=FALSE) 
	
	file <- "/auto/cmb-00/rr/engie/RNA/hets2_noB_noMultim.tsv"
	hets2 <- fread(file,fill=TRUE,data.table=FALSE,skip=1,header=TRUE) 
	hets2frag <- 56.315336 
	colnames(hets2)[7] <- "Counts"
	hets2 <- mutate(hets2,FPKM = Counts/(hets2frag*Length))
	hets2 <- mutate(hets2,RPK = Counts/Length)
	scaling2 <- sum(hets2$RPK)*hets2frag
	#[1] 251975.5
	hets2 <- mutate(hets2, TPM = Counts/scaling2)
	write.table(hets2,file="hets2_FPKMTPM.csv",row.names=FALSE) 
	
	file <- "/auto/cmb-00/rr/engie/RNA/hets2_noB.tsv"
	hets2_mm <- fread(file,fill=TRUE,data.table=FALSE,skip=1,header=TRUE) 
	colnames(hets2_mm)[7] <- "Counts"
	hets2_mm <- mutate(hets2_mm,FPKM = Counts/(hets2frag*Length))
	hets2_mm <- mutate(hets2_mm,RPK = Counts/Length)
	scaling2 <- sum(hets2_mm$RPK)*hets2frag
	#[1] 402954.9
	hets2_mm <- mutate(hets2_mm, TPM = Counts/scaling2)
	write.table(hets2_mm,file="hets2_mm_FPKMTPM.csv",row.names=FALSE) 
	
	
	

	
	#brute force test on chr 21 RAI14-RXFP3
	rai <- GRanges(seqnames="chr21",ranges=IRanges(start=seq(19492351,19497301,by=50),end=seq(19492401,19497351,by=50)),strand="-")
	rai_counts <- summarizeOverlaps(features=rai,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	rai4k_minus <- assay(rai_counts)
	rai <- GRanges(seqnames="chr21",ranges=IRanges(start=seq(19492351,19497301,by=50),end=seq(19492401,19497351,by=50)),strand="+")
	rai_counts <- summarizeOverlaps(features=rai,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	rai4k_plus <- assay(rai_counts)
	qplot(1:length(rai4k_minus),rai4k_plus)
	qplot(1:length(rai4k_minus),rai4k_minus)
	
	rai <- GRanges(seqnames="chr21",ranges=IRanges(start=seq(19492351,19497301,by=50),end=seq(19492401,19497351,by=50)),strand="-")
	rai_counts <- summarizeOverlaps(features=rai,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	rai4k_minus2 <- assay(rai_counts)
	rai <- GRanges(seqnames="chr21",ranges=IRanges(start=seq(19492351,19497301,by=50),end=seq(19492401,19497351,by=50)),strand="+")
	rai_counts <- summarizeOverlaps(features=rai,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	rai4k_plus2 <- assay(rai_counts)
	qplot(1:length(rai4k_minus2),rai4k_plus2)
	qplot(1:length(rai4k_minus2),rai4k_minus2)
	
	#brute force test on chr 6 BRINP3-RGS8 
	FAM <- GRanges(seqnames="chr6",ranges=IRanges(start=seq(35573980,35575020,by=10),end=seq(35573990,35575030,by=10)),strand="+")
	fam_counts <- summarizeOverlaps(features=FAM,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	fam75k_plus <- assay(fam_counts)
	FAM <- GRanges(seqnames="chr6",ranges=IRanges(start=seq(35573980,35575020,by=10),end=seq(35573990,35575030,by=10)),strand="-")
	fam_counts <- summarizeOverlaps(features=FAM,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	fam75k_minus <- assay(fam_counts)
	qplot(1:length(fam75k_minus),fam75k_plus)
	qplot(1:length(fam75k_minus),fam75k_minus)	
	FAM <- GRanges(seqnames="chr6",ranges=IRanges(start=c(35573980,35573980),end=c(35574460,35574460)),strand=c("-","+"))
	fam_counts <- summarizeOverlaps(features=FAM,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	read_count <- assay(fam_counts)
	fam_FPKM <- read_count/(hets1frag*580) 
	#Actual sequence is 6:35573980..35574560
		
	FAM <- GRanges(seqnames="chr6",ranges=IRanges(start=seq(35573980,35575020,by=10),end=seq(35573990,35575030,by=10)),strand="+")
	fam_counts <- summarizeOverlaps(features=FAM,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	fam75k_plus2 <- assay(fam_counts)
	FAM <- GRanges(seqnames="chr6",ranges=IRanges(start=seq(35573980,35575020,by=10),end=seq(35573990,35575030,by=10)),strand="-")
	fam_counts <- summarizeOverlaps(features=FAM,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	fam75k_minus2 <- assay(fam_counts)
	qplot(1:length(fam75k_minus2),fam75k_plus2)
	qplot(1:length(fam75k_minus2),fam75k_minus2)
	
	#try getting counts when I include multimapping?
	hets1 <- "/auto/cmb-00/rr/engie/RNA/hets1.bam" 
	flag <- scanBamFlag(isSecondaryAlignment=TRUE,isDuplicate=FALSE)
	hetsread1 <- readGAlignmentPairs(hets1,param=ScanBamParam(flag=flag))
	
	#Calculating FPKM of the subreads tsv file from the original merge
	#Let's try reading the multipmapped .tsv file and calculating FPKMs from that?
	file <- '/auto/cmb-00/rr/engie/RNA/output_strand_noB.tsv'
	tsv <- fread(file,skip=1,header=TRUE,data.table=FALSE)
	tsv <- mutate(tsv,fpkm = hets1.bam/(hets1frag*Length))
	
	#geom_freqpoly(mapping = NULL, data = NULL, stat = "bin",position = "identity", ..., na.rm = FALSE, show.legend = NA,inherit.aes = TRUE)

	ggplot(tsv,aes(fpkm)) + geom_histogram(bins = 100, show.legend = NA, inherit.aes = TRUE)
	#position = "stack", ..., binwidth = NULL, na.rm = FALSE, stat = "bin", 
	
	
	#Using grep to remove multimapped reads
	file <- "/auto/cmb-00/rr/engie/RNA/merged1_nob_noMultim.bed"
	bed1 <- fread(file,fill=TRUE,verbose=TRUE,data.table=FALSE) 
	bed1R <- GRanges(seqnames=bed1$V1,ranges=IRanges(start=bed1$V2, end=bed1$V3),strand=bed1$V4)
	#length(start(bed1R)) == 277105
	hets1 <- "/auto/cmb-00/rr/engie/RNA/hets1_noMultim.bam" 
	hets1frag <- 50.731011 
	flag <- scanBamFlag(isSecondaryAlignment=FALSE, isDuplicate=FALSE)
	#galp <- readGAlignmentPairs(hets1,param=ScanBamParam(flag=flag))	
	hetsread1 <- readGAlignmentPairs(hets1,param=ScanBamParam(flag=flag))
	counts <- summarizeOverlaps(features=bed1R,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	read_counts <- assay(counts) 
	TPM <- data.frame(seqnames(bed1R),start(bed1R),width(bed1R),read_counts)
	colnames(TPM) <- c("chr","start","width","counts")
	TPM <- mutate(TPM,fpkm = counts/(hets1frag*width))
	TPM <- mutate(TPM,rpk = counts/width)
	scaling <- sum(TPM$rpk)*hets1frag
	[1] 172468
	TPM <- mutate(TPM, tpm = counts/172468)
	ggplot(TPM,aes(tpm)) + geom_histogram(bins = 100, show.legend = NA, inherit.aes = TRUE)
	
	
	int_test <- GRanges(seqnames="chr2",ranges=IRanges(start=seq(14822009,14842009,by=10),end=seq(14822019,14842019,by=10)),strand="+")
	test_counts <- summarizeOverlaps(features=int_test,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	pde_3_plus <- assay(test_counts)
	int_test <- GRanges(seqnames="chr2",ranges=IRanges(start=seq(14822009,14842009,by=10),end=seq(14822019,14842019,by=10)),strand="-")
	test_counts <- summarizeOverlaps(features=int_test,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	pde_3_minus <- assay(test_counts)
	qplot(1:length(pde_3_plus),pde_3_plus)
	qplot(1:length(pde_3_plus),pde_3_minus)
	int_test <- GRanges(seqnames="chr2",ranges=IRanges(start=seq(14822009,14842009,by=10),end=seq(14822019,14842019,by=10)),strand="+")
	test_counts <- summarizeOverlaps(features=int_test,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	pde_3_plus2 <- assay(test_counts)
	int_test <- GRanges(seqnames="chr2",ranges=IRanges(start=seq(14822009,14842009,by=10),end=seq(14822019,14842019,by=10)),strand="-")
	test_counts <- summarizeOverlaps(features=int_test,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	pde_3_minus2 <- assay(test_counts)
	qplot(1:length(pde_3_plus2),pde_3_plus2)
	qplot(1:length(pde_3_plus2),pde_3_minus2)	
	
	pde4b <- GRanges(seqnames="chr6",ranges=IRanges(start=seq(30939636,30940376,by=10),end=seq(30939646,30940386,by=10)),strand="-")
	pde_counts <- summarizeOverlaps(features=pde4b,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	pde6_3_minus <- assay(pde_counts)
	pde4b <- GRanges(seqnames="chr6",ranges=IRanges(start=seq(30939636,30940376,by=10),end=seq(30939646,30940386,by=10)),strand="+")
	pde_counts <- summarizeOverlaps(features=pde4b,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	pde6_3_plus <- assay(pde_counts)
	qplot(1:length(pde6_3_minus),pde6_3_minus)
	qplot(1:length(pde6_3_minus),pde6_3_plus)
	pde4b <- GRanges(seqnames="chr6",ranges=IRanges(start=seq(30939636,30940376,by=10),end=seq(30939646,30940386,by=10)),strand="-")
	pde_counts <- summarizeOverlaps(features=pde4b,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	pde6_3_minus2 <- assay(pde_counts)
	pde4b <- GRanges(seqnames="chr6",ranges=IRanges(start=seq(30939636,30940376,by=10),end=seq(30939646,30940386,by=10)),strand="+")
	pde_counts <- summarizeOverlaps(features=pde4b,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	pde6_3_plus2 <- assay(pde_counts)
	qplot(1:length(pde6_3_minus2),pde6_3_minus2)
	qplot(1:length(pde6_3_minus2),pde6_3_plus2)
	
	gpc <-GRanges(seqnames="chr1",ranges=IRanges(start=seq(2863290,2864080,by=10),end=seq(2863300,2864090,by=10)),strand="-")
	gpc_counts <- summarizeOverlaps(features=gpc,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	gpca_3_minus <- assay(gpc_counts)
	gpc <-GRanges(seqnames="chr1",ranges=IRanges(start=seq(2863290,2864080,by=10),end=seq(2863300,2864090,by=10)),strand="+")
	gpc_counts <- summarizeOverlaps(features=gpc,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	gpca_3_plus <- assay(gpc_counts)
	qplot(1:length(gpca_3_minus),gpca_3_minus)
	qplot(1:length(gpca_3_minus),gpca_3_plus)
	gpc <-GRanges(seqnames="chr1",ranges=IRanges(start=seq(2863290,2864080,by=10),end=seq(2863300,2864090,by=10)),strand="-")
	gpc_counts <- summarizeOverlaps(features=gpc,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	gpca_3_minus2 <- assay(gpc_counts)
	gpc <-GRanges(seqnames="chr1",ranges=IRanges(start=seq(2863290,2864080,by=10),end=seq(2863300,2864090,by=10)),strand="+")
	gpc_counts <- summarizeOverlaps(features=gpc,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	gpca_3_plus2 <- assay(gpc_counts)
	qplot(1:length(gpca_3_minus2),gpca_3_minus2)
	qplot(1:length(gpca_3_minus2),gpca_3_plus2)
		
	tbx <- GRanges(seqnames="chr9",ranges=IRanges(start=seq(21101600,21107600,by=50),end=seq(21101650,21107650,by=50)),strand="+")
	tbx_counts <- summarizeOverlaps(features=tbx,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx_plus <- assay(tbx_counts)
	tbx <- GRanges(seqnames="chr9",ranges=IRanges(start=seq(21101600,21107600,by=50),end=seq(21101650,21107650,by=50)),strand="-")
	tbx_counts <- summarizeOverlaps(features=tbx,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx_minus <- assay(tbx_counts)
	qplot(1:length(tbx_minus),tbx_minus)
	qplot(1:length(tbx_minus),tbx_plus)
	tbx <- GRanges(seqnames="chr9",ranges=IRanges(start=seq(21101600,21107600,by=50),end=seq(21101650,21107650,by=50)),strand="+")
	tbx_counts <- summarizeOverlaps(features=tbx,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx_plus2 <- assay(tbx_counts)
	tbx <- GRanges(seqnames="chr9",ranges=IRanges(start=seq(21101600,21107600,by=50),end=seq(21101650,21107650,by=50)),strand="-")
	tbx_counts <- summarizeOverlaps(features=tbx,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx_minus2 <- assay(tbx_counts)
	qplot(1:length(tbx_minus2),tbx_minus2)
	qplot(1:length(tbx_minus2),tbx_plus2)
	
	#intergenic
	Stx <- GRanges(seqnames="chr14",ranges=IRanges(start=seq(233990,241590,by=50),end=seq(234040,241640,by=50)),strand="+")
	stx_counts <- summarizeOverlaps(features=Stx,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	stx_plus <- assay(stx_counts)
	Stx <- GRanges(seqnames="chr14",ranges=IRanges(start=seq(233990,241590,by=50),end=seq(234040,241640,by=50)),strand="-")
	stx_counts <- summarizeOverlaps(features=Stx,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	stx_minus <- assay(stx_counts)
	qplot(1:length(stx_minus),stx_plus)
	qplot(1:length(stx_minus),stx_minus)
	Stx <- GRanges(seqnames="chr14",ranges=IRanges(start=seq(233990,241590,by=50),end=seq(234040,241640,by=50)),strand="+")
	stx_counts <- summarizeOverlaps(features=Stx,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	stx_plus2 <- assay(stx_counts)
	Stx <- GRanges(seqnames="chr14",ranges=IRanges(start=seq(233990,241590,by=50),end=seq(234040,241640,by=50)),strand="-")
	stx_counts <- summarizeOverlaps(features=Stx,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	stx_minus2 <- assay(stx_counts)
	qplot(1:length(stx_minus2),stx_plus2)
	qplot(1:length(stx_minus2),stx_minus2)
	
	#4th intron
	int_test <- GRanges(seqnames="chr1",ranges=IRanges(start=seq(23142190,23143830,by=10),end=seq(23142200,23143840,by=10)),strand="-")
	test_counts <- summarizeOverlaps(features=int_test,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	slit_4_minus <- assay(test_counts)
	int_test <- GRanges(seqnames="chr1",ranges=IRanges(start=seq(23142190,23143830,by=10),end=seq(23142200,23143840,by=10)),strand="+")
	test_counts <- summarizeOverlaps(features=int_test,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	slit_4_plus <- assay(test_counts)
	qplot(1:length(slit_4_minus),slit_4_minus)
	qplot(1:length(slit_4_minus),slit_4_plus)
	int_test <- GRanges(seqnames="chr1",ranges=IRanges(start=seq(23142190,23143830,by=10),end=seq(23142200,23143840,by=10)),strand="-")
	test_counts <- summarizeOverlaps(features=int_test,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	slit_4_minus2 <- assay(test_counts)
	int_test <- GRanges(seqnames="chr1",ranges=IRanges(start=seq(23142190,23143830,by=10),end=seq(23142200,23143840,by=10)),strand="+")
	test_counts <- summarizeOverlaps(features=int_test,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	slit_4_plus2 <- assay(test_counts)
	qplot(1:length(slit_4_minus2),slit_4_minus2)
	qplot(1:length(slit_4_minus2),slit_4_plus2)

	#intergenic
	ATXN <- GRanges(seqnames="chr10",ranges=IRanges(start=seq(3360303,3363303,by=50), end=seq(3360353,3363353,by=50)),strand="+")
	atxn_counts <- summarizeOverlaps(features=ATXN,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	atxn_plus <- assay(atxn_counts)
	ATXN <- GRanges(seqnames="chr10",ranges=IRanges(start=seq(3360303,3363303,by=50), end=seq(3360353,3363353,by=50)),strand="-")
	atxn_counts <- summarizeOverlaps(features=ATXN,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	atxn_minus <- assay(atxn_counts)
	qplot(1:length(atxn_minus),atxn_minus)
	qplot(1:length(atxn_plus),atxn_plus)
	ATXN <- GRanges(seqnames="chr10",ranges=IRanges(start=seq(3360303,3363303,by=50), end=seq(3360353,3363353,by=50)),strand="+")
	atxn_counts <- summarizeOverlaps(features=ATXN,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	atxn_plus2 <- assay(atxn_counts)
	ATXN <- GRanges(seqnames="chr10",ranges=IRanges(start=seq(3360303,3363303,by=50), end=seq(3360353,3363353,by=50)),strand="-")
	atxn_counts <- summarizeOverlaps(features=ATXN,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	atxn_minus2 <- assay(atxn_counts)
	qplot(1:length(atxn_minus2),atxn_minus2)
	qplot(1:length(atxn_plus2),atxn_plus2)
	
	#9th intron
	hspa <- GRanges(seqnames="chr10",ranges=IRanges(start=seq(29966864,29967054,by=10),end=seq(29966874,29967064,by=10)),strand="-")
	hspa_counts <- summarizeOverlaps(features=hspa,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	hspa_plus <- assay(hspa_counts)
	hspa <- GRanges(seqnames="chr10",ranges=IRanges(start=seq(29966864,29967054,by=10),end=seq(29966874,29967064,by=10)),strand="+")
	hspa_counts <- summarizeOverlaps(features=hspa,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	hspa_minus <- assay(hspa_counts)
	qplot(1:length(hspa_minus),hspa_plus)
	qplot(1:length(hspa_minus),hspa_minus)
	hspa <- GRanges(seqnames="chr10",ranges=IRanges(start=seq(29966864,29967054,by=10),end=seq(29966874,29967064,by=10)),strand="-")
	hspa_counts <- summarizeOverlaps(features=hspa,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	hspa_plus2 <- assay(hspa_counts)
	hspa <- GRanges(seqnames="chr10",ranges=IRanges(start=seq(29966864,29967054,by=10),end=seq(29966874,29967064,by=10)),strand="+")
	hspa_counts <- summarizeOverlaps(features=hspa,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	hspa_minus2 <- assay(hspa_counts)
	qplot(1:length(hspa_minus2),hspa_plus2)
	qplot(1:length(hspa_minus2),hspa_minus2)
	
	#first intron
	ptpn <- GRanges(seqnames="chr10",ranges=IRanges(start=seq(3428000,3448650,by=50),end=seq(3428050,3448700,by=50)),strand="-")
	ptpn_counts <- summarizeOverlaps(features=ptpn,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	ptpn_minus2 <- assay(ptpn_counts)
	ptpn <- GRanges(seqnames="chr10",ranges=IRanges(start=seq(3428000,3448650,by=50),end=seq(3428050,3448700,by=50)),strand="+")
	ptpn_counts <- summarizeOverlaps(features=ptpn,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	ptpn_plus2 <- assay(ptpn_counts)
	qplot(1:length(ptpn_plus2),ptpn_minus2)
	qplot(1:length(ptpn_minus2),ptpn_plus2)
	
	#Testing from Smemo and Dickel papers: TBX5
	tbx5 <- GRanges(seqnames="chr5",ranges=IRanges(start=seq(71394233,71599334,by=50),end=seq(71394283,71599384,by=50)),strand="+")
	tbx_counts <- summarizeOverlaps(features=tbx5,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx5_plus <- assay(tbx_counts)
	tbx5 <- GRanges(seqnames="chr5",ranges=IRanges(start=seq(71394233,71599334,by=50),end=seq(71394283,71599384,by=50)),strand="-")
	tbx_counts <- summarizeOverlaps(features=tbx5,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx5_minus <- assay(tbx_counts)
	qplot(1:length(tbx5_minus),tbx5_minus)
	qplot(1:length(tbx5_minus),tbx5_plus)
	tbx5 <- GRanges(seqnames="chr5",ranges=IRanges(start=seq(71394233,71599334,by=50),end=seq(71394283,71599384,by=50)),strand="+")
	tbx_counts <- summarizeOverlaps(features=tbx5,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx5_plus2 <- assay(tbx_counts)
	tbx5 <- GRanges(seqnames="chr5",ranges=IRanges(start=seq(71394233,71599334,by=50),end=seq(71394283,71599384,by=50)),strand="-")
	tbx_counts <- summarizeOverlaps(features=tbx5,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx5_minus2 <- assay(tbx_counts)
	qplot(1:length(tbx5_minus2),tbx5_minus2)
	qplot(1:length(tbx5_minus2),tbx5_plus2)
	
	tbx5 <- GRanges(seqnames="chr5",ranges=IRanges(start=seq(71505142,71532242,by=50),end=seq(71505192,71532292,by=50)),strand="+")
	tbx_counts <- summarizeOverlaps(features=tbx5,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx5_plus <- assay(tbx_counts)
	tbx5 <- GRanges(seqnames="chr5",ranges=IRanges(start=seq(71505142,71532242,by=50),end=seq(71505192,71532292,by=50)),strand="-")
	tbx_counts <- summarizeOverlaps(features=tbx5,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx5_minus <- assay(tbx_counts)
	qplot(1:length(tbx5_minus),tbx5_minus)
	qplot(1:length(tbx5_minus),tbx5_plus)
	tbx5 <- GRanges(seqnames="chr5",ranges=IRanges(start=seq(71505142,71532242,by=50),end=seq(71505192,71532292,by=50)),strand="+")
	tbx_counts <- summarizeOverlaps(features=tbx5,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx5_plus2 <- assay(tbx_counts)
	tbx5 <- GRanges(seqnames="chr5",ranges=IRanges(start=seq(71505142,71532242,by=50),end=seq(71505192,71532292,by=50)),strand="-")
	tbx_counts <- summarizeOverlaps(features=tbx5,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx5_minus2 <- assay(tbx_counts)
	qplot(1:length(tbx5_minus2),tbx5_minus2)
	qplot(1:length(tbx5_minus2),tbx5_plus2) 
	
	tbx5b <- GRanges(seqnames="chr5",ranges=IRanges(start=seq(22643860,22672660,by=50),end=seq(22643910,22672710,by=50)),strand="+")
	tbx_counts <- summarizeOverlaps(features=tbx5b,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx5_plus <- assay(tbx_counts)
	tbx5 <- GRanges(seqnames="chr5",ranges=IRanges(start=seq(22643860,22672660,by=50),end=seq(22643910,22672710,by=50)),strand="-")
	tbx_counts <- summarizeOverlaps(features=tbx5,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx5_minus <- assay(tbx_counts)
	qplot(1:length(tbx5_minus),tbx5_minus)
	qplot(1:length(tbx5_minus),tbx5_plus)
	tbx5 <- GRanges(seqnames="chr5",ranges=IRanges(start=seq(22643860,22672660,by=50),end=seq(22643910,22672710,by=50)),strand="+")
	tbx_counts <- summarizeOverlaps(features=tbx5,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx5_plus2 <- assay(tbx_counts)
	tbx5b <- GRanges(seqnames="chr5",ranges=IRanges(start=seq(22643860,22672660,by=50),end=seq(22643910,22672710,by=50)),strand="-")
	tbx_counts <- summarizeOverlaps(features=tbx5b,reads=hetsread2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	tbx5_minus2 <- assay(tbx_counts)
	qplot(1:length(tbx5_minus2),tbx5_minus2)
	qplot(1:length(tbx5_minus2),tbx5_plus2)	
	
	#using the above functions
	TPMscale <- TPMScaleFac(hetsread1,bed1R,hets1frag)
	
	clusters <- findReadRegions("chr1",39693870,39800883,"-",bed1R)
	MAML_TPM <- TPMCalc(clusters,hetsread1,TPMscale)
	MAML_TPM2 <- TPMCalc(clusters,hetsread2,TPMscale)

	bach2 <- findReadRegions("chr20",24073750,24223950, "+",bed1R)
	bach2FPKM <- FPKMCalc(bach2,hetsread1,hets1frag)
	bach2FPKMcheck <- FPKMCalc(bach2,hetsread3,hets2frag)

	
}

	HSHeart <- "/auto/cmb-00/rr/engie/RNA/GRCh38.illumina.heart.1.bam" 
	#this is paired
	flag <- scanBamFlag(isSecondaryAlignment=FALSE, isDuplicate=FALSE)
	HumHeart <- readGAlignments(HSHeart,param=ScanBamParam(flag=flag))
	
	file <- '/auto/cmb-00/rr/engie/RNA/26hr1.bg'
	bg26 <- fread(file,fill=TRUE,verbose=TRUE,data.table=FALSE)
	head(bg26)
	unique(bg26$1)
	bg26_1R <- GRanges(seqnames=bed1$V1,ranges=IRanges(start=bed1$V2, end=bed1$V3),score=bed1$V4)
	
	
