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
