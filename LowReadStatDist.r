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
library(ballgown)
library(BSgenome.Drerio.UCSC.danRer7)
library(TFBSTools) 
library(motifmatchr)
options(scipen=999)

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
gtffile <- "/auto/cmb-00/rr/engie/GCF_000002035.4_Zv9_genomic.gtf"
txdb <- makeTxDbFromGFF(gtffile,
                        format="gtf",
                        circ_seqs = character()
                        )
seqlevelsStyle(txdb) <- "UCSC"
transcripts <- transcripts(txdb)
exons <- fread("Zv9exons.bed") #from UCSC Table Browser
colnames(exons) <- c("seqnames","start","end","ID","score","strand")
utr5 <- fread("Zv95UTR.bed")
utr3 <- fread("Zv93UTR.bed")
colnames(utr5) <- c("seqnames","start","end","ID","score","strand")
colnames(utr3) <- c("seqnames","start","end","ID","score","strand")

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


## extracting bidirectional regions in RNA-seq
utr5minus <- dplyr::filter(utr5,strand=="-")
utr5plus <- dplyr::filter(utr5,strand=="+")
utr5minus <- data.frame(seqnames=utr5minus$seqnames,
			start=utr5minus$start,
			end=utr5minus$end+500,
			strand=utr5minus$strand)
utr5plus <- data.frame(seqnames=utr5plus$seqnames,
		       start=utr5plus$start-500,
		       end=utr5plus$end,
		       strand=utr5plus$strand)

utr3minus <- dplyr::filter(utr3,strand=="-")
utr3plus <- dplyr::filter(utr3,strand=="+")
utr3minus <- data.frame(seqnames=utr3minus$seqnames,
			start=utr3minus$start-500,
			end=utr3minus$end,
			strand=utr3minus$strand)
utr3plus <- data.frame(seqnames=utr3minus$seqnames,
		       start=utr3minus$start,
		       end=utr3minus$end+500,
		       strand=utr3minus$strand)

coding <- c(GRanges(exons),GRanges(utr5plus),GRanges(utr5minus),GRanges(utr3plus),GRanges(utr3minus))

overlaps <- findOverlaps(sox10_nuc1,coding,ignore.strand=FALSE) 
hits1 <- sox10_nuc1[-queryHits(overlaps),]
underten1 <- hits1[width(hits1)<10000]
ten1 <- hits1[width(hits1)>10000]
write.table(ten1,file="sox10_Zv9nuc1FlankedNoncodingOver10k.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

overlaps2 <- findOverlaps(sox10_nuc2,coding,ignore.strand=FALSE) 
hits2 <- sox10_nuc2[-queryHits(overlaps2),]
underten2 <- hits2[width(hits2)<10000]
ten2 <- hits2[width(hits2)>10000]
write.table(ten2,file="sox10_Zv9nuc2FlankedNoncodingOver10k.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

underten1df <- data.frame(chr = as.character(seqnames(underten1)),
                    	  start = start(underten1),
                    	  end = end(underten1),
               		  name = "n/a",
			  score = 0,
			  strand = strand(underten1)
			  )
write.table(underten1df,file="sox10_Zv9nuc1FlankedNoncodingUnder10k.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

underten2df <- data.frame(chr = as.character(seqnames(underten2)),
                    	  start = start(underten2),
                    	  end = end(underten2),
                	  name = "n/a",
			  score = 0,
			  strand = strand(underten2)
			  )
write.table(underten2df,file="sox10_Zv9nuc2FlankedNoncodingUnder10k.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

# exit R
bedtools sort -i sox10_Zv9nuc1FlankedNoncodingUnder10k.bed > sox10_Zv9nuc1FlankedNoncodingUnder10kSorted.bed
bedtools merge -i sox10_Zv9nuc1FlankedNoncodingUnder10kSorted.bed -d 0 -S - > sox10_Zv9nuc1FlankedNoncodingUnder10kMinusOnlyMerged.bed
bedtools merge -i sox10_Zv9nuc1FlankedNoncodingUnder10kSorted.bed -d 0 -S + > sox10_Zv9nuc1FlankedNoncodingUnder10kPlusOnlyMerged.bed

bedtools sort -i sox10_Zv9nuc2FlankedNoncodingUnder10k.bed > sox10_Zv9nuc2FlankedNoncodingUnder10kSorted.bed
bedtools merge -i sox10_Zv9nuc2FlankedNoncodingUnder10kSorted.bed -d 0 -S - > sox10_Zv9nuc2FlankedNoncodingUnder10kMinusOnlyMerged.bed
bedtools merge -i sox10_Zv9nuc2FlankedNoncodingUnder10kSorted.bed -d 0 -S + > sox10_Zv9nuc2FlankedNoncodingUnder10kPlusOnlyMerged.bed

# back in R
plus <- fread("sox10_Zv9nuc1FlankedNoncodingUnder10kPlusOnlyMerged.bed",data.table=TRUE,fill=TRUE)
minus <- fread("sox10_Zv9nuc1FlankedNoncodingUnder10kMinusOnlyMerged.bed",data.table=TRUE,fill=TRUE)
plus2 <- fread("sox10_Zv9nuc2FlankedNoncodingUnder10kPlusOnlyMerged.bed",data.table=TRUE,fill=TRUE)
minus2 <- fread("sox10_Zv9nuc2FlankedNoncodingUnder10kMinusOnlyMerged.bed",data.table=TRUE,fill=TRUE)
colnames(minus) <- c("seqnames","start","end","strand")
colnames(plus) <- c("seqnames","start","end","strand")
colnames(minus2) <- c("seqnames","start","end","strand")
colnames(plus2) <- c("seqnames","start","end","strand")
plus <- GRanges(plus)
minus <- GRanges(minus)
plus2 <- GRanges(plus2)
minus2 <- GRanges(minus2)
nuc1merge <- mergeByOverlaps(minus,plus,ignore.strand=TRUE)
nuc2merge <- mergeByOverlaps(minus2,plus2,ignore.strand=TRUE)

nuc1merged <- data.frame(chr=c(seqnames(nuc1merge$plus),seqnames(nuc1merge$minus)),
                        start=c(start(nuc1merge$plus),start(nuc1merge$minus)),
                        end=c(end(nuc1merge$plus),end(nuc1merge$minus)),
                        ID="n/a",
                        score=0,
                        strand=c(strand(nuc1merge$plus),strand(nuc1merge$minus)))
nuc2merged <- data.frame(chr=c(seqnames(nuc2merge$plus2),seqnames(nuc2merge$minus2)),
                        start=c(start(nuc2merge$plus2),start(nuc2merge$minus2)),
                        end=c(end(nuc2merge$plus2),end(nuc2merge$minus2)),
                        ID="n/a",
                        score=0,
                        strand=c(strand(nuc2merge$plus2),strand(nuc2merge$minus2)))

write.table(nuc1merged,file="sox10_Zv9nuc1FlankedNoncodingUnder10kOverlapsToMerge.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
write.table(nuc2merged,file="sox10_Zv9nuc2FlankedNoncodingUnder10kOverlapsToMerge.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

# outside of R
bedtools sort -i sox10_Zv9nuc1FlankedNoncodingUnder10kOverlapsToMerge.bed > sox10_Zv9nuc1FlankedNoncodingUnder10kOverlapsToMergeSorted.bed
bedtools merge -i sox10_Zv9nuc1FlankedNoncodingUnder10kOverlapsToMergeSorted.bed -d 0 > sox10_Zv9nuc1FlankedNoncodingUnder10kOverlapsMerged.bed

bedtools sort -i sox10_Zv9nuc2FlankedNoncodingUnder10kOverlapsToMerge.bed > sox10_Zv9nuc2FlankedNoncodingUnder10kOverlapsToMergeSorted.bed
bedtools merge -i sox10_Zv9nuc2FlankedNoncodingUnder10kOverlapsToMergeSorted.bed -d 0 > sox10_Zv9nuc2FlankedNoncodingUnder10kOverlapsMerged.bed

# back in R
nuc1bidir <- fread("sox10_Zv9nuc1FlankedNoncodingUnder10kOverlapsMerged.bed")
nuc2bidir <- fread("sox10_Zv9nuc2FlankedNoncodingUnder10kOverlapsMerged.bed")
colnames(nuc1bidir) <- c("chr","start","end")
colnames(nuc2bidir) <- c("chr","start","end")	

feat <- GRanges(seqnames=nuc1bidir$chr,ranges=IRanges(start=nuc1bidir$start,end=nuc1bidir$end))  
bidir_read1 <- summarizeOverlaps(features=feat,reads=sox10_nuc1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=TRUE) 
bidir_read_counts1 <- assay(bidir_read1)
nuc1bidir <- dplyr::mutate(nuc1bidir, size = end-start)
nuc1bidir <- cbind(nuc1bidir,bidir_read_counts1,strand="+")
write.table(nuc1bidir,"sox10_Zv9nuc1FlankedWiderBidirRegionsIDisScore.bed",quote=FALSE,row.names=FALSE,col.names=FALSE)

feat2 <- GRanges(seqnames=nuc2bidir$chr,ranges=IRanges(start=nuc2bidir$start,end=nuc2bidir$end))  
bidir_read2 <- summarizeOverlaps(features=feat2,reads=sox10_nuc2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=TRUE) 
bidir_read_counts2 <- assay(bidir_read2)
nuc2bidir <- dplyr::mutate(nuc2bidir, size = end-start)
nuc2bidir <- cbind(nuc2bidir,bidir_read_counts2,strand="+")
write.table(nuc2bidir,"sox10_Zv9nuc2FlankedBidirRegionsIDisScore.bed",quote=FALSE,row.names=FALSE,col.names=FALSE)

#plus strand
feat1 <- GRanges(seqnames=nuc1bidir$chr,ranges=IRanges(start=nuc1bidir$start,end=nuc1bidir$end),strand="+")
feat2 <- GRanges(seqnames=nuc2bidir$chr,ranges=IRanges(start=nuc2bidir$start,end=nuc2bidir$end),strand="+")  

bidir_sox1plus <- summarizeOverlaps(features=feat1,reads=sox10_nuc1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=FALSE) 
bidir_sox10_counts1plus <- assay(bidir_sox1plus)
bidir_sox2plus <- summarizeOverlaps(features=feat2,reads=sox10_nuc2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=FALSE) 
bidir_sox10_counts2plus <- assay(bidir_sox2plus)

#minus strand
feat1 <- GRanges(seqnames=nuc1bidir$chr,ranges=IRanges(start=nuc1bidir$start,end=nuc1bidir$end),strand="-")  
feat2 <- GRanges(seqnames=nuc2bidir$chr,ranges=IRanges(start=nuc2bidir$start,end=nuc2bidir$end),strand="-")  

bidir_sox1minus <- summarizeOverlaps(features=feat1,reads=sox10_nuc1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=FALSE) 
bidir_sox10_counts1minus <- assay(bidir_sox1minus)
bidir_sox2minus <- summarizeOverlaps(features=feat2,reads=sox10_nuc2,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE,ignore.strand=FALSE) 
bidir_sox10_counts2minus <- assay(bidir_sox2minus)

counts1 <- cbind(nuc1bidir[,c(1,2,3,8,7)],plus=bidir_sox10_counts1plus, minus=bidir_sox10_counts1minus)
counts1 <- dplyr::mutate(counts1, size = end-start)
counts2 <- cbind(nuc2bidir[,c(1,2,3,8,7)],plus=bidir_sox10_counts2plus, minus=bidir_sox10_counts2minus)
counts2 <- dplyr::mutate(counts2, size = end-start)

write.table(counts1,file="sox10_Zv9nuc1FlankedBidirRegionsStrandedCounts.bed",quote=FALSE, row.names=FALSE,col.names=FALSE)
write.table(counts2,file="sox10_Zv9nuc2FlankedBidirRegionsStrandedCounts.bed",quote=FALSE, row.names=FALSE,col.names=FALSE)

#overlaps comparison
allten <- fread("ATACCrossRef/sox10_BiotaggingAllClusters.bed")
colnames(allten) <- c("chr","start","end","size","strand","counts","size")
allclusters <- GRanges(allten)

sox10_1 <- fread('ATACCrossRef/Sox10nuclear_cluster1.counts.txt')
sox10_2 <- fread('ATACCrossRef/Sox10nuclear_cluster2.counts.txt')
sox10_3 <- fread('ATACCrossRef/Sox10nuclear_cluster3.counts.txt')
sox10_4 <- fread('ATACCrossRef/Sox10nuclear_cluster4.counts.txt')
sox10_5 <- fread('ATACCrossRef/Sox10nuclear_cluster5.counts.txt')
sox10_6 <- fread('ATACCrossRef/Sox10nuclear_cluster6.counts.txt')
sox10_7 <- fread('ATACCrossRef/Sox10nuclear_cluster7.counts.txt')
sox10_8 <- fread('ATACCrossRef/Sox10nuclear_cluster8.counts.txt')
sox10_9 <- fread('ATACCrossRef/Sox10nuclear_cluster9.counts.txt')
sox10_10 <- fread('ATACCrossRef/Sox10nuclear_cluster10.counts.txt')

cluster1 <- GRanges(sox10_1)
cluster2 <- GRanges(sox10_2)
cluster3 <- GRanges(sox10_3)
cluster4 <- GRanges(sox10_4)
cluster5 <- GRanges(sox10_5)
cluster6 <- GRanges(sox10_6)
cluster7 <- GRanges(sox10_7)
cluster8 <- GRanges(sox10_8)
cluster9 <- GRanges(sox10_9)
cluster10 <- GRanges(sox10_10)

chrplus <- dplyr::filter(plusbig,grepl("chr",chr))
chrminus <- dplyr::filter(minusbig,grepl("chr",chr))
chrequal <- dplyr::filter(equal,grepl("chr",chr))
chrplus2 <- dplyr::filter(plusbig2,grepl("chr",chr))
chrminus2 <- dplyr::filter(minusbig2,grepl("chr",chr))
chrequal2 <- dplyr::filter(equal2,grepl("chr",chr))

gplus <- GRanges(chrplus)
gminus <- GRanges(chrminus)
gequal <- GRanges(chrequal)
gplus2 <- GRanges(chrplus2)
gminus2 <- GRanges(chrminus2)
gequal2 <- GRanges(chrequal2)

pctOverlap(gplus,allclusters)
pctOverlap(gequal,allclusters)
pctOverlap(gminus,allclusters)
pctOverlap(gplus2,allclusters)
pctOverlap(gequal2,allclusters)
pctOverlap(gminus2,allclusters)

## TFBSs for the bidirectional regions
motifs <- readJASPARMatrix("../JASPAR2020CoreTransfac/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt", matrixClass="PFM")
plusmatchscorecalc <- matchMotifs(motifs,gplus,genome="danRer7",out="score")
plusmatches <- motifMatches(plusmatchscorecalc)
plusscores <- assay(plusmatchscorecalc)

minusmatchscorecalc <- matchMotifs(motifs,gminus,genome="danRer7",out="score")
minusmatches <- motifMatches(minusmatchscorecalc)
minusscores <- assay(minusmatchscorecalc)

equalmatchscorecalc <- matchMotifs(motifs,gequal,genome="danRer7",out="score")
eqmatches <- motifMatches(equalmatchscorecalc)
eqscores <- assay(equalmatchscorecalc)


plusmatchscorecalc2 <- matchMotifs(motifs,gplus2,genome="danRer7",out="score")
plusmatches2 <- motifMatches(plusmatchscorecalc2)
plusscores2 <- assay(plusmatchscorecalc2)

minusmatchscorecalc2 <- matchMotifs(motifs,gminus2,genome="danRer7",out="score")
minusmatches2 <- motifMatches(minusmatchscorecalc2)
minusscores2 <- assay(minusmatchscorecalc2)

equalmatchscorecalc2 <- matchMotifs(motifs,gequal2,genome="danRer7",out="score")
eqmatches2 <- motifMatches(equalmatchscorecalc2)
eqscores2 <- assay(equalmatchscorecalc2)

sums <- rowSums(plusmatches)
zero <- which(sums==0)
rem <- chrplus[-zero,]
sumsm <- rowSums(minusmatches)
zerom <- which(sumsm==0)
remminus <- chrminus[-zerom,]
sumeq <- rowSums(eqmatches)
zeroeq <- which(sumeq==0)
remeq <- chrequal[-zeroeq,]

sums2 <- rowSums(plusmatches2)
zero2 <- which(sums2==0)
remp2 <- chrplus2[-zero2,]
sumsm2 <- rowSums(minusmatches2)
zerom2 <- which(sumsm2==0)
remminus2 <- chrminus2[-zerom2,]
sumeq2 <- rowSums(eqmatches2)
zeroeq2 <- which(sumeq2==0)
remeq2 <- chrequal[-zeroeq2,]

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
exon.norm.params <- fitdistr(exonNorm$sox10nuc1,"normal")$estimate
exon.poisson.params <- fitdistr(exonNorm$sox10nuc1,"poisson")$estimate
exon.negbinom.params <- fitdistr(exonNorm$sox10nuc1,"negative binomial", lower=c(0,0), method = "SANN")$estimate
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

tx.norm.params <- fitdistr(txNorm$sox10nuc1,"normal")$estimate
tx.poisson.params <- fitdistr(txNorm$sox10nuc1,"poisson")$estimate
tx.negbinom.params <- fitdistr(txNorm$sox10nuc1,"negative binomial", lower=c(0,0), method = "SANN")$estimate
tx.dist.params <- map(list(Normal = tx.norm.params,Poisson = tx.poisson.params,`Negative Binomial` = tx.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

tx2.norm.params <- fitdistr(txNorm$sox10nuc2,"normal")$estimate
tx2.poisson.params <- fitdistr(txNorm$sox10nuc2,"poisson")$estimate
tx2.negbinom.params <- fitdistr(txNorm$sox10nuc2,"negative binomial", lower=c(0,0), method = "SANN")$estimate
tx2.dist.params <- map(list(Normal = tx2.norm.params,Poisson = tx2.poisson.params,`Negative Binomial` = tx2.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

txA.norm.params <- fitdistr(txNorm$sox10polyA,"normal")$estimate
txA.poisson.params <- fitdistr(txNorm$sox10polyA,"poisson")$estimate
txA.negbinom.params <- fitdistr(txNorm$sox10polyA,"negative binomial", lower=c(0,0), method = "SANN")$estimate
txA.dist.params <- map(list(Normal = txA.norm.params,Poisson = txA.poisson.params,`Negative Binomial` = txA.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

nc.norm.params <- fitdistr(normnctb$sox10nuc2,"normal")$estimate
nc.poisson.params <- fitdistr(normnctb$sox10nuc1,"poisson")$estimate
nc.negbinom.params <- fitdistr(normnctb$sox10nuc1,"negative binomial", lower=c(0,0), method = "SANN")$estimate
nc.dist.params <- map(list(Normal = nc.norm.params,Poisson = nc.poisson.params,`Negative Binomial` = nc.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

nc2.norm.params <- fitdistr(normnctb$sox10nuc2,"normal")$estimate
nc2.poisson.params <- fitdistr(normnctb$sox10nuc2,"poisson")$estimate
nc2.negbinom.params <- fitdistr(normnctb$sox10nuc2,"negative binomial", lower=c(0,0), method = "SANN")$estimate
nc2.dist.params <- map(list(Normal = nc2.norm.params,Poisson = nc2.poisson.params,`Negative Binomial` = nc2.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

ncA.norm.params <- fitdistr(normnctb$sox10polyA,"normal")$estimate
ncA.poisson.params <- fitdistr(normnctb$sox10polyA,"poisson")$estimate
ncA.negbinom.params <- fitdistr(normnctb$sox10polyA,"negative binomial", lower=c(0,0), method = "SANN")$estimate
ncA.dist.params <- map(list(Normal = ncA.norm.params,Poisson = ncA.poisson.params,`Negative Binomial` = ncA.negbinom.params),
    ~ map2(names(.),.,~ paste0(.x," = ",round(.y,2))) %>% unlist %>% paste0(.,collapse = ", ")) %>% 
    map2_chr(names(.),., ~ paste(.x,.y,sep=":\n"))

# plots
mybinwidth <- 1

exon1 <- ggplot(data=exonNorm, aes(x=sox10nuc1)) + geom_histogram(binwidth=1) + xlim(0,100)+labs(x="sox10 nuc1 exon read counts") + ylim(0,30000)+
            stat_function(aes(color = "black"),
                          fun=function(x,mean,sd) mybinwidth * nrow(exonNorm) * dnorm(x,mean, sd),
                          args=fitdistr(exonNorm$sox10polyA,"normal")$estimate) +
            stat_function(aes(color = "blue"),
                          fun=function(x,lambda) mybinwidth * nrow(exonNorm) * dpois(x,lambda),
                          args=fitdistr(exonNorm$sox10polyA,"poisson")$estimate,
                          xlim=c(1,100), n=100) + 
            stat_function(aes(color = "orange"),
                          fun=function(x,size, mu) mybinwidth * nrow(exonNorm) * dnbinom(x,size = size, mu = mu),
                          args=fitdistr(exonNorm$sox10polyA,"negative binomial", method="SANN")$estimate,
                          lower=c(0,0),n=100) + 
            scale_color_manual("Distribution", 
                                values=c(black="black",blue="blue",orange="orange"), labels=exon.dist.params)
exon2 <- ggplot(data=exonNorm, aes(x=sox10nuc2)) + geom_histogram(binwidth=1) + xlim(0,100)+labs(x="sox10 nuc2 exon read counts") + ylim(0,30000)+
            stat_function(aes(color = "black"),
                          fun=function(x,mean,sd) mybinwidth * nrow(exonNorm) * dnorm(x,mean, sd),
                          args=fitdistr(exonNorm$sox10polyA,"normal")$estimate) +
            stat_function(aes(color = "blue"),
                          fun=function(x,lambda) mybinwidth * nrow(exonNorm) * dpois(x,lambda),
                          args=fitdistr(exonNorm$sox10polyA,"poisson")$estimate,
                          xlim=c(1,100), n=100) + 
            stat_function(aes(color = "orange"),
                          fun=function(x,size, mu) mybinwidth * nrow(exonNorm) * dnbinom(x,size = size, mu = mu),
                          args=fitdistr(exonNorm$sox10polyA,"negative binomial", method="SANN")$estimate,
                          lower=c(0,0),n=100) + 
            scale_color_manual("Distribution", 
                                values=c(black="black",blue="blue",orange="orange"), labels=exon2.dist.params)
exonA <- ggplot(data=exonNorm, aes(x=sox10polyA)) + geom_histogram(binwidth=1) + xlim(0,100)+labs(x="sox10 polyA exon read counts") + ylim(0,30000)+
            stat_function(aes(color = "black"),
                          fun=function(x,mean,sd) mybinwidth * nrow(exonNorm) * dnorm(x,mean, sd),
                          args=fitdistr(exonNorm$sox10polyA,"normal")$estimate) +
            stat_function(aes(color = "blue"),
                          fun=function(x,lambda) mybinwidth * nrow(exonNorm) * dpois(x,lambda),
                          args=fitdistr(exonNorm$sox10polyA,"poisson")$estimate,
                          xlim=c(1,100), n=100) + 
            stat_function(aes(color = "orange"),
                          fun=function(x,size, mu) mybinwidth * nrow(exonNorm) * dnbinom(x,size = size, mu = mu),
                          args=fitdistr(exonNorm$sox10polyA,"negative binomial", method="SANN")$estimate,
                          lower=c(0,0),n=100) + 
            scale_color_manual("Distribution", 
                                values=c(black="black",blue="blue",orange="orange"), labels=exonA.dist.params)

tx1 <- ggplot(data=txNorm, aes(x=sox10nuc1)) + geom_histogram(binwidth=1) + xlim(0,1000)+labs(x="sox10 nuc1 transcript read counts") + ylim(0,1000)+
            stat_function(aes(color = "black"),
                          fun=function(x,mean,sd) mybinwidth * nrow(txNorm) * dnorm(x,mean, sd),
                          args=fitdistr(txNorm$sox10polyA,"normal")$estimate) +
            stat_function(aes(color = "blue"),
                          fun=function(x,lambda) mybinwidth * nrow(txNorm) * dpois(x,lambda),
                          args=fitdistr(txNorm$sox10polyA,"poisson")$estimate,
                          xlim=c(1,100), n=100) + 
            stat_function(aes(color = "orange"),
                          fun=function(x,size, mu) mybinwidth * nrow(txNorm) * dnbinom(x,size = size, mu = mu),
                          args=fitdistr(txNorm$sox10polyA,"negative binomial", method="SANN")$estimate,
                          lower=c(0,0),n=100) + 
            scale_color_manual("Distribution", values=c(black="black",blue="blue",orange="orange"), labels=tx.dist.params)
tx2 <- ggplot(data=txNorm, aes(x=sox10nuc2)) + geom_histogram(binwidth=1) + xlim(0,1000)+labs(x="sox10 nuc2 transcript read counts") + ylim(0,1000)+ 
            stat_function(aes(color = "black"),
                          fun=function(x,mean,sd) mybinwidth * nrow(txNorm) * dnorm(x,mean, sd),
                          args=fitdistr(txNorm$sox10polyA,"normal")$estimate) +
            stat_function(aes(color = "blue"),
                          fun=function(x,lambda) mybinwidth * nrow(txNorm) * dpois(x,lambda),
                          args=fitdistr(txNorm$sox10polyA,"poisson")$estimate,
                          xlim=c(1,100), n=100) + 
            stat_function(aes(color = "orange"),
                          fun=function(x,size, mu) mybinwidth * nrow(txNorm) * dnbinom(x,size = size, mu = mu),
                          args=fitdistr(txNorm$sox10polyA,"negative binomial", method="SANN")$estimate,
                          lower=c(0,0),n=100) + 
            scale_color_manual("Distribution", values=c(black="black",blue="blue",orange="orange"), labels=tx2.dist.params)
txA <- ggplot(data=txNorm, aes(x=sox10polyA)) + geom_histogram(binwidth=1) + xlim(0,1000)+labs(x="sox10 polyA transcript read counts") + ylim(0,1000)+
            stat_function(aes(color = "black"),
                          fun=function(x,mean,sd) mybinwidth * nrow(txNorm) * dnorm(x,mean, sd),
                          args=fitdistr(txNorm$sox10polyA,"normal")$estimate) +
            stat_function(aes(color = "blue"),
                          fun=function(x,lambda) mybinwidth * nrow(txNorm) * dpois(x,lambda),
                          args=fitdistr(txNorm$sox10polyA,"poisson")$estimate,
                          xlim=c(1,100), n=100) + 
            stat_function(aes(color = "orange"),
                          fun=function(x,size, mu) mybinwidth * nrow(txNorm) * dnbinom(x,size = size, mu = mu),
                          args=fitdistr(txNorm$sox10polyA,"negative binomial", method="SANN")$estimate,
                          lower=c(0,0),n=100) + 
            scale_color_manual("Distribution",values=c(black="black",blue="blue",orange="orange"), labels=txA.dist.params)

lnc1 <- ggplot(data=lncNorm, aes(x=sox10nuc1)) + geom_histogram(binwidth=1) + xlim(0,100)+labs(x="sox10 nuc1 lncNRA read counts") + ylim(0,9000)+
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
lnc2 <- ggplot(data=lncNorm, aes(x=sox10nuc2)) + geom_histogram(binwidth=1) + xlim(0,100)+labs(x="sox10 nuc2 lncNRA read counts") + ylim(0,9000)+
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
chisq.test(exonNorm$sox10nuc1, p=fit)
fit <- dnbinom(0:length,size=0.1616598,mu=6.1709287)
chisq.test(exonNorm$sox10nuc2, p=fit)
fit <- dnbinom(0:length,size=0.0653513,mu=3.1287361)
chisq.test(exonNorm$sox10polyA, p=fit)

length <- nrow(txNorm)-1
fit <- dnbinom(0:length,size=0.3650183,mu=728.9520676)
chisq.test(txNorm$sox10nuc1, p=fit)
fit <- dnbinom(0:length,size=0.3067151,mu=525.8956677)
chisq.test(txNorm$sox10nuc2, p=fit)
fit <- dnbinom(0:length,size=0.1818739,mu=264.9459049)
chisq.test(txNorm$sox10polyA, p=fit)

length <- nrow(lncNorm)-1
fit <- dnbinom(0:length,size=0.1921967,mu=47.5004140)
chisq.test(normnctb$sox10nuc1, p=fit)
chisq.test(normnctb$sox10nuc2, p=fit)
fit <- dnbinom(0:length, size=0.06592968, mu=42.50703446)
chisq.test(normnctb$sox10polyA, p=fit)

