### ATAC-seq_Processing.sh
###
###
### Code written for analysis of ATAC-seq FASTQ files from Trinh et al 2017
### Last updated: July 2021
###
### Input: fastq files (from DanioCode)
### Output: ATAC-seq peaks called by MACS2
###
### Trinh et al 2017 paper protocols say The paper's supplemental methods says that ATAC-seq data was sequenced, reads were trimmed for quality using sickle (v1.33) 
###     and mapped using bowtie. Bigwig files were generated, and smoothed by a Perl script courtesy of Jim Hughes. Only paired reads with insert sizes larger than 
###     100bp were selected and reads were displaced by +4bp and -5bp and extended to a read length of 100bp. Peak calling was preformed using MACS2 with -nomodel and 
###     -slocal 1000 parameters. Zebrafish Ensembl gene models were extended by 100bp in 5' ofthe TSS to account for gene mis-annotation. Then ATAC-seq peaks not 
###     overlapping with Ensembl-annotated promoter regions or exons were put into the putative cis-regulatory element set. However, the bigwig files have no strand
###     so not sure how to do said analysis. Restarting from their raw files

# obtain the data
wget https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD006430BS/ATAC-seq_Sauka-Spengler_Lab_0007AS.DCD000883SQ.USERvanessachong.R1.fastq.gz --no-check-certificate 
wget https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD006430BS/ATAC-seq_Sauka-Spengler_Lab_0007AS.DCD000883SQ.USERvanessachong.R2.fastq.gz --no-check-certificate
wget https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD006430BS/ATAC-seq_Sauka-Spengler_Lab_0007AS.DCD000882SQ.USERvanessachong.R1.fastq.gz --no-check-certificate
wget https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD006430BS/ATAC-seq_Sauka-Spengler_Lab_0007AS.DCD000882SQ.USERvanessachong.R2.fastq.gz --no-check-certificate

# adapter trimming
NGmerge  -a  -1 ATAC-seq_Sauka-Spengler_Lab_0007AS.DCD000883SQ.USERvanessachong.R1.fastq.gz  -2 ATAC-seq_Sauka-Spengler_Lab_0007AS.DCD000883SQ.USERvanessachong.R2.fastq.gz  \
    -e 40 -o sox10ATAC883_noadapters
NGmerge  -a  -1 ATAC-seq_Sauka-Spengler_Lab_0007AS.DCD000882SQ.USERvanessachong.R1.fastq.gz  -2 ATAC-seq_Sauka-Spengler_Lab_0007AS.DCD000882SQ.USERvanessachong.R2.fastq.gz  \
    -e 40 -o sox10ATAC882_noadapters

## align with Bowtie2
# build the genome index in the folder where the .fa genome file is
bowtie2-build Danio_rerio.Zv9.74.dna.toplevel.fa danRer7
mv *.bt2 ../NCReadDist/Zv9Index

# align -- here we have reads that are 40bp runs
bowtie2 -q -1 sox10ATAC882_noadapters_1.fastq -2 sox10ATAC882_noadapters_2.fastq -x Zv9Index/danRer7 -I 80 -X 2000 --very-sensitive -p 20 \
    | samtools view -bS | samtools sort -o sox10ATAC882.sorted.bam

bowtie -q -1 sox10ATAC883_noadapters_1.fastq -2 sox10ATAC883_noadapters_2.fastq -x Zv9Index/danRer7 -I 80 -X 2000 --very-sensitive -p 20 \
    | samtools view -bS | samtools sort -o sox10ATAC883.sorted.bam

# removing mitochondrial reads
samtools view -@ 20 -h sox10ATAC883.sorted.bam | grep -v chrM | samtools sort -@ 20 -O bam -o sox10ATAC883.rmChrM.bam

samtools view -@ 20 -h sox10ATAC882.sorted.bam | grep -v chrM | samtools sort -@ 20 -O bam -o sox10ATAC882.rmChrM.bam

# removing PCR duplicates
java -XX:ParallelGCThreads=20 -Djava.io.tmpdir=/tmp -jar picard MarkDuplicates \
  QUIET=true INPUT=sox10ATAC883.rmChrM.bam OUTPUT=sox10ATAC883.rmChrM.marked.bam METRICS_FILE=sox10ATAC883.pcrdup.sorted.metrics \
  REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp

samtools view -h -b -F 1024 sox10ATAC883.rmChrM.bam > sox10ATAC883.rmChrM.rmDup.bam


picard MarkDuplicates QUIET=true INPUT=sox10ATAC882.rmChrM.bam OUTPUT=sox10ATAC882.rmChrM.marked.bam METRICS_FILE=sox10ATAC882.pcrdup.sorted.metrics \
  REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp

samtools view -h -b -F 1024 sox10ATAC882.rmChrM.bam > sox10ATAC882.rmChrM.rmDup.bam


# Remove multi-mapped reads (i.e. those with MAPQ < 30, using -q in SAMtools)
samtools view -h -q 30 sox10ATAC882.rmChrM.rmDup.bam > sox10ATAC882.rmChrM.rmDup.rmMulti.bam
samtools view -h -q 30 sox10ATAC883.rmChrM.rmDup.bam > sox10ATAC883.rmChrM.rmDup.rmMulti.bam

# Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804), retain properly paired reads -f 2
samtools view -h -b -F 1804 -f 2 sox10ATAC882.rmChrM.rmDup.rmMulti.bam > sox10ATAC882.rmChrM.rmDup.rmMulti.filtered.bam
samtools view -h -b -F 1804 -f 2 sox10ATAC883.rmChrM.rmDup.rmMulti.bam > sox10ATAC883.rmChrM.rmDup.rmMulti.filtered.bam

# shift reads to account for Tn5 insertion/duplication
# Not critical for the work we are doing here as we do not need super precise boundaries
alignmentSieve --numberOfProcessors 20 --ATACshift --bam sox10ATAC882.rmChrM.rmDup.rmMulti.filtered.bam -o sox10ATAC882.shiftedtmp.bam

# the bam file needs to be sorted again
samtools sort -@ 20 -O bam -o sox10ATAC882.rmChrM.rmDup.rmMulti.filtered.shifted.bam sox10ATAC882.shiftedtmp.bam
samtools index -@ 20 sox10ATAC882.rmChrM.rmDup.rmMulti.filtered.shifted.bam
rm sox10ATAC882.shiftedtmp.bam

# effective genome size was calculated with khmer

# -f BAMPE, use paired-end information
# --keep-dup all, keep all duplicate reads.
macs2 callpeak --nomodel --extsize 100 -g 1267788788 --keep-dup all --cutoff-analysis -n sox10ATAC882 \
  -t sox10ATAC882.rmChrM.rmDup.rmMulti.filtered.shifted.bam 
  
# using genrich as an alternative to MACS  
samtools sort -n sox10ATAC882.rmChrM.rmDup.rmMulti.filtered.shifted.bam -o sox10ATAC882.rmall.sorted.bam -@ 20
samtools sort -n sox10ATAC883.rmChrM.rmDup.rmMulti.filtered.shifted.bam -o sox10ATAC883.rmall.sorted.bam -@ 20

Genrich -t sox10ATAC883.rmall.sorted.bam,sox10ATAC882.rmall.sorted.bam -j -d 100 -v -o sox10ATAC_genrich.narrowPeaks

## some basic comparisons between the ATAC peaks files in R
#library(remotes) 
#remotes::install_github("jokergoo/cotools")
library(cotools) 
library(data.table)
library(GenomicRanges)
library(gridExtra)
library(dplyr) 
library(parallel)
library(factoextra)
library(cluster)
library(GenomicAlignments)
options(scipen=999)

# load the data
genrich <- fread("sox10ATAC_genrich.narrowPeaks") 
dualMACs <- fread("sox10ATAC_MACS2_narrowPeaks.bed")
colnames(dualMACs) <- c("chr","start","end","name","score","strand","fold_enrichment","pval","qval","blockCount")
colnames(genrich) <- c("chr","start","end","name","score","strand","signalValue","pval","qval","peak")
ggen <- GRanges(genrich)
gmacs <- GRanges(dualMACs)

# jaccard distance between the two genomic ranges
genomicCorr.jaccard(ggen,gmacs)

# get some summary statistics
sumstatdf <- function(index){
    width <- index$end - index$start
    print(min(width))
    print(median(width))
    print(mean(width))
    print(max(width))
    print(nrow(index))}
sumstatcol <- function(index){
    print(min(index))
    print(median(index))
    print(mean(index))
    print(max(index))}
    
sumstatdf(genrich)
sumstatdf(dualMACs)
sumstatcol(genrich$pval)
sumstatcol(dualMACs$pval)

compare <- findOverlaps(ggen,gmacs,ignore.strand=TRUE) 
length(compare)
length(unique(queryHits(compare)))
length(unique(subjectHits(compare)))

# k-means on the ATAC-seq regions, features being plus and minus strand read counts

nuc2 <- fread("sox10nuc873_primary.bed")
colnames(nuc2) <- c("chr","start","end","ID","score","strand")
sox10_nuc2 <- GRanges(nuc2)
#seqlevelsStyle(sox10_nuc2) <- "UCSC"


gen_reads <- data.frame(plus=strand_counts(ggen,sox10_nuc2,"+"),
                        minus=strand_counts(ggen,sox10_nuc2,"-"))

macs_reads <- data.frame(plus=strand_counts(gmacs,sox10_nuc2,"+"),
                         minus=strand_counts(gmacs,sox10_nuc2,"-"))

## How many clusters?
# assuming that the columns are features and rows are samples; variance is being calculated over 
var_mean_plot <- function(input,filename){  
    under_var <- apply(input, 1, var)
    under_mean <- apply(input, 1, mean)
    
    png(paste(filename,".png",sep=""))
        plot(log2(under_mean), log2(under_var), pch='.')
        abline(h=log2(50), col='red')
        abline(v=log2(20), col='red')
        text(x=9,y=24, labels="variance > 50 &\n mean > 20", col='red')
    dev.off()
    return(list(under_var,under_mean))
}

gen_vm <- var_mean_plot(gen_reads,"sox10ATAC_genRichVarMeanStrandedReadPlot")
macs_vm <- var_mean_plot(macs_reads,"sox10ATAC_MACS2VarMeanStrandedReadPlot")

sub20g <- gen_reads[which(gen_vm[[1]] > 50 & gen_vm[[2]] > 20),]
sub20m <- macs_reads[which(macs_vm[[1]]> 50 & macs_vm[[2]] > 20),]

for (i in 2:15) wss20[i] <- sum(kmeans(sub20g,centers=i)$withinss)

png("sox10ATAC_WithinGSS_genrich.png")
plot(1:15, wss20, type="b", xlab="Number of Clusters",
  ylab="Within groups sum of squares",main="GenRich ATAC Var>50 MeanReads>20")
dev.off()

wss20m <- (nrow(sub20m)-1)*sum(apply(sub20m,2,var))

for (i in 2:15) wss20m[i] <- sum(kmeans(sub20m,centers=i)$withinss)

png("sox10ATAC_WithinGSS_MACS2.png")
plot(1:15, wss20m, type="b", xlab="Number of Clusters",
  ylab="Within groups sum of squares",main="MACS2 ATAC Var>50 MeanReads>20")
dev.off()

## clustering
fit_cluster_g <- clara(gen_reads,k=3,metric="euclidean")
fit_cluster_m <- clara(macs_reads,k=3,metric="euclidean")
fit_cluster_m6 <- clara(macs_reads,k=6,metric="euclidean")

# read out the clusters
gen_reads$cluster <- fit_cluster_g$clustering
macs_reads$cluster <- fit_cluster_m$clustering
macs_reads$cluster6 <- fit_cluster_m6$clustering
 
# plotting
g <- ggplot(gen_reads,aes(x=reads,y=reads.1,color=cluster))+geom_point() +
      ggtitle("GenRich ATAC peaks, k-means k=3") + ylab("plus strand read counts") + xlab("minus strand read counts") + 
      scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10')
m <- ggplot(macs_reads,aes(x=reads,y=reads.1,color=cluster)) + geom_point() +
      ggtitle("MACS2 ATAC peaks, k-means k=3") + ylab("plus strand read counts") + 
      xlab("minus strand read counts")+ scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10')
m6 <- ggplot(macs_reads,aes(x=reads,y=reads.1,color=cluster6)) + geom_point()+ggtitle("MACS2 ATAC peaks, k-means k=6") + ylab("plus strand read counts") + xlab("minus strand read counts")+ scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10')

png("sox10ATAC_kmeans.png",width=860,height=540)
    grid.arrange(g,m,m6,ncol=3)
dev.off()


