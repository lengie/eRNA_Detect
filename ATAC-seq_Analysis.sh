### ATAC-seq_Analysis.sh
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
  
samtools sort -n sox10ATAC882.rmChrM.rmDup.rmMulti.filtered.shifted.bam -o sox10ATAC882.rmall.sorted.bam -@ 20
samtools sort -n sox10ATAC883.rmChrM.rmDup.rmMulti.filtered.shifted.bam -o sox10ATAC883.rmall.sorted.bam -@ 20

Genrich -t sox10ATAC883.rmall.sorted.bam,sox10ATAC882.rmall.sorted.bam -j -d 100 -v -o sox10ATAC_genrich #this should be sox10ATAC_genrich.narrowPeaks
