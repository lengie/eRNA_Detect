## Build the genome index
hisat2-build -p 20 ../Danio_rerio.Zv9.74.dna.toplevel.fa Zv9

## Alignment of sox10 nuclear datasets
first=FASTQfiles_logs/RNA-seq*71*.R1.fastq
second=FASTQfiles_logs/RNA-seq*71*.R2.fastq
index=/panfs/qcb-panasas/engie/GRCz9or10Analysis/HISAT2Genome/Zv9

hisat2 -p 20 -x $index -1 $first -2 $second --rna-strandness FR > sox10_871_HisatDefault.SAM

first=FASTQfiles_logs/RNA-seq*73*.R1.fastq
second=FASTQfiles_logs/RNA-seq*73*.R2.fastq
hisat2 -p 20 -x $index -1 $first -2 $second --rna-strandness FR > sox10_873_HisatDefault.SAM

## Prepping bigwig files for visualization
samtools view -S -b sox10_871_HisatDefault.SAM > sox10_871_HisatDefault.bam
samtools view -S -@ 20 -b sox10_873_HisatDefault.SAM > sox10_873_HisatDefault.bam

(samtools view -@ 20 -H sox10_871_HisatDefault.bam; samtools view -@ 20 -F 2308 sox10_871_HisatDefault.bam | grep -w 'NH:i:1') | samtools view -@ 20 -bS - > sox10_871_HisatDefault.primary.bam
(samtools view -@ 20 -H sox10_873_HisatDefault.bam; samtools view -@ 20 -F 2308 sox10_873_HisatDefault.bam | grep -w 'NH:i:1') | samtools view -@ 20 -bS - > sox10_873_HisatDefault.primary.bam

# files must be indexed first
# make sure have the right conda env loaded
samtools sort -@ 20 -O bam -o sox10_873_HisatDefault.primarySorted.bam sox10_873_HisatDefault.primary.bam 
samtools sort -@ 20 -O bam -o sox10_871_HisatDefault.primarySorted.bam sox10_871_HisatDefault.primary.bam 
samtools sort -@ 20 -O bam -o sox10_871_HisatDefault.Sorted.bam sox10_871_HisatDefault.bam 

samtools index -@ 20 sox10_871_HisatDefault.primarySorted.bam
samtools index -@ 20 sox10_873_HisatDefault.primarySorted.bam
samtools index -@ 20 sox10_871_HisatDefault.Sorted.bam

bamCoverage -b sox10_871_HisatDefault.Sorted.bam -o sox10_871_HisatDefault.plus.bw --binSize 1 --filterRNAstrand forward --numberOfProcessors max 
bamCoverage -b sox10_871_HisatDefault.Sorted.bam -o sox10_871_HisatDefault.minus.bw --binSize 1 --filterRNAstrand reverse --numberOfProcessors max 

bamCoverage -b sox10_871_HisatDefault.primarySorted.bam -o sox10_871_HisatDefault.primary.plus.bw --binSize 1 --filterRNAstrand forward --numberOfProcessors max 
bamCoverage -b sox10_871_HisatDefault.primarySorted.bam -o sox10_871_HisatDefault.primary.minus.bw --binSize 1 --filterRNAstrand reverse --numberOfProcessors max 
bamCoverage -b sox10_873_HisatDefault.primarySorted.bam -o sox10_873_HisatDefault.primary.plus.bw --binSize 1 --filterRNAstrand forward --numberOfProcessors max 

bamCoverage -b sox10_873_HisatDefault.primarySorted.bam -o sox10_873_HisatDefault.primary.minus.bw --binSize 1 --filterRNAstrand reverse --numberOfProcessors max 

## bactin data sets
first=FASTQfiles_logs/RNA-seq*69*.R1.fastq
second=FASTQfiles_logs/RNA-seq*69*.R2.fastq
index=/panfs/qcb-panasas/engie/GRCz9or10Analysis/HISAT2Genome/Zv9

hisat2 -p 20 -x $index -1 $first -2 $second --rna-strandness FR > bactin_869_HisatDefault.SAM

first=FASTQfiles_logs/RNA-seq*70*.R1.fastq
second=FASTQfiles_logs/RNA-seq*70*.R2.fastq
hisat2 -p 20 -x $index -1 $first -2 $second --rna-strandness FR > bactin_870_HisatDefault.SAM

tobam=bactin*.SAM
for file in $tobam
do
    LABEL=$(echo $file | sed 's/\_HisatDefault.SAM//');
    echo $LABEL
    samtools view -S -b $LABEL'_HisatDefault.SAM' > $LABEL'_HisatDefault.bam'
    (samtools view -@ 20 -H $LABEL'_HisatDefault.bam'; samtools view -@ 20 -F 2308 $LABEL'_HisatDefault.bam' | grep -w 'NH:i:1') | samtools view -@ 20 -bS - > $LABEL'_HisatDefault.primary.bam'
    samtools sort -@ 20 -O bam -o $LABEL'_HisatDefault.primarySorted.bam' $LABEL'_HisatDefault.primary.bam'
    samtools index -@ 20 $LABEL'_HisatDefault.primarySorted.bam'
done


list=*_HisatDefault.primarySorted.bam

for file in $list
do 
    LABEL=$(echo $file | sed 's/\_HisatDefault.primarySorted.bam//');
    echo $LABEL
    bamCoverage -b $file -o $LABEL'_HisatDefault.primaryBPM.plus.bw' \
        --binSize 1 --filterRNAstrand forward --numberOfProcessors max --effectiveGenomeSize 1267788788 --normalizeUsing BPM
    bamCoverage -b $file -o $LABEL'_HisatDefault.primaryBPM.minus.bw' \
        --binSize 1 --filterRNAstrand reverse --numberOfProcessors max --effectiveGenomeSize 1267788788 --normalizeUsing BPM
done
