## Written by Liana Engie, Dr. Scott E Fraser lab
## Last updated 1/25/2021

# Investigating the statistical distribution of nuclear RNA reads
# Uses sox10 nuclear RNA pulled from Trinh et al 2017 Biotagging paper: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89670
# Loads the data, combines both strands into single GRanges object, loads annotations, makes count tables

library(data.table)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicAlignments)

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

# combining into GRangesLists
grl <- GRangesList(sox10_nuc1minus,sox10_nuc1plus)
grl2 <- GRangesList(sox10_nuc2minus,sox10_nuc2plus)

# turn the GRangesLists into GRanges
sox10_nuc1 <- unlist(grl)
sox10_nuc2 <- unlist(grl2)

# retrieve the annotations to generate count tables
gtffile <- "/panfs/qcb-panasas/engie/GRCz11EnhDet/Danio_rerio.GRCz11.99.gtf"
txdb <- makeTxDbFromGFF(gtffile,
                        format="gtf",
                        circ_seqs = character()
                        )
seqlevelsStyle(txdb) <- "UCSC"
transcripts <- transcripts(txdb)
exons <- exons(txdb)

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
