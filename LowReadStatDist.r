## Written by Liana Engie, Dr. Scott E Fraser lab
## Last updated 1/25/2021

# Investigating the statistical distribution of nuclear RNA reads
# Uses sox10 nuclear RNA pulled from Trinh et al 2017 Biotagging paper: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89670

library(Rsamtools)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)

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
