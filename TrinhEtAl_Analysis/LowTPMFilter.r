### LowTPMFilter.r
### Plotting average values per row of tabular data (plus or minus strand of RNA-seq counts TPM values), with desired labels
###
### Purpose: Compare average TPM across region of plus strand or minus strand of three different species' replicates
### Input: Tabular files with TPM binned 50bp across 3kb span, plus then minus strand of a single region represented per row
### Output: Scatter plot with average TPM, plus strand y-axis minus strand x-axis, every enhancer region a separate dot, colors by species

library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(data.table)
library(dplyr)
library(ggplot2)
library(fitdistrplus)
library(ballgown)
options(scipen=999)

# if there were no reads, I had deepTools output an NA so that would appear white in the heatmaps. Convert to 0 for the program
c1_1tab <- fread("Prescott2015/chimp_1_1_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
c1_2tab <- fread("Prescott2015/chimp_1_2_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
c2_1tab <- fread("Prescott2015/chimp_2_1_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
c2_2tab <- fread("Prescott2015/chimp_2_2_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h1_1tab <- fread("Prescott2015/human_1_1_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h1_2tab <- fread("Prescott2015/human_1_2_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h2_1tab <- fread("Prescott2015/human_2_1_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h2_2tab <- fread("Prescott2015/human_2_2_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h3_1tab <- fread("Prescott2015/human_3_1_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h3_2tab <- fread("Prescott2015/human_3_2_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
sox1_enh <- fread("DLTesting/sox10_871_bidirATACnoribo_pt1rpb_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
sox2_enh <- fread("DLTesting/sox10_873_bidirATACnoribo_pt1rpb_3kb_50bpTEST.tabular",header=FALSE) %>% replace(is.na(.), 0)
sox1_no <- fread("DLTesting/sox10_871_bidirUnder500bpMerged800bpNoATACOverlap.tabular",header=FALSE) %>% replace(is.na(.), 0)
sox2_no <- fread("DLTesting/sox10_873_bidirUnder500bpMerged800bpNoATACOverlap.tabular",header=FALSE) %>% replace(is.na(.), 0)
c1_1no <- fread("Prescott2015/chimp_1_1_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
c1_2no <- fread("Prescott2015/chimp_1_2_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
c2_1no <- fread("Prescott2015/chimp_2_1_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
c2_2no <- fread("Prescott2015/chimp_2_2_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h1_1no <- fread("Prescott2015/human_1_1_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h1_2no <- fread("Prescott2015/human_1_2_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h2_1no <- fread("Prescott2015/human_2_1_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h2_2no <- fread("Prescott2015/human_2_2_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h3_1no <- fread("Prescott2015/human_3_1_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h3_2no <- fread("Prescott2015/human_3_2_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)

# I want every single row to be a new point
plus <- sapply(1:nrow(c1_1tab), function(x){average(c1_1tab[x,1;120])}
minus <- sapply(1:nrow(c1_1tab), function(x){average(c1_1tab[x,121;240])}

plot <- data.frame(plus=plus,
                   minus=minus
                   species='chimp')

addRepPlot <- function(tab,species){
    
    new <- cbind(plot,
    return(new)
}

plot <- addRepPlot(,'human')
