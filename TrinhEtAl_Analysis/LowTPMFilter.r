### LowTPMFilter.r
### Plotting average values per row of tabular data (plus or minus strand of RNA-seq counts TPM values), with desired labels
###
### Purpose: Compare average TPM across region of plus strand or minus strand of three different species' replicates
### Input: Tabular files with TPM binned 50bp across 3kb span, plus then minus strand of a single region represented per row
### Output: Scatter plot with average TPM, plus strand x-axis minus strand y-axis, every enhancer region a separate dot, colors by species

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
plus <- sapply(1:nrow(c1_1tab), function(x){mean(as.matrix(c1_1tab[x,1:60]))})
minus <- sapply(1:nrow(c1_1tab), function(x){mean(as.matrix(c1_1tab[x,61:120]))})

toplot <- data.frame(plus=plus,
                     minus=minus,
                     species='chimp')

addRepPlot <- function(tab,sp,plotdf=toplot){
    plusn <- sapply(1:nrow(tab), function(x){mean(as.matrix(tab[x,1:60]))})
    minusn <- sapply(1:nrow(tab), function(x){mean(as.matrix(tab[x,61:120]))})
    new <- rbind(plotdf, data.frame(plus=plusn,minus=minusn,species=sp))  
    return(new)
}

toplot <- addRepPlot(c1_2tab,'chimp',toplot)
toplot <- addRepPlot(c2_1tab,'chimp',toplot)
toplot <- addRepPlot(c2_2tab,'chimp',toplot)
toplot <- addRepPlot(h1_1tab,'human',toplot)
toplot <- addRepPlot(h1_2tab,'human',toplot)
toplot <- addRepPlot(h2_1tab,'human',toplot)
toplot <- addRepPlot(h2_2tab,'human',toplot)
toplot <- addRepPlot(h3_1tab,'human',toplot)
toplot <- addRepPlot(h3_2tab,'human',toplot)
toplot <- addRepPlot(sox1_enh,'zebrafish',toplot)
toplot <- addRepPlot(sox2_enh,'zebrafish',toplot)

                 
                 
png("PrescottTrinh_TPMByStrandBySpecies",width=1200,height=1080)
    ggplot(toplot, aes(x=minus, y=plus,color=species)) +
      geom_point(size=0.25) + xlab("plus") + ylab("minus") +
      ggtitle("Prescott et al Trinh et al avg TPM by strand by species") +
      scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10")
dev.off() 
