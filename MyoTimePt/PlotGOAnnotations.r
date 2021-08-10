### PlotGOAnnotations.r
###
### Purpose: Gene ontology analyses between conditions: 48hpf, 72hpf, 96hpf
###
###
### Written by Liana Engie
### Last updated: June 2021
###
### runPairGO(bamfile,gtffile)
### Input: 
### Output: bar plot of enriched GO categories 

library(data.table)
library(GOplot)

plotGObar <- function(data,title){
    ggplot(arrange(data[1:20,],over_represented_pvalue), 
                 aes(x = reorder(term,over_represented_pvalue),
                 y = numDEInCat, 
                 fill = over_represented_pvalue)) +
    geom_bar(stat = "identity") + coord_flip() + 
    ggtitle(title) + ylab("Number of DE genes in category") +labs(fill="p-value")
}
 
foursev <- fread("tpm4_devtimepts4872_goseqGOWall.csv")
fournine <- fread("tpm4_devtimepts4896_goseqGOWall.csv")
sevnine <- fread("tpm4_devtimepts7296_goseqGOWall.csv")

plotGObar(foursev,"48hpf vs 72hpf") 
plotGObar(fournine,"48hpf vs 96hpf") 
plotGObar(sevnine,"72hpf vs 96hpf") 
