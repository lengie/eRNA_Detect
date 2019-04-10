### BidirNoncodRNA.r
###
### Purpose: Load RNA-seq data as BAM, remove coding regions via GTF annotations, merge strands into overlapped regions (get TPM of these regions?), explore these. 
###
###
### Written by Liana Engie
### Last updated: April 2019
###
### bidirncRNA(bamfile,gtffile)
### Input: string chromosome number, int input_start, int input_end, string strand (either "+" or "-")
### Output:

library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments) 
library(ggplot2)
library(dplyr) 
library(data.table)
library(bedr)

bamfile <- "/auto/cmb-00/rr/engie/RNA/hets1.bam" 
gtffile <- "/auto/cmb-00/rr/engie/RNA/Danio_rerio.GRCz11.96.gtf" 
fragno <- 50.731011 

bidirncRNAwGTF{
    txdb <- makeTxDbFromGFF(gtffile,
                            format="gtf",
                            circ_seqs = character()
                           )
	coding <- cds(txdb)
	#class(coding)

	flag <- scanBamFlag(isSecondaryAlignment=FALSE, isDuplicate=FALSE)
    bamread <- readGAlignmentPairs(bamfile,
                                   param=ScanBamParam(flag=flag)
                                  )
    #class(bamread)
	
	ncbam <- setdiff(bamread,coding,
                     ignore.strand=FALSE)
	#ignoring for now the * strand reads
    #ncbam then creates a list of bamread that does not include any regions in coding
    
    ncbed <- data.table(chr <- as.character(seqnames(ncbam)),
                        start <- start(ncbam)-1,
                        end <- end(ncbam),
                        strand <- strand(ncbam)
                       )

    
}

