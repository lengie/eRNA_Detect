### FalsePosNeg_Lookup.r
###
### Written by Liana Engie
### 
###
### Looking up if a bed file of genomic locations has enhancer functionality or is associated with a GWAS SNP

## FANTOM5
fantom <- fread("mouse_permissive_enhancers_phase_1_and_2.bed")
colnames(fantom) <- c("chr","start","end","region","ID","strand","initiationPlus","initiationMinus","thickStart","thickEnd","poolExpPlus","poolExpMinus") #I'm guessing on these a bit based on the readme but I don't need the later columns anyway
fantomgr <- GRanges(fantom)
seqlevelsStyle(fantomgr) <- "ensembl"

# load the original bed files
loadbd <- function(file){
    bd <- fread(file)
    colnames(bd) <- c("chr","start","end","ID","score","strand")
    return(bd)
}
bed <- loadbd(YOUR FILE)

num_enh_bed <- function(bed,enhgr){
	gr <- GRanges(bed)
	seqlevelsStyle(gr) <- "ensembl"
	ov <- findOverlaps(gr,enhgr)
	print(paste("Number of overlaps: ",length(ov),sep=""))
	all <- enhgr[subjectHits(ov),]                         # this is the list of enhancer regions
	return(all)
}

bed_FANTOM <- num_enh_bed(bed,fantomgr)

write.table(bed_FANTOM,YOUR FILE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

## ENSEMBL
reg <- fread("mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20201021.gff")
colnames(reg) <- c("chr","build","feature","start","end","V6","V7","V8","comments")
enhgr <- GRanges(reg)
seqlevelsStyle(enhgr) <- "ensembl" # I got this dataset from ensembl, but just making sure

num_enh <- function(filename){
	report <- fread(filename)
	clean <- report[,c(5,13,14)]
	colnames(clean) <- c("chr","start","end")
	gr <- GRanges(clean)
	seqlevelsStyle(gr) <- "ensembl"
	ov <- findOverlaps(gr,enhgr)
	print(paste("Number of overlaps: ",length(ov),sep=""))
	enh <- reg[subjectHits(ov),] %>% dplyr::filter(feature=="enhancer")
        print(paste("Number of enhancers: ",nrow(enh),sep=""))
	all <- reg[subjectHits(ov),]
	return(all)
}

h0_fn <- num_enh("report_h0_falseneg_12layersBinaryBed3.bed")
write.table(h0_fn,"h0_falseneg_12layersBinary_EnsemblRegOverlap.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
