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

fantomhg <- fread("human_permissive_enhancers_phase_1_and_2.bed")
colnames(fantomhg) <- c("chr","start","end","region","ID","strand","initiationPlus","initiationMinus","thickStart","thickEnd","poolExpPlus","poolExpMinus") #I'm guessing on these a bit based on the readme but I don't need the later columns anyway
fantomgr <- GRanges(fantomhg)
seqlevelsStyle(fantomgrh) <- "ensembl"

# load the original bed files
loadbd <- function(file){
    bd <- fread(file)
    colnames(bd) <- c("chr","start","end","ID","score","strand")
    return(bd)
}
bed <- loadbd(YOUR FILE)

num_enh_bed <- function(file,enhgr){
        bed <- loadbd(paste(file,".bed",sep=""))
	gr <- GRanges(bed)
	seqlevelsStyle(gr) <- "ensembl"
	ov <- findOverlaps(gr,enhgr)
	print(paste("Number of overlaps: ",length(ov),sep=""))
        print(paste("Percentage of false positives:,",length(ov)/nrow(bed),sep=" "))
	all <- enhgr[subjectHits(ov),]                         # this is the list of enhancer regions
        filename <- paste(file,"_FANTOMfalsepos.bed",sep="")
	write.table(all,filename,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
        return(all)
}
	      
bed_FANTOM <- num_enh_bed(bed,fantomgr)

## ENSEMBL
regmm <- fread("mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20201021.gff")
colnames(regmm) <- c("chr","build","feature","start","end","V6","V7","V8","comments")
enhgr <- GRanges(regmm)
seqlevelsStyle(enhgr) <- "ensembl" # I got this dataset from ensembl, but just making sure

reghg <- fread("homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff")
colnames(reghg) <- c("chr","build","feature","start","end","V6","V7","V8","comments")
enhgrh <- GRanges(reghg)
seqlevelsStyle(enhgrh) <- "ensembl" 

num_enh_ensembl <- function(filename,enhgr){
	report <- fread(paste(filename,".bed",sep=""))
	clean <- report[,c(5,13,14)]
	colnames(clean) <- c("chr","start","end")
	gr <- GRanges(clean)
	seqlevelsStyle(gr) <- "ensembl"
	ov <- findOverlaps(gr,enhgr)
	print(paste("Number of overlaps: ",length(ov),sep=""))
	enh <- reg[subjectHits(ov),] %>% dplyr::filter(feature=="enhancer")
        print(paste("Number of enhancers: ",nrow(enh),sep=""))
        print(paste("Percentage of false positives:",nrow(enh)/nrow(clean),sep=" "))
	all <- reg[subjectHits(ov),]
        write.table(all,paste(filename,"_EnsemblFPOverlap.bed",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
	return(all)
}
	      
h0_fn <- num_enh("report_h0_falseneg_12layersBinaryBed3")
