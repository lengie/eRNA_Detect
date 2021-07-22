

comp <- function(gr1,gr2){
    overlaps <- findOverlaps(gr1,gr2,ignore.strand=TRUE)
    hits <- gr1[unique(queryHits(overlaps)),]
    print(length(hits))
    return(hits)
}

gr_to_df <- function(gr){
    df <- data.frame(chr=seqnames(gr),
                    start = start(gr),
                    end = end(gr),
                    ID = 1:length(gr),
                    score=0,
                    strand=strand(gr))
    return(df)
}

sumstatcol <- function(index){
    print(min(index))
    print(median(index))
    print(mean(index))
    print(max(index))}
    
sumstatdf <- function(index){
    width <- index$end - index$start
    print(min(width))
    print(median(width))
    print(mean(width))
    print(max(width))
    print(nrow(index))}
   
genrich <- fread("sox10ATAC_genrich.narrowPeaks")
dualMACs <- fread("sox10ATAC_MACS2_narrowPeaks.bed")

colnames(dualMACs) <- c("chr","start","end","name","score","strand","fold_enrichment","pval","qval","blockCount")
colnames(genrich) <- c("chr","start","end","name","score","strand","signalValue","pval","qval","peak")


genp <- dplyr::filter(genrich,pval<=5)
macsp <- dplyr::filter(dualMACs,pval<=5) 

reprod <- fread("sox10_Zv9FromFASTQNoncodingUnder10kReprodOnly.bed")
colnames(reprod) <- c("chr","start","end")

greprod <- GRanges(reprod)
ggen <- GRanges(genrich)
gmacs <- GRanges(dualMACs)
ggenp <- GRanges(genp)
gmacsp <- GRanges(macsp)

