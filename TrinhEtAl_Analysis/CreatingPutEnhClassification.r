remOverlap <- function(gr,toRem){
    compare <- findOverlaps(toRem,gr,ignore.strand=TRUE)
    keep <- gr[-subjectHits(compare),]
    return(keep)
}

sumstatgr <- function(index){
    print(min(width(index)))
    print(median(width(index)))
    print(mean(width(index)))
    print(max(width(index)))
    print(length(index))}

testov <- function(gr,minus_gr1,minus_gr2){
    print("findOverlaps rep1...")
    min <- remOverlap(gr,minus_gr1)
    print("summary stats and basepair coverage rep1")  
    sumstatgr(min)
    sum(width(min))
    print("findOverlaps rep2...")
    min2 <- remOverlap(gr,minus_gr2)
    print("summary stats and basepair coverage rep2")  
    sumstatgr(min2)
    sum(width(min2))
    print("Jaccard distance between them...")
    genomicCorr.jaccard(min,min2)

    print("setdiff rep1...")
    sd1 <- GenomicRanges::setdiff(gr,minus_gr1,ignore.strand=TRUE)
    print("summary stats and basepair coverage rep1 setdiff")  
    sumstatgr(sd1)
    sum(width(sd1))

    print("setdiff rep2...")
    sd2 <- GenomicRanges::setdiff(gr,minus_gr2,ignore.strand=TRUE)
    print("summary stats and basepair coverage rep2 setdiff")  
    sumstatgr(sd2)
    sum(width(sd2))
    
    print("Jaccard distance between setdiffs...")
    genomicCorr.jaccard(sd2,sd2)
}
