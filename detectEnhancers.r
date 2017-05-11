### detectEnhancers.R
###
### Purpose: Take RNA-seq data and output locations and FPKMs of regions with bidirectional reads within a given region 
###
###
### Written by Liana Engie
### Last updated: April 2017
###
### detectEnhancers(chromosome, input_start, input_end,strand)
### Input: string chromosome number, 
### Output:

ranges <- findReadRegions(chromosome,input_start,input_end,strand,bed){
	interval <- IRanges(input_start,input_end)
	list <- subsetByOverlaps(bed,interval)
	
	start <- c()
	end <- c()
	j <- 1
	
	if(length(ranges(list))==0){
		print('No clusters of reads in range.')
	} else{
		while(j < length(ranges(list))){
			start <- c(start, start(list)[j])
			end <- c(end, end(list)[j])
			j <- j+1
		}
	}
	
	ranges <- IRanges(start,end)
}

removeZeroCounts(chromosome, ranges, strand="+", sense_counts){

	if(strand!="-"||strand!="+"){
		print('Strand should be "-" or "+"')
		break 
	}

	j <- 0
	newranges <- IRanges()
	i <- length(min_counts)
	if(sense_counts[i]==0){
		newranges <- append(newranges,ranges[i])
		j <- j+1
	}if(j==0){break}

#	if(strand=="-"){
#		anti <- "+"
#	}else{anti <- "-"}

	int <- GRanges(seqnames=chromosome,ranges=newranges,strand=anti)
	counts <- summarizeOverlaps(features=int,reads=hetsread1,singleEnd=FALSE,fragments=FALSE,inter.feature=FALSE)
	anti_counts <- assay(counts)

}

