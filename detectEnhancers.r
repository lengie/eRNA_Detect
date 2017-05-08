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
	
	#search the bed file for a region that includes the input_start, if any
	if(bed$chr[i]==chromosome && bed$strand[i]=strand && bed$start[i] <= input_start <= bed$end[i]){ #just the idea of it so far  
		start <- bed$start[i]
		end <- bed$end[i]
	}else{
		end <- input_start 
	}
	
	while(end < input_end){
	
	
	
		#Look for the bed file line that includes start, if any
		#Look for next start, and where it ends
		#If end < input.end, then keep going #but have to change end if there aren't any reads in the region
		
	
	start <- c()
	end <- c()
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

