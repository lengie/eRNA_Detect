### DLFalsePosNegExtract.r
###
### Written by Liana Engie
### Nov 18 2021
###
### Written to extract the false positives and negatives extracted from training, testing, evaluation of a tensorflow keras deep learning model in R

## Data had all been processed in one data frame then split into testing and eval
# isolate which were mislabeled. Training was top half, eval was bottom half
sidebyside <- rbind(cbind(as.matrix(binary),y_all_training[,2]),
                    cbind(as.matrix(classes),y_all_eval) )

diff <- which(sidebyside[,1] != sidebyside[,2])
length(diff)  # sanity check

enh_false_neg <- which(sidebyside[,1] != sidebyside[,2] & sidebyside[,2]==1) # enhancers are labeled as 1
enh_false_pos <- which(sidebyside[,1] != sidebyside[,2] & sidebyside[,2]==0)

# location in shuffled all_data
shuffle_enh_false_neg <- shuffle[enh_false_neg]
shuffle_enh_false_pos <- shuffle[enh_false_pos]
samp_enh_fpos <- noise_pool_ind[shuffle_enh_false_pos - train_putenh_no]

# sanity check. All the false negatives should be in first 239626 rows
min(shuffle_enh_false_neg)
max(shuffle_enh_false_neg)
min(shuffle_enh_false_pos)
max(shuffle_enh_false_pos)

length(shuffle_enh_false_neg)
length(enh_false_pos)

# we could just use the transcriptional information here to make new heatmaps directly instead of backtracking to the original regions then making new matrices and then heatmaps, but deeptools takes the compressed data and I'm not sure how to generate them from text matrices. Anyway this let's me remap the genomic locations to other genome builds if I want to, which I do need to for the mouse data.
## load genomic locations of putative enhancers and not-eRNA-bidirectional-regions
loadbd <- function(file){
    bd <- fread(file)
    colnames(bd) <- c("chr","start","end","ID","score","strand")
    return(bd)
}

PUT <- loadbd(FILENAME)
NOISE <- loadbd(FILENAME)

## so we have their locations in all_data. We need to split putenh and all_noise
# split putenh into composite parts
splitAndSave <- function(posneg,start,end,putnoise){
  if(posneg=='pos'){
    whichind <- which(samp_enh_fpos > start & samp_enh_fpos <= end)
    fp_index <- samp_enh_fpos[whichind] - start
    falsepn <- putnoise[fp_index]
    print(nrow(falsepn))
    bedname <- paste(deparse(substitute(putnoise)),"_FalsePos_5layersBinary.bed",sep="")
  }else if(posneg=='neg'){
    # putative enhancers are in the first half of all_data and are not shuffled
    whichind <- which(shuffle_enh_false_neg > start & shuffle_enh_false_neg <= end)
    fn_index <- shuffle_enh_false_neg[whichind] - start
    falsepn <- putnoise[fn_index]
    print(nrow(falsepn))
    bedname <- paste(deparse(substitute(putnoise)),"_FalseNeg_5layersBinary.bed",sep="") #NOTE I CHANGED THE FILE NAME BASED ON MODEL
  }else{print("Please set whether false positive or negative with 'pos' or 'neg'!")}
  bedname <- paste(deparse(substitute(putnoise)),"_5layersBinary.bed",sep="")
  print(bedname)
  write.table(falsepn,bedname,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
}
