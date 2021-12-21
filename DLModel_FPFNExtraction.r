### DLModel_FPFNExtraction.r
### Testing tensorflow - keras model by extracting false positive and false negatives and original genomic regions
###
### Purpose: isolate genomic locations of putative eRNAs or non-enhancer-associated bidirectional regions
### Input: 
### Output: bed files with genomic locations corresponding to transcriptional profiles flagged as false positives or negatives

## load model then load weights
model %>% load_model_weights_hdf5('[weights].h5')

## input data
putenh <- rbind(YOUR DATASETS)
noise <- rbind(YOUR DATASETS)

set.seed(333)
noise_pool_ind <- sample(dim(noise)[1], train_noise_no)
putenh_pool_ind <- sample(dim(putenh)[1],train_putenh_no)
train_noise <- noise[noise_pool_ind,]
train_enh <- putenh[putenh_pool_ind,]

# shuffle and put it all into one dataset
training <- rbind(train_enh,train_noise)
shuffle <- sample(nrow(training))
shuffled <- training[shuffle,]
reshape <- array_reshape(as.matrix(shuffled),list(length(shuffle),120,1))
#reshape_k <- k_reshape(as.matrix(shuffled),list(length(shuffle),120,1))

## Note: 0 is enhancer 1 is noise in this categorization
categories <- c(rep(0,times=nrow(train_enh)),rep(1,times=nrow(train_noise)))
categories <- categories[shuffle]
shuffled <- as.matrix(shuffled)
y_all_training <- to_categorical(categories)

## test data
evals <- model %>% evaluate(all_reshape, y_training, batch_size = 512)

# Get confusion matrix for predictions
#classes <- model %>% predict_classes(all_reshape, batch_size=512)
classes <- model %>% predict(reshape) %>% k_argmax()
ct <- table(round(as.matrix(classes)),y_all_training[,2])
cm <- as.matrix(ct)
cm


sidebyside <- cbind(as.matrix(classes),y_all_training[,2])
diff <- which(sidebyside[,1] != sidebyside[,2])
length(diff)
enh_false_neg <- which(sidebyside[,1] != sidebyside[,2] & sidebyside[,2]==0)
enh_false_pos <- which(sidebyside[,1] != sidebyside[,2] & sidebyside[,2]==1)
shuffle_enh_false_neg <- shuffle[enh_false_neg]
shuffle_enh_false_pos <- shuffle[enh_false_pos]
put_false_neg <- which(shuffle_enh_false_neg <= train_putenh_no)
noise_false_neg <- which(shuffle_enh_false_neg > train_putenh_no)

put_false_pos <- which(shuffle_enh_false_pos <= train_putenh_no)
noise_false_pos <- which(shuffle_enh_false_pos > train_putenh_no)
samp_enh_fneg <- putenh_pool_ind[shuffle_enh_false_neg]
samp_enh_fpos <- noise_pool_ind[shuffle_enh_false_pos - train_putenh_no]


splitAndSave <- function(posneg,start,end,putnoise){
  if(posneg=='pos'){
    whichind <- which(samp_enh_fpos > start & samp_enh_fpos <= end)
    fp_index <- samp_enh_fpos[whichind] - start
    falsepn <- putnoise[fp_index]
    print(nrow(falsepn))
    bedname <- paste(deparse(substitute(putnoise)),"_FalsePos_5layersBinary.bed",sep="")
  }else if(posneg=='neg'){
    whichind <- which(shuffle_enh_false_neg > start & shuffle_enh_false_neg <= end)
    fn_index <- shuffle_enh_false_neg[whichind] - start
    falsepn <- putnoise[fn_index]
    print(nrow(falsepn))
    bedname <- paste(deparse(substitute(putnoise)),"_FalseNeg_5layersBinary.bed",sep="")
  }else{print("Please set whether false positive or negative with 'pos' or 'neg'!")}
  bedname <- paste(deparse(substitute(putnoise)),"_5layersBinary.bed",sep="")
  print(bedname)
  write.table(falsepn,bedname,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
}

######### Binary categories, 1 is enhancer 0 is non-eRNA, training and evaluation data processed simultaneously ############
set.seed(123)
train_putenh_no <- floor(nrow(putenh)/2)
noise_pool_ind <- sample(dim(noise)[1], (train_putenh_no*2))
all_noise <- noise[noise_pool_ind,]

# shuffle and put it all into one dataset
all_data <- rbind(putenh,all_noise)
shuffle <- sample(nrow(all_data))
shuffled <- all_data[shuffle,]

# split into training vs eval
categories <- c(rep(1,times=nrow(putenh)),rep(0,times=nrow(all_noise)))
categories <- categories[shuffle]

training <- shuffled[1:(train_putenh_no*2),]
train_cat <- categories[1:(train_putenh_no*2)]
y_all_training <- to_categorical(train_cat)
reshape <- array_reshape(as.matrix(training),list(nrow(training),120,1))

# evaluation data
evaluation <- shuffled[(train_putenh_no*2 + 1):nrow(shuffled),]
reshape_eval <- array_reshape(as.matrix(evaluation),list(nrow(evaluation),120,1))
eval_cat <- categories[(train_putenh_no*2 + 1):nrow(shuffled)]
y_all_eval <- as.matrix(eval_cat)

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
splitAndSave <- function(posneg,start,end,putnoise,suffix){
  if(posneg=='pos'){
    whichind <- which(samp_enh_fpos > start & samp_enh_fpos <= end)
    fp_index <- samp_enh_fpos[whichind] - start
    falsepn <- putnoise[fp_index]
    print(nrow(falsepn))
  }else if(posneg=='neg'){
    # putative enhancers are in the first half of all_data and are not shuffled
    whichind <- which(shuffle_enh_false_neg > start & shuffle_enh_false_neg <= end)
    fn_index <- shuffle_enh_false_neg[whichind] - start
    falsepn <- putnoise[fn_index]
    print(nrow(falsepn))
  }else{print("Please set whether false positive or negative with 'pos' or 'neg'!")}
  bedname <- paste(deparse(substitute(putnoise)),suffix,sep="")
  print(bedname)
  write.table(falsepn,bedname,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
}
