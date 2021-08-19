### SequentialDLModel.r
### Initial code for basic tensor flow - keras model
###
### Purpose: separate ATAC-seq peak overlapping bidirectional regions from noise
### Input: transcription across 50bp bins spanning 3kb on both strands of genomic regions
### Output: predictions

library(data.table)
library(keras)
library(tensorflow)
library(dplyr)
options(scipen=999)
use_condaenv("/project/rohs_102/keras")

## data that has been combined together already
data <- fread("sox10_DLtestingBPM.csv")
y_all <- data$V121
x_all <- data[,1:120]
y_cat <- to_categorical(y_all,2)

## OR load raw datasets
# the way Galaxy had uploaded the bigwig files, minus is used first in createMatrix
bins <- c(paste(seq(-30,-1,by=1),"minus",sep=""),
            paste(seq(1,30,by=1),"minus",sep=""),
            paste(seq(-30,-1,by=1),"plus",sep=""),
            paste(seq(1,30,by=1),"plus",sep=""))
# each sample is a row and each feature/bin is a column so it looks like an image file that's been reformatted into a single line


## BPM-scaled datasets
BPMovfile <- fread("sox10_BidirOverlapATAC3kb50bpbinsBPM.tabular")
BPMnotov_file <- fread("sox10_BidirNotOverlapATAC3kb50bpbinsBPM.tabular")
BPMoverlaps <- as.matrix(BPMovfile[,1:120])
BPMnot_ov <- as.matrix(BPMnotov_file[,1:120])
colnames(BPMoverlaps) <- bins
colnames(BPMnot_ov) <- bins

# bactin data sets: both replicates
data <- fread("bactnuc869and70_bothBW3kb50bins.tabular")
# convert the NAs to 0 before converting to a matrix, so it is definitely a numeric matrix
data <- data %>% replace(is.na(.), 0)
x_bact <- as.matrix(data[,1:120])
colnames(x_bact) <- bins

## Pull out 1/8th of data for validation
set.seed(111)

# n and n2 are the number of each condition to put in the VALIDATION training data set
set_training_valid <- function(cond1,cond2,n,n2){  
    overlap_ind <- sample(dim(cond1)[1], n) 
    noov_ind <- sample(dim(cond2)[1], n2)
    
    x_training <- rbind( cond1[-overlap_ind,],
                        cond2[-noov_ind,] )
    x_valid <- rbind(cond1[overlap_ind,],
                     cond2[noov_ind,] )
    
    # categories should be a binary class matrix
    cond1left <- dim(cond1)[1]-n
    cond2left <- dim(cond2)[1]-n2

    y_training <- matrix(c(rep(1,times=cond1left),rep(0,times=cond2left),rep(0,times=cond1left),rep(1,times=cond2left)),
                         nrow=(cond1left+cond2left), ncol=2,
                         byrow=FALSE,
                         dimnames= list(c(1:(cond1left+cond2left)),c("putenh","noise")) )
    y_valid <- matrix(c(rep(1,times=n),rep(0,times=n2),rep(0,times=n),rep(1,times=n2)),
                    nrow=(n+n2), ncol=2,
                    byrow=FALSE,
                    dimnames= list(c(1:(n+n2)),c("putenh","noise")) )

    return(list(x_training,y_training, x_valid,y_valid))
}
## simplified version with no numbers
set_training <- function(cond1,cond2){     
    x_training <- rbind(cond1,cond2)

    # categories should be a binary class matrix
    len <- nrow(cond1)
    len2 <- nrow(cond2)
    y_training <- matrix(c(rep(c(1,0),times=len),rep(c(0,1),times=len2)),
                         nrow=(len+len2), ncol=2,
                         byrow=TRUE,
                         dimnames= list(c(1:(len+len2)),c("putenh","noise")) )
    return(list(x_training,y_training))
}


datasets <- set_training_valid(BPMoverlaps,BPMnot_ov,5860,62900)


## the actual DL model
model <- keras_model_sequential()
model %>%
    layer_dense(units= 512, activation="relu", input_shape=c(120)) %>%
    layer_dense(units= 512, activation="relu") %>%
    layer_dropout(rate=0.1) %>%
    layer_dense(units=2, activation="softmax")

summary(model)
model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
)

history <- model %>% fit(
    datasets[[1]], datasets[[2]],
    epochs=20, batch_size = 10000,
    validation_split = 0.2,
    verbose = 1
)
