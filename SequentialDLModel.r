### SequentialDLModel.r
### Initial code for basic tensor flow - keras model
###
### Purpose: separate ATAC-seq peak overlapping bidirectional regions from noise
### Input: transcription across 50bp bins spanning 3kb on both strands of genomic regions
### Output: predictions

## loading required packages
library(data.table)
library(keras)
library(tensorflow)
library(dplyr)
options(scipen=999)
# use a conda environment with tensorflow and keras loaded
use_condaenv("/project/rohs_102/keras")

## column names to keep track of everything
# the way Galaxy had uploaded the bigwig files, minus is used first in createMatrix
bins <- c(paste(seq(-30,-1,by=1),"minus",sep=""),
            paste(seq(1,30,by=1),"minus",sep=""),
            paste(seq(-30,-1,by=1),"plus",sep=""),
            paste(seq(1,30,by=1),"plus",sep=""))
# each sample is a row and each feature/bin is a column so it looks like an image file that's been reformatted into a single line

## data that has been combined together already. Transcription per bin in columns and last column is category
# category is 0,1 where 0 is noise and 1 is a putative enhancer
data <- fread("sox10_DLtestingBPM.csv")
y_all <- data$V121
x_all <- data[,1:120]
colnames(x_all) <- bins
y_cat <- to_categorical(y_all,2)

## OR load raw datasets
## BPM-scaled datasets
BPMovfile <- fread("sox10_BidirOverlapATAC3kb50bpbinsBPM.tabular")
BPMnotov_file <- fread("sox10_BidirNotOverlapATAC3kb50bpbinsBPM.tabular")
BPMoverlaps <- as.matrix(BPMovfile[,1:120])
BPMnot_ov <- as.matrix(BPMnotov_file[,1:120])
colnames(BPMoverlaps) <- bins
colnames(BPMnot_ov) <- bins

# bactin data sets: both replicates. Transcription per bin in columns and last column is category
data <- fread("bact_DLtestingBPM.csv")
x_bact <- data[,c(1:120,241)]
bact_ov <- dplyr::filter(x_bact,V241==1)
bact_no <- dplyr::filter(x_bact,V241==0)

# removing the column with categories
bact_ov <- bact_ov[,1:120]
colnames(bact_ov) <- bins
bact_no <- bact_no[,1:120]
colnames(bact_no) <- bins

# two category data sets
putenh <- rbind(BPMoverlaps,bact_ov)
noise <- rbind(BPMnot_ov,bact_no)


# splitting the datasets for training & testing
set_training_valid <- function(cond1,cond2,n,n2){  
       overlap_ind <- sample(dim(cond1)[1], n) 
       noov_ind <- sample(dim(cond2)[1], n2)
    
       # rbind the putenh large fraction and noise large fraction
       x_training <- rbind(as.matrix(cond1[-overlap_ind,]),
                           as.matrix(cond2[-noov_ind,] ))
       # rbind the smaller fractions, overlap then not-overlap below
       x_valid <- rbind(as.matrix(cond1[overlap_ind,]),
                        as.matrix(cond2[noov_ind,] ))
 
       # categories should be a binary class matrix. c("noise","putenh") c(0,1)
       cond1left <- dim(cond1)[1]-n
       cond2left <- dim(cond2)[1]-n2
        
       # putative enhancers are first and get number 1
       y_training <- matrix(c(rep(c(1,0),times=cond1left),rep(c(0,1),times=cond2left)),
                            nrow=(cond1left+cond2left), ncol=2,
                            byrow=TRUE,
                            dimnames= list(c(1:(cond1left+cond2left)),c("putenh","noise")) )
       y_valid <- matrix(c(rep(c(1,0),times=n),rep(c(0,1),times=n2)),
                        nrow=(n+n2), ncol=2,
                        byrow=TRUE,
                        dimnames= list(c(1:(n+n2)),c("putenh","noise")) )
       return(list(x_training,y_training, x_valid,y_valid))
}

putenhno <- floor(nrow(BPMoverlaps)/8)
noiseno <- floor(nrow(BPMnot_ov)/8)
datasets <- set_training_valid(BPMoverlaps,BPMnot_ov,putenhno,noiseno)


## testing (predicting) data
# data from other replicates for testing. Minus is still before plus
soxnuc1Ov_file <- <- fread("sox10_BidirOverlapATAC3kb50bpbinsNuc1BPM.tabular")
soxnuc1NoOv_file <- fread("sox10_BidirNotOverlapATAC3kb50bpbinsNuc1BPM.tabular")
nuc1Overlaps <- as.matrix(soxnuc1Ov_file[,1:120])
nuc1Not_ov <- as.matrix(soxnuc1NoOv_file[,1:120])
colnames(nuc1Overlaps) <- bins
colnames(nuc1Not_ov) <- bins
x_bact_test <- data[,121:240]
colnames(x_bact_test) <- bins

x_test <- rbind(nuc1Overlaps,nuc1Not_ov,x_bact_test)
y_test_cat <- c(rep(0,times=nrow(BPMoverlaps)),rep(1,times=nrow(BPMnot_ov)),data$V241) 
y_test <- to_categorical(y_test_cat,2)



## setting up combinations of epoch numbers and data splits to examine the model
ep <- c(5,10,20) 
seeds <- c(111,123,222,333,233)


for(i in 1:length(ep)){
    for(j in 1:length(seeds)){
        ## Pull out 1/8th of data for validation
        set.seed(seeds[i])
        print(paste("seed is ",seeds[j]," and number of epochs will be ",ep[i],sep=""))

        # n and n2 are the number of each condition to put in the VALIDATION training data set
    
        #ATAC-overlap goes first, so putative enhancers = 0, noise = 1
        datasets <- set_training_valid(putenh,noise,putenhno,noiseno)
    
        model <- keras_model_sequential()
        model %>%
            layer_dense(units= 512, activation="relu", input_shape=c(120)) %>%
            layer_dense(units= 512, activation="relu") %>%
            layer_dropout(rate=0.1) %>%
            layer_dense(units=2, activation="softmax")

        print(summary(model))

        model %>% compile(
            loss = 'categorical_crossentropy',
            optimizer = 'adam',
            metrics = c('accuracy')
        )

        history <- model %>% fit(
            datasets[[1]], datasets[[2]],
            epochs=ep[i], batch_size = 10000,
            validation_data=list(datasets[[3]],datasets[[4]]),
            verbose = 1)

        # Get confusion matrix for predictions
        classes <- model %>% predict_classes(x_test, batch_size=100)
        ct <- table(round(classes),y_test)
        cm <- as.matrix(cttest)
        print(cm)
        print(multi_class_rates(cm)

        # evaluate the model
        evals <- model %>% evaluate(x_test, t_test, batch_size = 100)
        accuracy = evals[2][[1]]* 100
        
    }
}
