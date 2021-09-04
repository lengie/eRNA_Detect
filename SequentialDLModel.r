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

# for evaluation later
multi_class_rates <- function(confusion_matrix) {
    true_positives  <- diag(confusion_matrix)
    false_positives <- colSums(confusion_matrix) - true_positives
    false_negatives <- rowSums(confusion_matrix) - true_positives
    true_negatives  <- sum(confusion_matrix) - true_positives -
        false_positives - false_negatives
    return(data.frame(true_positives, false_positives, true_negatives,
                      false_negatives, row.names = names(true_positives)))
}

## column names to keep track of everything
# plus strand came first, then minus in createMatrix
bins <- c(paste(seq(-30,-1,by=1),"plus",sep=""),
            paste(seq(1,30,by=1),"plus",sep=""),
            paste(seq(-30,-1,by=1),"minus",sep=""),
            paste(seq(1,30,by=1),"minus",sep=""))
# each sample is a row and each feature/bin is a column so it looks like an image file that's been reformatted into a single line

			
Hovfile <- fread("sox10_871873_HisatATACOverlaps3kb50bp.tabular")
Hnotov_file <- fread("sox10_871873_HisatATACNoOverlaps3kb50bp.tabular")
# overlaps file has first replicate then second
Hov <- as.matrix(Hovfile[,121:240])
Hnot_ov <- as.matrix(Hnotov_file[,121:240])
colnames(Hov) <- bins
colnames(Hnot_ov) <- bins

bnot_file <- fread("bactin_869-70_HisatATACNoOverlaps3kb50bp.tabular")
# convert the NAs to 0 before converting to a matrix, so it is definitely a numeric matrix
bnot_file <- bnot_file %>% replace(is.na(.), 0)
bnot_ov <- as.matrix(bnot_file[,1:120])
colnames(bnot_ov) <- bins
bovfile <- fread("bactin_869-70_HisatATACOverlaps3kb50bp.tabular")
bovfile <- bovfile %>% replace(is.na(.),0)
bov <- as.matrix(bovfile[,1:120])

# because I have already separated ATAC overlapping from non-overlapping, no need to get a matrix of 1s and 0s

#second replicate for testing
nuc1Overlaps <- as.matrix(Hovfile[,1:120]) #the first rep
nuc1Not_ov <- as.matrix(Hnotov_file[,1:120])
colnames(nuc1Overlaps) <- bins
colnames(nuc1Not_ov) <- bins
bov_train <- as.matrix(bovfile[,121:240])
bno_train <- as.matrix(bnot_file[,121:240])
colnames(bov_train) <- bins
colnames(bno_train) <- bins
x_test <- rbind(nuc1Overlaps,nuc1Not_ov,bov_train,bno_train)
y_test <- as.matrix( c(rep(1,times=nrow(nuc1Overlaps)),rep(0,times=nrow(nuc1Not_ov)),
                    rep(1,times=nrow(bov_train)),rep(0,times=nrow(bno_train)) ))
y_test_cat <- to_categorical(y_test,2)
y_test_cat <- y_test_cat[,c(2,1)]

putenh <- rbind(Hov,bov)
noise <- rbind(Hnot_ov,bnot_ov)

putenhno <- floor(nrow(putenh)/8)
noiseno <- floor(nrow(noise)/8)

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
            layer_dense(units=1, activation="softmax")

        print(summary(model))

        model %>% compile(
            loss = 'binary_crossentropy',
            optimizer = 'adam',
            metrics = c('binary_accuracy','accuracy')
        )

        history <- model %>% fit(
            as.matrix(datasets[[1]]), datasets[[2]],
            epochs=ep[i], batch_size = 100,
            validation_data=list(as.matrix(datasets[[3]]),datasets[[4]]),
            verbose = 1)
                
        # save the model weights
        modname <- paste("sox10bactinTrained_4layers",seeds[j],"epoch",ep[i],".h5",sep="")
        save_model_hdf5(model,modname)
                
        # evaluate the model
        evals <- model %>% evaluate(x_test, y_test_cat, batch_size = 100)
        accuracy = evals[2][[1]]* 100
        print(accuracy)
                        
        # Get confusion matrix for predictions
        classes <- model %>% predict_classes(x_test, batch_size=100)
        ct <- table(round(classes),y_test)
        cm <- as.matrix(ct)
        print(cm)
        print(multi_class_rates(cm))
    }
}
