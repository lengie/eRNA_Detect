### DLModel_FPFNExtraction.r
### Testing tensorflow - keras model by extracting false positive and false negatives and original genomic regions
###
### Purpose: isolate genomic locations of putative eRNAs or non-enhancer-associated bidirectional regions
### Input: 
### Output: bed files with genomic locations corresponding to transcriptional profiles flagged as false positives or negatives

## load model then load weights
model %>% load_model_weights_hdf5('[weights].h5')

## input data
putenh <- rbind(h0_enh,h6_enh,sox1_enh,sox2_enh)
noise <- rbind(h0_no,h6_no,sox1_no,sox2_no)

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


whichh0n <- which(samp_enh_fneg <= 17137) 
h0_fn_index <- samp_enh_fneg[whichh0n]
h0_false_neg <- h0_put[h0_fn_index,]
write.table(h0_false_neg,"h0_falseneg_ConvDensePoolDenseDense79Acc333.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

whichh6 <- which(samp_enh_fneg > 17137 & samp_enh_fneg <= (17137+4429))
h6_fn_index <- samp_enh_fneg[whichh6]- 17137
h6_false_neg <- h6_put[h6_fn_index,] 
write.table(h6_false_neg,"h6_falseneg_ConvDensePoolDenseDense79Acc333.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

whichsox1 <- which(samp_enh_fneg > (17137+4429) & samp_enh_fneg <= (17137+4429+82811))
sox1_fn_index <- samp_enh_fneg[whichsox1]- (17137+4429)
sox1_false_neg <- sox1_put[sox1_fn_index,] 
write.table(sox1_false_neg,"sox1_falseneg_ConvDensePoolDenseDense79Acc333.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

whichsox2 <- which(samp_enh_fneg > (17137+4429+82811) & samp_enh_fneg <= (17137+4429+82811+135249))
sox2_fn_index <- samp_enh_fneg[whichsox2]- (17137+4429+82811)
sox2_false_neg <- sox2_put[sox2_fn_index,] 
write.table(sox2_false_neg,"sox2_falseneg_ConvDensePoolDenseDense79Acc333.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')


whichh0 <- which(samp_enh_fpos <= 509598)
h0_fp_index <- samp_enh_fpos[whichh0]
h0_false_pos <- h0_noise[h0_fp_index,] 
write.table(h0_false_pos,"h0_falsepos_ConvDensePoolDenseDense79Acc333.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

whichh6 <- which(samp_enh_fpos > 509598 & samp_enh_fpos <= (509598+112698))
h6_fp_index <- samp_enh_fpos[whichh6]  - 509598
h6_false_pos <- h6_noise[h6_fp_index,] 
write.table(h6_false_pos,"h6_falsepos_ConvDensePoolDenseDense79Acc333.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

whichsox1 <- which(samp_enh_fpos > (509598+112698) & samp_enh_fpos <= (509598+112698))
sox1_fp_index <- samp_enh_fpos[whichsox1]- - (509598+112698)
sox1_false_pos <- sox1_noise[sox1_fp_index,] 
write.table(sox1_false_pos,"sox1_falsepos_ConvDensePoolDenseDense79Acc333.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

whichsox2 <- which(samp_enh_fpos > (509598+112698+452236))
sox2_fp_index <- samp_enh_fpos[whichsox2] - (509598+112698+452236)
sox2_false_pos <- testsox2_no[sox2_fp_index,] 
write.table(sox2_false_pos,"sox2_falsepos_ConvDensePoolDenseDense79Acc333.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')