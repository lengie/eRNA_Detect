# use a conda environment with tensorflow and keras loaded
use_condaenv("keras")
library(data.table)
library(keras)
library(tensorflow)
library(dplyr)
#library(kerasR)
options(scipen=999)

bins <- c(paste("plus",seq(-30,-1,by=1),sep=""),
            paste("plus",seq(1,30,by=1),sep=""),
            paste("minus",seq(-30,-1,by=1),sep=""),
            paste("minus",seq(1,30,by=1),sep=""))

# if there were no reads, I had deepTools output an NA so that would appear white in the heatmaps. Convert to 0 for the program
c1_1tab <- fread("Prescott2015/chimp_1_1_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
c1_2tab <- fread("Prescott2015/chimp_1_2_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
c2_1tab <- fread("Prescott2015/chimp_2_1_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
c2_2tab <- fread("Prescott2015/chimp_2_2_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h1_1tab <- fread("Prescott2015/human_1_1_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h1_2tab <- fread("Prescott2015/human_1_2_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h2_1tab <- fread("Prescott2015/human_2_1_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h2_2tab <- fread("Prescott2015/human_2_2_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h3_1tab <- fread("Prescott2015/human_3_1_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h3_2tab <- fread("Prescott2015/human_3_2_enh_overlaps_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)

sox1_enh <- fread("DLTesting/sox10_871_bidirATACnoribo_pt1rpb_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
sox2_enh <- fread("DLTesting/sox10_873_bidirATACnoribo_pt1rpb_3kb_50bpTEST.tabular",header=FALSE) %>% replace(is.na(.), 0)

sox1_no <- fread("DLTesting/sox10_871_bidirUnder500bpMerged800bpNoATACOverlap.tabular",header=FALSE) %>% replace(is.na(.), 0)
sox2_no <- fread("DLTesting/sox10_873_bidirUnder500bpMerged800bpNoATACOverlap.tabular",header=FALSE) %>% replace(is.na(.), 0)
c1_1no <- fread("Prescott2015/chimp_1_1_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
c1_2no <- fread("Prescott2015/chimp_1_2_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
c2_1no <- fread("Prescott2015/chimp_2_1_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
c2_2no <- fread("Prescott2015/chimp_2_2_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h1_1no <- fread("Prescott2015/human_1_1_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h1_2no <- fread("Prescott2015/human_1_2_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h2_1no <- fread("Prescott2015/human_2_1_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h2_2no <- fread("Prescott2015/human_2_2_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h3_1no <- fread("Prescott2015/human_3_1_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)
h3_2no <- fread("Prescott2015/human_3_2_BidirNoEnh_3kb_50bp.tabular",header=FALSE) %>% replace(is.na(.), 0)

colnames(sox1_enh) <- bins
colnames(sox1_no) <- bins	
colnames(sox2_enh) <- bins
colnames(sox2_no) <- bins
colnames(c1_1tab) <- bins
colnames(c1_2tab) <- bins
colnames(c1_1no) <- bins
colnames(c1_2no) <- bins
colnames(c2_1tab) <- bins
colnames(c2_2tab) <- bins
colnames(c2_1no) <- bins
colnames(c2_2no) <- bins
colnames(h1_1tab) <- bins
colnames(h1_2tab) <- bins
colnames(h1_1no) <- bins
colnames(h1_2no) <- bins
colnames(h2_1tab) <- bins
colnames(h2_2tab) <- bins
colnames(h2_1no) <- bins
colnames(h2_2no) <- bins
colnames(h3_1tab) <- bins
colnames(h3_2tab) <- bins
colnames(h3_1no) <- bins
colnames(h3_2no) <- bins

putenh <- rbind(h1_1tab,h1_2tab,h2_1tab,h2_2tab,h3_1tab,h3_2tab,c1_1tab,c1_2tab,c2_1tab,c2_2tab,sox1_enh,sox2_enh)
noise <- rbind(h1_1no,h1_2no,h2_1no,h2_2no,h3_1no,h3_2no,c1_1no,c1_2no,c2_1no,c2_2no,sox1_no,sox2_no)

# create the training, validation, and evaluation pools
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


## FIRST MODEL
model <- keras_model_sequential()
model %>%
         layer_conv_1d(filters = 16, kernel_size = 32, activation = "relu", input_shape = c(120,1)) %>%
         layer_dense(units=8, activation="relu") %>%
         layer_max_pooling_1d(pool_size=8) %>%
         layer_dense(units=24, activation="relu") %>%
         layer_dense(units=8, activation="relu") %>%
         layer_flatten() %>%
         layer_dense(units=2, activation="softmax")
summary(model)

model %>% compile(loss = 'binary_crossentropy',
            optimizer = 'adam',
            metrics = c('accuracy','AUC')
        )

history <- model %>% fit(
            reshape,y_all_training,
            epochs=200, batch_size = 512,
            validation_split=0.2,
            verbose = 1,shuffle=TRUE)
 
png("PrescottTrinh_ConvDensePoolDenseDenseBinary200ep_plot.png")
	plot(history)
dev.off()

save_model_hdf5(model,filepath="PrescottTrinh_Trained_5layersBinary200epochs.h5")

likelihood <- model %>% predict(reshape)
binary <- likelihood %>% k_argmax()
table(as.matrix(binary),y_all_training[,2])
classes <- model %>% predict(reshape_eval) %>% k_argmax()
table(round(as.matrix(classes)),y_all_eval)
head(likelihood)


# SECOND MODEL
model2 <- keras_model_sequential()
model2 %>%
         layer_conv_1d(filters = 32, kernel_size = 8, activation = "relu", input_shape = c(120,1)) %>%
         layer_batch_normalization() %>%
         layer_conv_1d(filters = 32, kernel_size = 8, activation = "relu") %>%
         layer_batch_normalization() %>%
         layer_max_pooling_1d(pool_size=2) %>%
         layer_conv_1d(filters = 16, kernel_size = 3, activation = "relu") %>%
         layer_batch_normalization() %>%
         layer_conv_1d(filters = 16, kernel_size = 3, activation = "relu") %>%
         layer_batch_normalization() %>%
         layer_max_pooling_1d(pool_size=2) %>%
         layer_dense(units= 16, activation="relu") %>%
         layer_dense(units=8, activation="relu") %>%
	 layer_flatten() %>%
         layer_dense(units=2, activation="softmax")
summary(model2)

model2 %>% compile(loss = 'binary_crossentropy',
            optimizer = 'adam',
            metrics = c('accuracy','AUC')
        )
history <- model2 %>% fit(
            reshape,y_all_training,
            epochs=200, batch_size = 512,
            validation_split=0.2,
            verbose = 1,shuffle=TRUE)


png("PrescottTrinh_ConvBNConvBNPoolConvBNConvBNPoolDenseDenseBinary200ep_plot.png")
	plot(history)
dev.off()

likelihood2 <- model2 %>% predict(reshape)
binary2 <- likelihood2 %>% k_argmax()
table(as.matrix(binary2),y_all_training[,2])
classes2 <- model2 %>% predict(reshape_eval) %>% k_argmax()
table(round(as.matrix(classes2)),y_all_eval)
head(likelihood2)

save_model_hdf5(model,filepath="PrescottTrinh_Trained_12layersBinary200epochs.h5")
