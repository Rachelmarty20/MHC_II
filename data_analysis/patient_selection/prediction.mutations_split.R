library(parallel)
library(lme4)
library(mgcv)
library(xtable)
library(pROC)

PATH_TO_DATA = '/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/'

args = commandArgs(trailingOnly=TRUE)
mutation_threshold = as.integer(args[1])

#Format data
tissue <- read.csv(paste(PATH_TO_DATA, 'patient_tissues.csv', sep=""),header=TRUE)
mut <- read.csv(paste(PATH_TO_DATA, 'combined_classes/patient_mutations.csv', sep=""),header=TRUE)
aff1 <- read.csv(paste(PATH_TO_DATA, 'combined_classes/patient_affinities.class_i.csv', sep=""),header=TRUE)
aff2 <- read.csv(paste(PATH_TO_DATA, 'combined_classes/patient_affinities.class_ii.csv', sep=""),header=TRUE)
patient <- as.character(mut[,1])
mut <- as.matrix(mut[,-1])
aff1 <- as.matrix(aff1[,-1])
aff2 <- as.matrix(aff2[,-1])
rownames(mut) <- rownames(aff1) <- rownames(aff2) <- patient

# probably need to update
y= as.vector(mut); x= as.vector(aff1); z= as.vector(aff2)
gene= rep(colnames(mut),each=nrow(mut))
pat= rep(rownames(mut),ncol(mut))
nmut= colSums(mut)
sel= gene %in% names(nmut[nmut>=mutation_threshold])

# to select a smaller data subset
#sampled_pats = head(unique(pat), 500)
sampled_pats = unique(pat)
sel_pat = pat %in% sampled_pats

df = data.frame(y[sel&sel_pat], x[sel&sel_pat], z[sel&sel_pat], pat[sel&sel_pat])
colnames(df)<-c('y', 'x', 'z', 'pat')


# both MHC-I and MHC-II
all_labels=NULL
all_predictions=NULL
for (i in 1:3)
{
    print(i)
    # sample indices
    sample_rows = sample(nrow(df), round(nrow(df)/10))
    # to test the model
    DataC1=df[sample_rows, ]
    # to train the model
    DataCV=df[-sample_rows, ]
    # train the model
    M1 <- glmer(y ~ log(z) + log(x) + (1|pat), data=DataCV, family='binomial')
    # predict mutation probabilities
    P1=predict(M1, DataC1)
    names(P1)=NULL
    all_predictions= c(all_predictions, P1)
    all_labels = c(all_labels, DataC1$y)
}
# put in dataframe and transform
results_df = data.frame(all_labels, all_predictions)
colnames(results_df)<-c('label', 'predicted')
# probabilities
results_df$predicted_prob<-exp(results_df$predicted)
# labels as factors
results_df$label_fact <- factor(results_df$label)
# create the roc object
roc_obj <- roc(results_df$label_fact, results_df$predicted_prob)


# only MHC-I
all_labels=NULL
all_predictions=NULL
for (i in 1:3)
{
    print(i)
    # sample indices
    sample_rows = sample(nrow(df), round(nrow(df)/10))
    # to test the model
    DataC1=df[sample_rows, ]
    # to train the model
    DataCV=df[-sample_rows, ]
    # train the model
    M1 <- glmer(y ~ log(x) + (1|pat), data=DataCV, family='binomial')
    # predict mutation probabilities
    P1=predict(M1, DataC1)
    names(P1)=NULL
    all_predictions= c(all_predictions, P1)
    all_labels = c(all_labels, DataC1$y)
}
# put in dataframe and transform
results_dfI = data.frame(all_labels, all_predictions)
colnames(results_dfI)<-c('label', 'predicted')
# probabilities
results_dfI$predicted_prob<-exp(results_dfI$predicted)
# labels as factors
results_dfI$label_fact <- factor(results_dfI$label)
# create the roc object
roc_objI <- roc(results_dfI$label_fact, results_dfI$predicted_prob)


# only MHC-II
all_labels=NULL
all_predictions=NULL
for (i in 1:3)
{
    print(i)
    # sample indices
    sample_rows = sample(nrow(df), round(nrow(df)/10))
    # to test the model
    DataC1=df[sample_rows, ]
    # to train the model
    DataCV=df[-sample_rows, ]
    # train the model
    M1 <- glmer(y ~ log(z) + (1|pat), data=DataCV, family='binomial')
    # predict mutation probabilities
    P1=predict(M1, DataC1)
    names(P1)=NULL
    all_predictions= c(all_predictions, P1)
    all_labels = c(all_labels, DataC1$y)
}
# put in dataframe and transform
results_dfII = data.frame(all_labels, all_predictions)
colnames(results_dfII)<-c('label', 'predicted')
# probabilities
results_dfII$predicted_prob<-exp(results_dfII$predicted)
# labels as factors
results_dfII$label_fact <- factor(results_dfII$label)
# create the roc object
roc_objII <- roc(results_dfII$label_fact, results_dfII$predicted_prob)

# output results
auc_summary <-vector("list", 3)
auc_summary[[1]] <- c(ci(roc_obj))
auc_summary[[2]] <- c(ci(roc_objI))
auc_summary[[3]] <- c(ci(roc_objII))
auc_df <- data.frame(auc_summary)
colnames(auc_df) <- c('Both', 'Only_I', 'Only_II')
rownames(auc_df) <- c('low_CI', 'AUC', 'high_CI')
write.table(auc_df, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/predictions.mutations_split.", args[1], ".txt", sep=''))

# Plot the ROCs
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/predictions/ROC.mutations_split.threshold_', args[1], '.pdf', sep=''))
plot(roc_obj, col='red')
plot(roc_objI, col='darkgreen', add=TRUE)
plot(roc_objII, col='blue', add=TRUE)
legend("bottomright",
  legend = c("MHC-I and MHC-II", "MHC-I", "MHC-II"),
      col=c('red','darkgreen', 'blue'),
      lwd=c(2.5,2.5))
 dev.off()