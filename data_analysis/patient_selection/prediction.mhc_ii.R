library(parallel)
library(lme4)
library(mgcv)
library(xtable)
library(pROC)

PATH_TO_DATA = '/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/'

args = commandArgs(trailingOnly=TRUE)
mutation_threshold = as.integer(args[1])
iterations = as.integer(args[2])
model = as.integer(args[3])

#Format data
tissue <- read.csv(paste(PATH_TO_DATA, 'patient_tissues.csv', sep=""),header=TRUE)
mut <- read.csv(paste(PATH_TO_DATA, 'patient_mutations.cancer.all.csv', sep=""),header=TRUE)
aff2 <- read.csv(paste(PATH_TO_DATA, 'patient_affinities.cancer.all.csv', sep=""),header=TRUE)
patient <- as.character(mut[,1])
mut <- as.matrix(mut[,-1])
aff2 <- as.matrix(aff2[,-1])
rownames(mut) <- rownames(aff2) <- patient

# probably need to update
y= as.vector(mut); z= as.vector(aff2)
gene= rep(colnames(mut),each=nrow(mut))
pat= rep(rownames(mut),ncol(mut))
nmut= colSums(mut)
sel= gene %in% names(nmut[nmut>=5])

# to select a smaller data subset
sampled_pats = head(unique(pat), 100)
#sampled_pats = unique(pat)
sel_pat = pat %in% sampled_pats

df = data.frame(y[sel&sel_pat], z[sel&sel_pat], pat[sel&sel_pat])
colnames(df)<-c('y', 'z', 'pat')

# linear - PHBR
if (model == 0){
    # only MHC-II
    all_labels=NULL
    all_predictions=NULL
    for (i in 1:iterations)
    {
        print(i)
        # sample indices
        sample_rows = sample(nrow(df), round(nrow(df)/10))
        # to test the model
        DataC1=df[sample_rows, ]
        # to train the model
        DataCV=df[-sample_rows, ]
        # train the model
        M1 <- glmer(y ~ z + (1|pat), data=DataCV, family='binomial')
        # predict mutation probabilities
        P1=predict(M1, DataC1)
        names(P1)=NULL
        all_predictions= c(all_predictions, P1)
        all_labels = c(all_labels, DataC1$y)
    }

}

# linear - log(PHBR)
if (model == 1){
    all_labels=NULL
    all_predictions=NULL
    for (i in 1:iterations)
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
}

# smoothed - PHBR
if (model == 2){
    all_labels=NULL
    all_predictions=NULL
    for (i in 1:iterations)
    {
        print(i)
        # sample indices
        sample_rows = sample(nrow(df), round(nrow(df)/10))
        # to test the model
        DataC1=df[sample_rows, ]
        # to train the model
        DataCV=df[-sample_rows, ]
        # train the model
        gam1= gam(y ~ s(z), data=DataCV, family='binomial')
        # predict mutation probabilities
        P1=predict(gam1, DataC1)
        names(P1)=NULL
        all_predictions= c(all_predictions, P1)
        all_labels = c(all_labels, DataC1$y)
    }
}


# smoothed - log(PHBR)
if (model == 3){
    all_labels=NULL
    all_predictions=NULL
    for (i in 1:iterations)
    {
        print(i)
        # sample indices
        sample_rows = sample(nrow(df), round(nrow(df)/10))
        # to test the model
        DataC1=df[sample_rows, ]
        # to train the model
        DataCV=df[-sample_rows, ]
        # train the model
        gam1= gam(y ~ s(log(z)), data=DataCV, family='binomial')
        # predict mutation probabilities
        P1=predict(gam1, DataC1)
        names(P1)=NULL
        all_predictions= c(all_predictions, P1)
        all_labels = c(all_labels, DataC1$y)
    }
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
auc_summary <-vector("list", 1)
auc_summary[[1]] <- c(ci(roc_objII))
auc_df <- data.frame(auc_summary)
colnames(auc_df) <- c('Only_II')
rownames(auc_df) <- c('low_CI', 'AUC', 'high_CI')
write.table(auc_df, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/predictions.mhc_ii.model_", args[3], ".", args[1], '.', args[2], ".txt", sep=''))

# Plot the ROCs
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/predictions/ROC.mhc_ii.model_', args[3], '.threshold_', args[1], '.', args[2], '.pdf', sep=''))
plot(roc_obj, col='red')
plot(roc_objI, col='darkgreen', add=TRUE)
plot(roc_objII, col='blue', add=TRUE)
legend("bottomright",
  legend = c("MHC-I and MHC-II", "MHC-I", "MHC-II"),
      col=c('red','darkgreen', 'blue'),
      lwd=c(2.5,2.5))
 dev.off()



# only MHC-II
all_labels=NULL
all_predictions=NULL
for (i in 1:iterations)
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
