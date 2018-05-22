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
mut <- read.csv(paste(PATH_TO_DATA, 'patient_mutations.cancer.TCGA.inclusive.mut.csv', sep=""),header=TRUE)
aff1 <- read.csv(paste(PATH_TO_DATA, 'patient_affinities.cancer.TCGA.conservative.mut.ClassI.csv', sep=""),header=TRUE)
aff2 <- read.csv(paste(PATH_TO_DATA, 'patient_affinities.cancer.TCGA.conservative.mut.ClassII.csv', sep=""),header=TRUE)

#mut <- read.csv(paste(PATH_TO_DATA, 'combined_classes/patient_mutations.csv', sep=""),header=TRUE)
#aff1 <- read.csv(paste(PATH_TO_DATA, 'combined_classes/patient_affinities.class_i.csv', sep=""),header=TRUE)
#aff2 <- read.csv(paste(PATH_TO_DATA, 'combined_classes/patient_affinities.class_ii.csv', sep=""),header=TRUE)
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
sel= gene %in% names(nmut[nmut>mutation_threshold])

# to select a smaller data subset
#sampled_pats = head(unique(pat), 500)
sampled_pats = unique(pat)
sel_pat = pat %in% sampled_pats

df = data.frame(y[sel&sel_pat], x[sel&sel_pat], z[sel&sel_pat], pat[sel&sel_pat])
colnames(df)<-c('y', 'x', 'z', 'pat')
# randomize order for CV
df <- df[sample(1:nrow(df)), ]
split = round(dim(df)[1] / 10)

#  s(z, x)
all_labels=NULL
all_predictions=NULL
# linear - PHBR
if (model == 0){
    # only MHC-II
    all_labels=NULL
    all_predictions=NULL
    for (i in 1:10)
    {
        # sample indices
        sample_rows = seq((i-1)*split+1, (i)*split)
        # to test the model
        DataC1=df[sample_rows, ]
        # to train the model
        DataCV=df[-sample_rows, ]
        # train the model
        M1 <- glm(y ~ z + x, data=DataCV, family='binomial')
        # predict mutation probabilities
        P1=predict(M1, DataC1)
        names(P1)=NULL
        all_predictions= c(all_predictions, P1)
        all_labels = c(all_labels, DataC1$y)
    }

}

# s(log(z), x)
if (model == 1){
    all_labels=NULL
    all_predictions=NULL
    for (i in 1:10)
    {
        # sample indices
        sample_rows = seq((i-1)*split+1, (i)*split)
        # to test the model
        DataC1=df[sample_rows, ]
        # to train the model
        DataCV=df[-sample_rows, ]
        # train the model
        M1 <- glm(y ~ log(z) + log(x), data=DataCV, family='binomial')
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
    for (i in 1:10)
    {
        # sample indices
        sample_rows = seq((i-1)*split+1, (i)*split)
        # to test the model
        DataC1=df[sample_rows, ]
        # to train the model
        DataCV=df[-sample_rows, ]
        # train the model
        gam1= gam(y ~ s(z, x), data=DataCV, family='binomial')
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
    for (i in 1:10)
    {
        # sample indices
        sample_rows = seq((i-1)*split+1, (i)*split)
        # to test the model
        DataC1=df[sample_rows, ]
        # to train the model
        DataCV=df[-sample_rows, ]
        # train the model
        gam1= gam(y ~ s(log(z), log(x)), data=DataCV, family='binomial')
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

write.table(results_dfII, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/predictions.both_classes.model_", args[3], ".", args[1], '.', args[2], ".data.txt", sep=''))
