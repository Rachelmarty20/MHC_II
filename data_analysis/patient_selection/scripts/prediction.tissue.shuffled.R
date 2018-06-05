library(parallel)
library(lme4)
library(mgcv)
library(xtable)
library(pROC)

PATH_TO_DATA = '/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/'

args = commandArgs(trailingOnly=TRUE)
mutation_threshold = as.integer(args[1])
tissue_type = args[2]

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

# Subsetting tissue df to match mut df
patients = as.vector(row.names(mut))
tissue = subset(tissue, Sample %in% patients)

# These are added for the meaning change
tissue_patients = as.vector(tissue$Sample[tissue$Tissue==tissue_type])
tissue_mut = mut[tissue_patients,]

y= as.vector(mut); x= as.vector(aff1); z= as.vector(aff2)
gene= rep(colnames(mut),each=nrow(mut))
pat= rep(rownames(mut),ncol(mut))
# Changing this line changes the whole meaning!!!

nmut= colSums(tissue_mut)
genesel= gene %in% names(nmut[nmut>=mutation_threshold])
patsel= pat %in% as.character(tissue$Sample[tissue$Tissue==tissue_type])

sel= genesel & patsel
df = data.frame(y[sel], x[sel], z[sel], sample(pat[sel]))
colnames(df)<-c('y', 'x', 'z', 'pat')
# randomize patients
df <- transform( df, pat = sample(pat) )

# randomize order for CV
df <- df[sample(1:nrow(df)), ]
split = round(dim(df)[1] / 10)
print(dim(df))

#  MHC-I
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
    gam1= gam(y ~ s(x), data=DataCV, family='binomial')
    # predict mutation probabilities
    P1=predict(gam1, DataC1)
    names(P1)=NULL
    all_predictions= c(all_predictions, P1)
    all_labels = c(all_labels, DataC1$y)
}
results_df = data.frame(all_labels, all_predictions)
colnames(results_df)<-c('label', 'predicted')
results_df$predicted_prob<-exp(results_df$predicted)
results_df$label_fact <- factor(results_df$label)
write.table(results_df, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/tissues/predictions.muts_in_", tissue_type, ".MHC_I.", mutation_threshold, ".data.shuffled.txt", sep=''))


#  MHC-II
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
    gam1= gam(y ~ s(z), data=DataCV, family='binomial')
    # predict mutation probabilities
    P1=predict(gam1, DataC1)
    names(P1)=NULL
    all_predictions= c(all_predictions, P1)
    all_labels = c(all_labels, DataC1$y)
}
results_df = data.frame(all_labels, all_predictions)
colnames(results_df)<-c('label', 'predicted')
results_df$predicted_prob<-exp(results_df$predicted)
results_df$label_fact <- factor(results_df$label)
write.table(results_df, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/tissues/predictions.muts_in_", tissue_type, ".MHC_II.", mutation_threshold, ".data.shuffled.txt", sep=''))


#  Both
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
results_df = data.frame(all_labels, all_predictions)
colnames(results_df)<-c('label', 'predicted')
results_df$predicted_prob<-exp(results_df$predicted)
results_df$label_fact <- factor(results_df$label)
write.table(results_df, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/tissues/predictions.muts_in_", tissue_type, ".Both.", mutation_threshold, ".data.shuffled.txt", sep=''))

