library(parallel)
library(lme4)
library(mgcv)
library(xtable)
library(pROC)

PATH_TO_DATA = '/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/'

args = commandArgs(trailingOnly=TRUE)
input = args[1]
iteration = as.integer(args[2])


mutation_threshold = 2
model = 2

#Format data
tissue <- read.csv(paste(PATH_TO_DATA, 'patient_tissues.csv', sep=""),header=TRUE)
mut <- read.csv(paste(PATH_TO_DATA, 'patient_mutations.cancer.TCGA.inclusive.mut.csv', sep=""),header=TRUE)
aff1 <- read.csv(paste(PATH_TO_DATA, 'patient_affinities.cancer.TCGA.inclusive.mut.ClassI.csv', sep=""),header=TRUE)
aff2 <- read.csv(paste(PATH_TO_DATA, 'patient_affinities.cancer.TCGA.inclusive.mut.ClassII.csv', sep=""),header=TRUE)
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


# smoothed - PHBR
all_labels=NULL
all_predictions=NULL
i = 1
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


# put in dataframe and transform
results_dfII = data.frame(all_labels, all_predictions)
colnames(results_dfII)<-c('label', 'predicted')
# probabilities
results_dfII$predicted_prob<-exp(results_dfII$predicted)
# labels as factors
results_dfII$label_fact <- factor(results_dfII$label)

write.table(results_dfII, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/predictions/confidence_intervals/", args[1], "/iteration_", args[2], ".data.txt", sep=''))
