library(parallel)
library(lme4)
library(mgcv)
library(xtable)

PATH_TO_DATA = '/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/'

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
sel= gene %in% names(nmut[nmut>=2])

gamlog10 = gam(y[sel] ~ s(log10(x[sel]), log10(z[sel])), family='binomial')

pdf('/cellar/users/ramarty/Data/hla_ii/generated_figures/class_comparison/mutation_probability_density.pdf',
      width=2.5,height=2.5)

plot(gamlog10, rug=FALSE, all.terms=TRUE, scheme=c(2,1), ylim=c(-2,2), xlim=c(-2,2),
     xlab='Log(MHC-I PHBR)',ylab='Log(MHC-II PHBR)', main='Mutation Probability')
 dev.off()
