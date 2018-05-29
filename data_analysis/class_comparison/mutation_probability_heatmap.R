library(parallel)
library(lme4)
library(mgcv)
library(xtable)

PATH_TO_DATA = '/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/'

#Format data
tissue <- read.csv(paste(PATH_TO_DATA, 'patient_tissues.conservative.csv', sep=""),header=TRUE)
mut <- read.csv(paste(PATH_TO_DATA, 'patient_mutations.cancer.TCGA.conservative.mut.csv', sep=""),header=TRUE)
aff1 <- read.csv(paste(PATH_TO_DATA, 'patient_affinities.cancer.TCGA.conservative.mut.ClassI.csv', sep=""),header=TRUE)
aff2 <- read.csv(paste(PATH_TO_DATA, 'patient_affinities.cancer.TCGA.conservative.mut.ClassII.csv', sep=""),header=TRUE)
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

pdf('/cellar/users/ramarty/Data/hla_ii/generated_figures/class_comparison/mutation_probability_density.conservative.pdf',
      width=3,height=3)
par(mar=c(3.5,3,1.1,1.5)-0.3, mgp=c(1.25,0.25,0),las=1)
plot(gamlog10, rug=FALSE, all.terms=TRUE, scheme=c(2,1), ylim=c(-2,2), xlim=c(-2,2),
     xlab='Log(MHC-I)',ylab='Log(MHC-II)', main=' ',
     cex.axis=.75, cex.lab=0.75, tck=0)
dev.off()
