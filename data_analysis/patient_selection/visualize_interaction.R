library(parallel)
library(lme4)
library(mgcv)
library(xtable)
get_or <- function(fit) { c(exp(c(coef(fit)[2,1],coef(fit)[2,1]-1.96*coef(fit)[2,2],coef(fit)[2,1]+1.96*coef(fit)[2,2])),coef(fit)[2,4]) }

args = commandArgs(trailingOnly=TRUE)
mutation_threshold = as.integer(args[1])

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
print('Data formatted.')


y= as.vector(mut); x= as.vector(aff1); z= as.vector(aff2)
gene= rep(colnames(mut),each=nrow(mut))
pat= rep(rownames(mut),ncol(mut))
nmut= colSums(mut)
sel= gene %in% names(nmut[nmut>=mutation_threshold])

# MHC-I
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam.mhc_i.threshold_', mutation_threshold, '.pdf', sep=''))
gam1= gam(y[sel] ~ s(log(x[sel])), family='binomial')
ypred= predict(gam1,type='response',se.fit=TRUE)
o= order(x[sel])
plot(x[sel][o],ypred$fit[o],type='l',xlim=c(1,100),xlab='MHC-I PHBR',ylab='Mutation probability',main='Generalized additive model')
lines(x[sel][o],ypred$fit[o]-1.96*ypred$se.fit[o],lty=2)
lines(x[sel][o],ypred$fit[o]+1.96*ypred$se.fit[o],lty=2)
dev.off()
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam_log.mhc_i.threshold_', mutation_threshold, '.pdf', sep=''))
plot(gam1,rug=FALSE,xlab='log MHC-I PHBR',ylab='logit mutation probability',main='Estimated logit(prob) vs log(affinity)')
dev.off()


# MHC-II
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam.mhc_ii.threshold_', mutation_threshold, '.pdf', sep=''))
gam2= gam(y[sel] ~ s(log(z[sel])), family='binomial')
ypred= predict(gam2,type='response',se.fit=TRUE)
o= order(z[sel])
plot(z[sel][o],ypred$fit[o],type='l',xlim=c(1,100),xlab='MHC-II PHBR',ylab='Mutation probability',main='Generalized additive model')
lines(z[sel][o],ypred$fit[o]-1.96*ypred$se.fit[o],lty=2)
lines(z[sel][o],ypred$fit[o]+1.96*ypred$se.fit[o],lty=2)
dev.off()
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam_log.mhc_ii.threshold_', mutation_threshold, '.pdf', sep=''))
plot(gam2,rug=FALSE,xlab='log MHC-II PHBR',ylab='logit mutation probability',main='Estimated logit(prob) vs log(affinity)')
dev.off()

# MHC-I / MHC-II

## log-log
gam_interact= gam(y[sel] ~ s(log(x[sel]), log(z[sel])), family='binomial')
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam.log_mhc_i_log_mhc_ii.threshold_', mutation_threshold, '.surface.pdf', sep=''))
plot(gam_interact, all.terms=TRUE, scheme=1, ylim=c(0,4.5), xlim=c(-2,2),
     xlab='Log(MHC-I PHBR)',ylab='Log(MHC-II PHBR)', main='Mutation Probability')
dev.off()
'''
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam.log_mhc_i_log_mhc_ii.threshold_', mutation_threshold, '.heatmap_points.pdf', sep=''))
plot(gam_interact, all.terms=TRUE, scheme=c(2,1), ylim=c(0,4.5), xlim=c(-2,2),
     xlab='Log(MHC-I PHBR)',ylab='Log(MHC-II PHBR)', main='Mutation Probability')
dev.off()
'''
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam.log_mhc_i_log_mhc_ii.threshold_', mutation_threshold, '.heatmap.pdf', sep=''))
plot(gam_interact, rug=FALSE, all.terms=TRUE, scheme=c(2,1), ylim=c(0,4.5), xlim=c(-2,2),
     xlab='Log(MHC-I PHBR)',ylab='Log(MHC-II PHBR)', main='Mutation Probability')
dev.off()

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam.log_mhc_i_log_mhc_ii.threshold_', mutation_threshold, '.surface.whole.pdf', sep=''))
plot(gam_interact, all.terms=TRUE, scheme=1,
     xlab='Log(MHC-I PHBR)',ylab='Log(MHC-II PHBR)', main='Mutation Probability')
dev.off()
'''
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam.log_mhc_i_log_mhc_ii.threshold_', mutation_threshold, '.heatmap_points.whole.pdf', sep=''))
plot(gam_interact, all.terms=TRUE, scheme=c(2,1),
     xlab='Log(MHC-I PHBR)',ylab='Log(MHC-II PHBR)', main='Mutation Probability')
dev.off()
'''
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam.log_mhc_i_log_mhc_ii.threshold_', mutation_threshold, '.heatmap.whole.pdf', sep=''))
plot(gam_interact, rug=FALSE, all.terms=TRUE, scheme=c(2,1),
     xlab='Log(MHC-I PHBR)',ylab='Log(MHC-II PHBR)', main='Mutation Probability')
dev.off()


# nonlog-nonlog
gam_interact_nonlog= gam(y[sel] ~ s(x[sel], z[sel]), family='binomial')
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam.mhc_i_mhc_ii.threshold_', mutation_threshold, '.surface.pdf', sep=''))
plot(gam_interact_nonlog, all.terms=TRUE, scheme=1, ylim=c(0,60), xlim=c(0,5),
     xlab='MHC-I PHBR',ylab='MHC-II PHBR', main='Mutation Probability')
dev.off()
'''
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam.mhc_i_mhc_ii.threshold_', mutation_threshold, '.heatmap_points.pdf', sep=''))
plot(gam_interact_nonlog, all.terms=TRUE, scheme=c(2,1), ylim=c(0,60), xlim=c(0,5),
     xlab='MHC-I PHBR',ylab='MHC-II PHBR', main='Mutation Probability')
dev.off()
'''
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam.mhc_i_mhc_ii.threshold_', mutation_threshold, '.heatmap.pdf', sep=''))
plot(gam_interact_nonlog, rug=FALSE, se=TRUE, scheme=2, ylim=c(0,60), xlim=c(0,5),
     xlab='MHC-I PHBR',ylab='MHC-II PHBR', main='Mutation Probability')
dev.off()

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam.mhc_i_mhc_ii.threshold_', mutation_threshold, '.surface.whole.pdf', sep=''))
plot(gam_interact_nonlog, all.terms=TRUE, scheme=1,
     xlab='MHC-I PHBR',ylab='MHC-II PHBR', main='Mutation Probability')
dev.off()
'''
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam.mhc_i_mhc_ii.threshold_', mutation_threshold, '.heatmap_points.whole.pdf', sep=''))
plot(gam_interact_nonlog, all.terms=TRUE, scheme=c(2,1),
     xlab='MHC-I PHBR',ylab='MHC-II PHBR', main='Mutation Probability')
dev.off()
'''
pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/gam.mhc_i_mhc_ii.threshold_', mutation_threshold, '.heatmap.whole.pdf', sep=''))
plot(gam_interact_nonlog, rug=FALSE, se=TRUE, scheme=2,
     xlab='MHC-I PHBR',ylab='MHC-II PHBR', main='Mutation Probability')
dev.off()
