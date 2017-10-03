library(parallel)
library(lme4)
library(mgcv)
library(xtable)
get_or <- function(fit) { c(exp(c(coef(fit)[2,1],coef(fit)[2,1]-1.96*coef(fit)[2,2],coef(fit)[2,1]+1.96*coef(fit)[2,2])),coef(fit)[2,4]) }

args = commandArgs(trailingOnly=TRUE)
# args1 = class

#Format data
tissue <- read.csv('/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/patient_tissues.csv',header=TRUE)
mut <- read.csv('/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/combined_classes/patient_mutations.csv',header=TRUE)
aff <- read.csv(paste('/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/combined_classes/patient_affinities.', args[1], '.csv', sep=''),header=TRUE)
patient <- as.character(mut[,1])
mut <- as.matrix(mut[,-1])
aff <- as.matrix(aff[,-1])
rownames(mut) <- rownames(aff) <- patient
print('Data formatted.')

y= as.vector(mut); x= as.vector(aff)
gene= rep(colnames(mut),each=nrow(mut))
pat= rep(rownames(mut),ncol(mut))
nmut= colSums(mut)
genesel= (gene %in% names(nmut[nmut>=5]))

tissuetypes <- c('MESO', 'BRCA', 'UCS', 'LUSC', 'GBM', 'READ', 'KICH', 'COAD', 'SKCM', 'STAD', 'THCA', 'PRAD', 'CESC', 'BLCA', 'UVM', 'ACC', 'LGG', 'UCEC', 'TGCT', 'OV', 'LAML', 'LUAD', 'LIHC', 'HNSC', 'PCPG', 'KIRP', 'DLBC', 'KIRC', 'PAAD')
#tissuetypes <- as.character(unique(tissue[,2]))
mysummary0 <- mysummary1 <- mysummary2 <- vector("list",length(tissuetypes))
names(mysummary0) <- names(mysummary1) <- names(mysummary2) <- tissuetypes
for (i in 1:length(tissuetypes)) {
    cat("TISSUE",tissuetypes[i])
    #
    patsel= pat %in% as.character(tissue$Sample[tissue$Tissue==tissuetypes[i]])
    sel= genesel & patsel
    #
    lme0= glm(y[sel] ~ log(x[sel]), family='binomial')
    mysummary0[[i]] <- summary(lme0)
    #
    lme1= glmer(y[sel] ~ log(x[sel]) + (1|gene[sel]), family='binomial')
    mysummary1[[i]] <- summary(lme1)
    #
    lme2= glmer(y[sel] ~ log(x[sel]) + (1|pat[sel]), family='binomial')
    mysummary2[[i]] <- summary(lme2)
    cat("Done \n")
}

save(mysummary0,mysummary1,mysummary2,file='/cellar/users/ramarty/Data/hla_ii/generated_data/fit_cancertype.RData')

print('Tissues >= 5 completed.')

load('/cellar/users/ramarty/Data/hla_ii/generated_data/fit_cancertype.RData')
tabgene <- do.call(rbind,lapply(mysummary1,get_or))
tabpat <- do.call(rbind,lapply(mysummary2,get_or))

xtable(tabgene[order(rownames(tabgene)),],digits=c(0,3,3,3,4))
xtable(tabpat[order(rownames(tabpat)),],digits=c(0,3,3,3,4))

colnames(tabgene) <- c('OR', "conf_OR_low", 'conf_OR_high', 'P')
colnames(tabpat) <- c('OR', "conf_OR_low", 'conf_OR_high', 'P')
write.table(tabgene, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/matched_models/tissues.5mut.", args[1], ".txt", sep=''))
write.table(tabpat, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/matched_models/tissues.5mut.", args[1], ".txt", sep=''))


pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/', args[1], '/oddsratio_withingene_cancertype.pdf', sep=''))
x2plot <- tabgene[order(tabgene[,1]),]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/', args[1], '/oddsratio_withinpat_cancertype.pdf', sep=''))
x2plot <- tabpat[order(tabpat[,1]),]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()

npat <- table(tissue$Tissue)

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/', args[1], '/oddsratio_withingene_cancertype_100pat.pdf', sep=''))
x2plot <- tabgene[order(tabgene[,1]),]
x2plot <- x2plot[rownames(x2plot) %in% names(npat[npat>=100]), ]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/', args[1], '/oddsratio_withinpat_cancertype_100pat.pdf', sep=''))
x2plot <- tabpat[order(tabpat[,1]),]
x2plot <- x2plot[rownames(x2plot) %in% names(npat[npat>=100]), ]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()


y= as.vector(mut); x= as.vector(aff)
gene= rep(colnames(mut),each=nrow(mut))
pat= rep(rownames(mut),ncol(mut))
nmut= colSums(mut)
genesel= (gene %in% names(nmut[nmut>=20]))

#tissuetypes <- c('HNSC', 'LUAD', 'PRAD', 'LUSC', 'LGG', 'BRCA', 'GBM', 'STAD')
#tissuetypes <- as.character(unique(tissue[,2]))
mysummary0 <- mysummary1 <- mysummary2 <- vector("list",length(tissuetypes))
names(mysummary0) <- names(mysummary1) <- names(mysummary2) <- tissuetypes
for (i in 1:length(tissuetypes)) {
    cat("TISSUE",tissuetypes[i])
    #
    patsel= pat %in% as.character(tissue$Sample[tissue$Tissue==tissuetypes[i]])
    sel= genesel & patsel
    #
   if (length(unique(y[sel])) < 2) next
   if (min(table(y[sel])) == 1) next
	#
    lme0= glm(y[sel] ~ log(x[sel]), family='binomial')
    mysummary0[[i]] <- summary(lme0)
    #
    lme1= glmer(y[sel] ~ log(x[sel]) + (1|gene[sel]), family='binomial')
    mysummary1[[i]] <- summary(lme1)
    #
    lme2= glmer(y[sel] ~ log(x[sel]) + (1|pat[sel]), family='binomial')
    mysummary2[[i]] <- summary(lme2)
    cat("Done \n")
}

save(mysummary0,mysummary1,mysummary2,file='/cellar/users/ramarty/Data/hla_ii/generated_data/fit_cancertype.20.RData')

print('Tissues >= 20 completed.')

load('/cellar/users/ramarty/Data/hla_ii/generated_data/fit_cancertype.20.RData')
tabgene <- do.call(rbind,lapply(mysummary1,get_or))
tabpat <- do.call(rbind,lapply(mysummary2,get_or))

xtable(tabgene[order(rownames(tabgene)),],digits=c(0,3,3,3,4))
xtable(tabpat[order(rownames(tabpat)),],digits=c(0,3,3,3,4))

colnames(tabgene) <- c('OR', "conf_OR_low", 'conf_OR_high', 'P')
colnames(tabpat) <- c('OR', "conf_OR_low", 'conf_OR_high', 'P')
write.table(tabgene, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/matched_models/tissues.20mut.", args[1], ".txt", sep=''))
write.table(tabpat, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/matched_models/tissues.20mut.", args[1], ".txt", sep=''))

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/', args[1], '/oddsratio_withingene_cancertype.20.pdf', sep=''))
x2plot <- tabgene[order(tabgene[,1]),]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/', args[1], '/oddsratio_withinpat_cancertype.20.pdf', sep=''))
x2plot <- tabpat[order(tabpat[,1]),]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()

npat <- table(tissue$Tissue)


pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/', args[1], '/oddsratio_withingene_cancertype_100pat.20.pdf', sep=''))
x2plot <- tabgene[order(tabgene[,1]),]
x2plot <- x2plot[rownames(x2plot) %in% names(npat[npat>=100]), ]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/matched_models/', args[1], '/oddsratio_withinpat_cancertype_100pat.20.pdf', sep=''))
x2plot <- tabpat[order(tabpat[,1]),]
x2plot <- x2plot[rownames(x2plot) %in% names(npat[npat>=100]), ]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()
