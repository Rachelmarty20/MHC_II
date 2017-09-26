library(parallel)
library(lme4)
library(mgcv)
library(xtable)
get_or <- function(fit) { c(exp(c(coef(fit)[2,1],coef(fit)[2,1]-1.96*coef(fit)[2,2],coef(fit)[2,1]+1.96*coef(fit)[2,2])),coef(fit)[2,4]) }

args = commandArgs(trailingOnly=TRUE)

#Format data
tissue <- read.csv('/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/patient_tissues.csv',header=TRUE)
mut <- read.csv(paste('/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/patient_mutations.cancer.', args[1], '.csv', sep=''),header=TRUE)
aff <- read.csv(paste('/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/patient_affinities.cancer.', args[1], '.csv', sep=''),header=TRUE)
patient <- as.character(mut[,1])
mut <- as.matrix(mut[,-1])
aff <- as.matrix(aff[,-1])
rownames(mut) <- rownames(aff) <- patient
print('Data formatted.')

######################################################################################
# ANALYSIS 1. ALL CANCER TYPES COMBINED
######################################################################################

y= as.vector(mut); x= as.vector(aff)
gene= rep(colnames(mut),each=nrow(mut))
pat= rep(rownames(mut),ncol(mut))
nmut= colSums(mut)
sel= gene %in% names(nmut[nmut>=5])

lme1= glmer(y[sel] ~ log(x[sel]) + (1|gene[sel]), family='binomial')
summary(lme1)

lme2= glmer(y[sel] ~ log(x[sel]) + (1|pat[sel]), family='binomial')
summary(lme2)

print('Models created.')

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/globalassoc_gam_harm.', args[1], '.pdf', sep=''))
#thre <- c(0,0.5,1,1.5,2,2.5,3,4,5,Inf)
#xr= cut(x,breaks=thre)
#m <- tapply(y,xr,'mean')
thre <- c(seq(0,30,length=20),Inf)
xr= cut(x,breaks=thre)
m <- tapply(y,xr,'mean')
plot(log(tapply(x,xr,'mean')),log(m/(1-m)),xlab='log-affinity',ylab='logit probability of mutation')

#
gam1= gam(y[sel] ~ s(log(x[sel])), family='binomial')
ypred= predict(gam1,type='response',se.fit=TRUE)
o= order(x[sel])
plot(x[sel][o],ypred$fit[o],type='l',xlim=c(0,10),xlab='Affinity',ylab='Probability of mutation',main='Generalized additive model')
lines(x[sel][o],ypred$fit[o]-1.96*ypred$se.fit[o],lty=2)
lines(x[sel][o],ypred$fit[o]+1.96*ypred$se.fit[o],lty=2)
#
plot(gam1,rug=FALSE,xlab='log affinity',ylab='logit mutation probability',main='Estimated logit(prob) vs log(affinity)')

dev.off()

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/globalassoc_descriptive.', args[1], '.pdf', sep=''))
boxplot(x ~ y, outline=F, ylab='Affinity', xlab='Mutation')
#
hist(x[y==0],main='',xlim=c(0,10),ylim=c(0,1),prob=T,breaks=seq(0,100,by=.5),xlab='Affinity',ylab='Frequency',cex.lab=1.3,cex.axis=1.3)
par(new=TRUE)
hist(x[y==1],main='',xlim=c(0,10),ylim=c(0,1),prob=T,breaks=seq(0,100,by=.5),border=2,xaxt='n',yaxt='n',xlab='',ylab='')
legend('topright',c('No mutation','Mutation'),lty=1,col=1:2,cex=1.3)
dev.off()

print('Pan completed created.')



######################################################################################
# ANALYSIS 2. SEPARATELY FOR CANCER TYPES
######################################################################################

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

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/oddsratio_withingene_cancertype.', args[1], '.pdf', sep=''))
x2plot <- tabgene[order(tabgene[,1]),]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/oddsratio_withinpat_cancertype.', args[1], '.pdf', sep=''))
x2plot <- tabpat[order(tabpat[,1]),]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()

npat <- table(tissue$Tissue)

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/oddsratio_withingene_cancertype_100pat.', args[1], '.pdf', sep=''))
x2plot <- tabgene[order(tabgene[,1]),]
x2plot <- x2plot[rownames(x2plot) %in% names(npat[npat>=100]), ]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/oddsratio_withinpat_cancertype_100pat.', args[1], '.pdf', sep=''))
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

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/oddsratio_withingene_cancertype.20.', args[1], '.pdf', sep=''))
x2plot <- tabgene[order(tabgene[,1]),]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/oddsratio_withinpat_cancertype.20.', args[1], '.pdf', sep=''))
x2plot <- tabpat[order(tabpat[,1]),]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()

npat <- table(tissue$Tissue)

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/oddsratio_withingene_cancertype_100pat.20.', args[1], '.pdf', sep=''))
x2plot <- tabgene[order(tabgene[,1]),]
x2plot <- x2plot[rownames(x2plot) %in% names(npat[npat>=100]), ]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/oddsratio_withinpat_cancertype_100pat.20.', args[1], '.pdf', sep=''))
x2plot <- tabpat[order(tabpat[,1]),]
x2plot <- x2plot[rownames(x2plot) %in% names(npat[npat>=100]), ]
plot(x2plot[,1],1:nrow(x2plot),pch=15,xlim=c(0.25,4),yaxt='n',ylab='',xlab='Odds-ratio',log='x')
segments(x0=x2plot[,2],x1=x2plot[,3],1:nrow(x2plot),lty=2)
text(x2plot[,2],1:nrow(x2plot),rownames(x2plot),pos=2)
abline(v=1,col='gray',lwd=2)
dev.off()
