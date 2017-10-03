library(parallel)
library(lme4)
library(mgcv)
library(xtable)
get_or <- function(fit) { c(exp(c(coef(fit)[2,1],coef(fit)[2,1]-1.96*coef(fit)[2,2],coef(fit)[2,1]+1.96*coef(fit)[2,2])),coef(fit)[2,4]) }

args = commandArgs(trailingOnly=TRUE)
# args1 = genes
# args2 = category
# args3 = mutation matrix
# args4 = affinity matrix

print(args[1])
print(args[2])
print(args[3])
print(args[4])

#Format data
mut <- read.csv(args[3], header=TRUE)
aff <- read.csv(args[4], header=TRUE)
patient <- as.character(mut[,1])
mut <- as.matrix(mut[,-1])
aff <- as.matrix(aff[,-1])
rownames(mut) <- rownames(aff) <- patient
print('Data formatted.')


y= as.vector(mut); x= as.vector(aff)
gene= rep(colnames(mut),each=nrow(mut))
pat= rep(rownames(mut),ncol(mut))
nmut= colSums(mut)
sel= gene %in% names(nmut[nmut>=5])

lme1= glmer(y[sel] ~ log(x[sel]) + (1|gene[sel]), family='binomial')
summary(lme1)

lme2= glmer(y[sel] ~ log(x[sel]) + (1|pat[sel]), family='binomial')
summary(lme2)

mysummarypan <- vector("list",2)
mysummarypan[[1]] <- summary(lme1)
mysummarypan[[2]] <- summary(lme2)
tabgene <- do.call(rbind,lapply(mysummarypan,get_or))
rownames(tabgene) <- c('mutation', 'patient')
colnames(tabgene) <- c('OR', "conf_OR_low", 'conf_OR_high', 'P')
write.table(tabgene, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/pan.", args[1], ".", args[2], ".txt", sep=''))

print('Models created.')

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/', args[1], '/globalassoc_gam_harm.', args[2], '.pdf', sep=''))
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
plot(x[sel][o],ypred$fit[o],type='l',xlim=c(0,30),xlab='Affinity',ylab='Probability of mutation',main='Generalized additive model')
lines(x[sel][o],ypred$fit[o]-1.96*ypred$se.fit[o],lty=2)
lines(x[sel][o],ypred$fit[o]+1.96*ypred$se.fit[o],lty=2)
#
plot(gam1,rug=FALSE,xlab='log affinity',ylab='logit mutation probability',main='Estimated logit(prob) vs log(affinity)')

dev.off()

pdf(paste('/cellar/users/ramarty/Data/hla_ii/generated_figures/',  args[1], '/globalassoc_descriptive.', args[2], '.pdf', sep=''))
boxplot(x ~ y, outline=F, ylab='Affinity', xlab='Mutation')
#
hist(x[y==0],main='',xlim=c(0,100),ylim=c(0,0.05),prob=T,breaks=seq(0,100,by=5),xlab='Affinity',ylab='Frequency',cex.lab=1.3,cex.axis=1.3)
par(new=TRUE)
hist(x[y==1],main='',xlim=c(0,100),ylim=c(0,0.05),prob=T,breaks=seq(0,100,by=5),border=2,xaxt='n',yaxt='n',xlab='',ylab='')
legend('topright',c('No mutation','Mutation'),lty=1,col=1:2,cex=1.3)
dev.off()

print('Pan completed created.')
