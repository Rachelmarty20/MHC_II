# Imports
library(parallel)
library(lme4)
library(mgcv)
library(xtable)

PATH_TO_DATA = '/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/'

args = commandArgs(trailingOnly=TRUE)
mutation_threshold = as.integer(args[1])
class = args[2] # this is for whether it is MHC-II or both
tissue_file = args[3]
mut_file = args[4]
aff1_file = args[5]
name = args[6]
pan = as.integer(args[7])

get_or <- function(fit) { c(exp(c(coef(fit)[2,1],coef(fit)[2,1]-1.96*coef(fit)[2,2],coef(fit)[2,1]+1.96*coef(fit)[2,2])),coef(fit)[2,4]) }

print(paste(pan, class, name, mutation_threshold))


# Import data
tissue <- read.csv(paste(PATH_TO_DATA, tissue_file, sep=""),header=TRUE)
mut <- read.csv(paste(PATH_TO_DATA, mut_file, sep=""),header=TRUE)
aff <- read.csv(paste(PATH_TO_DATA, aff1_file, sep=""),header=TRUE)
print(dim(mut))
print(dim(aff))
patient <- as.character(mut[,1])
mut <- as.matrix(mut[,-1])
aff <- as.matrix(aff[,-1])
rownames(mut) <- rownames(aff) <- patient

if (pan == 1){
    # Format data
    y= as.vector(mut); x= as.vector(aff)
    gene= rep(colnames(mut),each=nrow(mut))
    pat= rep(rownames(mut),ncol(mut))
    nmut= colSums(mut)
    sel= gene %in% names(nmut[nmut>=mutation_threshold])

    lme2= glmer(y[sel] ~ log(x[sel]) + (1|pat[sel]), family='binomial')
    mysummarypan <- vector("list",1)
    mysummarypan[[1]] <- summary(lme2)
    tabgene <- do.call(rbind,lapply(mysummarypan,get_or))
    rownames(tabgene) <- c('mutation')
    colnames(tabgene) <- c('OR', "conf_OR_low", 'conf_OR_high', 'P')
    write.table(tabgene,
    file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/OR_clean/pan/", class, "/", name, ".thresh_", mutation_threshold, ".txt", sep=''))

}

if (pan == 0){
    ## Tissue-specific ##
    # Make dataframe
    y= as.vector(mut); x= as.vector(aff)
    gene= rep(colnames(mut),each=nrow(mut))
    pat= rep(rownames(mut),ncol(mut))
    nmut= colSums(mut)
    genesel= (gene %in% names(nmut[nmut>=mutation_threshold]))

    #tissuetypes <- as.character(unique(tissue[,2]))
    # 'OV','PRAD','BRCA','UCEC',
    tissuetypes <- c('GBM','LUAD','LUSC','BLCA','PAAD','COAD','STAD','SKCM',
                    'THCA','HNSC','READ','LGG')
    mysummary0 <- mysummary1 <- mysummary2 <- vector("list",length(tissuetypes))
    names(mysummary0) <- names(mysummary1) <- names(mysummary2) <- tissuetypes


    for (i in 1:length(tissuetypes)) {
        cat("TISSUE",tissuetypes[i])
        #
        patsel= pat %in% as.character(tissue$X[tissue$Tissue==tissuetypes[i]])
        sel= genesel & patsel
        #
        lme2= glmer(y[sel] ~ log(x[sel]) + (1|pat[sel]), family='binomial')
        mysummary2[[i]] <- summary(lme2)
        cat("Done \n")
    }

    tabpat <- do.call(rbind,lapply(mysummary2,get_or))
    colnames(tabpatI) <- c('OR', 'Lci', 'Hci', 'P')
    write.table(tabpat, sep=',',
                file=paste("/cellar/users/ramarty/Data/hla_ii/generated_data/OR_clean/tissue/", class, "/", name, ".thresh_", mutation_threshold, ".txt", sep=''))


}

