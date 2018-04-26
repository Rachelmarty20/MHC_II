# Imports
library(parallel)
library(lme4)
library(mgcv)
library(xtable)

PATH_TO_DATA = '/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/'

args = commandArgs(trailingOnly=TRUE)
pan = as.integer(args[1])
class = args[2]
name = args[3]
mutation_threshold = as.integer(args[4])
tissue_file = args[5]
mut_file = args[6]
aff_file = args[7]


get_or <- function(fit) { c(exp(c(coef(fit)[2,1],coef(fit)[2,1]-1.96*coef(fit)[2,2],coef(fit)[2,1]+1.96*coef(fit)[2,2])),coef(fit)[2,4]) }

print(paste(pan, class, name, mutation_threshold))

# Import data
tissue <- read.csv(paste(PATH_TO_DATA, tissue_file, sep=""),header=TRUE)
mut <- read.csv(paste(PATH_TO_DATA, mut_file, sep=""),header=TRUE)
aff <- read.csv(paste(PATH_TO_DATA, aff_file, sep=""),header=TRUE)
print(dim(mut))
print(dim(aff))


if (pan == 1){
    tissue <- read.csv(paste(PATH_TO_DATA, tissue_file, sep=""),header=TRUE)
    mut <- read.csv(paste(PATH_TO_DATA, mut_file, sep=""),header=TRUE)
    aff <- read.csv(paste(PATH_TO_DATA, aff_file, sep=""),header=TRUE)
    if (class == 'random'){
        aff <- aff[sample(nrow(aff)),]
    }
    print(dim(mut))
    print(dim(aff))

    patient <- as.character(mut[,1])
    mut <- as.matrix(mut[,-1])
    aff <- as.matrix(aff[,-1])
    rownames(mut) <- rownames(aff) <- patient

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
    #tissuetypes <- as.character(unique(tissue[,2]))
    tissuetypes <- c('GBM','LUAD','LUSC','BLCA','PAAD','COAD','STAD','SKCM',
                    'THCA','HNSC','READ','LGG', 'BRCA','OV','PRAD')
    mysummary0 <- mysummary1 <- mysummary2 <- vector("list",length(tissuetypes))
    names(mysummary0) <- names(mysummary1) <- names(mysummary2) <- tissuetypes


    for (i in 1:length(tissuetypes)) {
        cat("TISSUE",tissuetypes[i])

        # restrict by tissue
        tissue2 <- tissue[tissue$X %in% as.vector(tissue$X[tissue$Tissue==tissuetypes[i]]), ]
        aff2 <- aff[aff$X %in% as.vector(tissue$X[tissue$Tissue==tissuetypes[i]]), ]
        mut2 <- mut[mut$X %in% as.vector(tissue$X[tissue$Tissue==tissuetypes[i]]), ]

        # randomize
        if (class == 'random'){
            aff2 <- aff2[sample(nrow(aff2)),]
        }

        # getting everything in the right format
        patient <- as.character(mut2[,1])
        mut2 <- as.matrix(mut2[,-1])
        aff2 <- as.matrix(aff2[,-1])
        rownames(mut2) <- rownames(aff2) <- patient

        # make dataframe
        y= as.vector(mut2); x= as.vector(aff2)
        gene= rep(colnames(mut2),each=nrow(mut2))
        pat= rep(rownames(mut2),ncol(mut2))
        nmut= colSums(mut2)
        genesel= (gene %in% names(nmut[nmut>=mutation_threshold]))

        patsel= pat %in% as.character(tissue2$X[tissue2$Tissue==tissuetypes[i]])
        sel= genesel & patsel
        #
        lme2= glmer(y[sel] ~ log(x[sel]) + (1|pat[sel]), family='binomial')
        mysummary2[[i]] <- summary(lme2)
        cat("Done \n")
    }

    tabpat <- do.call(rbind,lapply(mysummary2,get_or))
    colnames(tabpat) <- c('OR', 'Lci', 'Hci', 'P')
    write.table(tabpat, sep=',',
                file=paste("/cellar/users/ramarty/Data/hla_ii/generated_data/OR_clean/tissue/", class, "/", name, ".thresh_", mutation_threshold, ".txt", sep=''))


}