# Imports
library(parallel)
library(lme4)
library(mgcv)
library(xtable)
library(pROC)
library(oddsratio)

PATH_TO_DATA = '/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/'

args = commandArgs(trailingOnly=TRUE)
mutation_threshold = as.integer(args[1])
model = as.integer(args[2]) # this is for whether it is MHC-II or both
tissue_file = args[3]
mut_file = args[4]
aff1_file = args[5]
aff2_file = args[6]
name = args[7]
pan = args[8]


### Model with MHC-II ###
if (model == 0){
    # Import data
    tissue <- read.csv(paste(PATH_TO_DATA, tissue_file, sep=""),header=TRUE)
    mut <- read.csv(paste(PATH_TO_DATA, mut_file, sep=""),header=TRUE)
    aff <- read.csv(paste(PATH_TO_DATA, aff1_file, sep=""),header=TRUE)
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

        ## Pan-cancer ##
        # Make dataframe
        sampled_pats = unique(pat)
        sel_pat = pat %in% sampled_pats
        df = data.frame(y[sel&sel_pat], x[sel&sel_pat], pat[sel&sel_pat])
        colnames(df)<-c('y', 'x', 'pat')

        # Train model
        gam = gam(y ~ s(x), data=df, family='binomial')
        low_x = quantile(x, 0.25, names=FALSE)
        high_x = quantile(x, 0.75, names=FALSE)
        results = or_gam(data = df, model = gam, pred = c("x"), values=c(low_x, high_x))

        # Format output
        write.table(results, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/OR/MHC_II.pan.thresh_", mutation_threshold, ".", name, ".txt", sep=''))
    }

    if (pan == 0){
        ## Tissue-specific ##
        # Make dataframe
        y= as.vector(mut); x= as.vector(aff)
        gene= rep(colnames(mut),each=nrow(mut))
        pat= rep(rownames(mut),ncol(mut))
        nmut= colSums(mut)
        genesel= (gene %in% names(nmut[nmut>=mutation_threshold]))

        # train models
        tissuetypes <- c('MESO', 'BRCA', 'UCS', 'LUSC', 'GBM', 'READ', 'KICH', 'COAD', 'SKCM', 'STAD', 'THCA', 'PRAD', 'CESC', 'BLCA', 'UVM', 'ACC', 'LGG', 'UCEC', 'TGCT', 'OV', 'LAML', 'LUAD', 'LIHC', 'HNSC', 'PCPG', 'KIRP', 'DLBC', 'KIRC', 'PAAD')
        OR <- CI_low <- CI_high <- predicted <- vector("list",length(tissuetypes))
        for (i in 1:length(tissuetypes)) {
            cat("TISSUE",tissuetypes[i])
            #
            patsel= pat %in% as.character(tissue$Sample[tissue$Tissue==tissuetypes[i]])
            sel= genesel & patsel

            df = data.frame(y[sel&patsel], x[sel&patsel], pat[sel&patsel])
            colnames(df)<-c('y', 'x', 'pat')

            gam = gam(y ~ s(x), data=df, family='binomial')
            low_x = quantile(x, 0.25, names=FALSE)
            high_x = quantile(x, 0.75, names=FALSE)
            results = or_gam(data = df, model = gam, pred = c("x"), values=c(low_x, high_x))
            OR[[i]] <- results[['oddsratio']]
            CI_low[[i]] <- results[['CI_low (2.5%)']]
            CI_high[[i]] <- results[['CI_high (97.5%)']]
            predicted[[i]] <- results[['predictor']]

            cat("Done \n")
        }

        # Format output
        results = cbind(cbind(cbind(cbind(OR, CI_low), CI_high), predicted), tissuetypes)
        write.table(results, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/OR/MHC_II.tissue.thresh_", mutation_threshold, ".", name, ".txt", sep=''))
    }
}


### Both MHC-I and MHC-II
if (model == 1){
    # Import data
    tissue <- read.csv(paste(PATH_TO_DATA, tissue_file, sep=""),header=TRUE)
    mut <- read.csv(paste(PATH_TO_DATA, mut_file, sep=""),header=TRUE)
    aff1 <- read.csv(paste(PATH_TO_DATA, aff1_file, sep=""),header=TRUE)
    aff2 <- read.csv(paste(PATH_TO_DATA, aff2_file, sep=""),header=TRUE)
    patient <- as.character(mut[,1])
    mut <- as.matrix(mut[,-1])
    aff1 <- as.matrix(aff1[,-1])
    aff2 <- as.matrix(aff2[,-1])
    rownames(mut) <- rownames(aff1) <- rownames(aff2) <- patient

    if (pan == 1){
        # Format data
        y= as.vector(mut); x= as.vector(aff1);  z= as.vector(aff2)
        gene= rep(colnames(mut),each=nrow(mut))
        pat= rep(rownames(mut),ncol(mut))
        nmut= colSums(mut)
        sel= gene %in% names(nmut[nmut>=mutation_threshold])

        ## Pan-cancer ##
        # Make dataframe
        sampled_pats = unique(pat)
        sel_pat = pat %in% sampled_pats
        df = data.frame(y[sel&sel_pat], x[sel&sel_pat], z[sel&sel_pat], pat[sel&sel_pat])
        colnames(df)<-c('y', 'x', 'z', 'pat')

        # Train model
        gam = gam(y ~ s(x, z), data=df, family='binomial')
        low_x = quantile(x, 0.25, names=FALSE)
        high_x = quantile(x, 0.75, names=FALSE)
        low_z = quantile(z, 0.25, names=FALSE)
        high_z = quantile(z, 0.75, names=FALSE)
        results1 = or_gam(data = df, model = gam, pred = c("x"), values=c(low_x, high_x))
        results2 = or_gam(data = df, model = gam, pred = c("z"), values=c(low_z, high_z))

        OR <- CI_low <- CI_high <- predicted <- tissue <- vector("list",2)

        OR[[1]] <- results1[['oddsratio']]
        CI_low[[1]] <- results1[['CI_low (2.5%)']]
        CI_high[[1]] <- results1[['CI_high (97.5%)']]
        predicted[[1]] <- results1[['predictor']]

        OR[[2]] <- results2[['oddsratio']]
        CI_low[[2]] <- results2[['CI_low (2.5%)']]
        CI_high[[2]] <- results2[['CI_high (97.5%)']]
        predicted[[2]] <- results2[['predictor']]

        # Format output
        results = cbind(cbind(cbind(OR, CI_low), CI_high), predicted)
        write.table(results, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/OR/Both.pan.thresh_", mutation_threshold, ".", name, ".txt", sep=''))
    }

    if (pan == 0){
        ## Tissue-specific ##
        # Make dataframe
        y= as.vector(mut); x= as.vector(aff)
        gene= rep(colnames(mut),each=nrow(mut))
        pat= rep(rownames(mut),ncol(mut))
        nmut= colSums(mut)
        genesel= (gene %in% names(nmut[nmut>=mutation_threshold]))

        # train models
        tissuetypes <- c('MESO', 'BRCA', 'UCS', 'LUSC', 'GBM', 'READ', 'KICH', 'COAD', 'SKCM', 'STAD', 'THCA', 'PRAD', 'CESC', 'BLCA', 'UVM', 'ACC', 'LGG', 'UCEC', 'TGCT', 'OV', 'LAML', 'LUAD', 'LIHC', 'HNSC', 'PCPG', 'KIRP', 'DLBC', 'KIRC', 'PAAD')
        OR <- CI_low <- CI_high <- predicted <- tissue <- vector("list",length(tissuetypes)*2)
        for (i in 1:length(tissuetypes)) {
            cat("TISSUE",tissuetypes[i])
            #
            patsel= pat %in% as.character(tissue$Sample[tissue$Tissue==tissuetypes[i]])
            sel= genesel & patsel

            df = data.frame(y[sel&patsel], x[sel&patsel], pat[sel&patsel])
            colnames(df)<-c('y', 'x', 'pat')

            gam = gam(y ~ s(x), data=df, family='binomial')
            low_x = quantile(df[['x']], 0.25, names=FALSE)
            high_x = quantile(df[['x']], 0.75, names=FALSE)
            results1 = or_gam(data = df, model = gam, pred = c("x"), values=c(low_x, high_x))

            OR[[i]] <- results1[['oddsratio']]
            CI_low[[i]] <- results1[['CI_low (2.5%)']]
            CI_high[[i]] <- results1[['CI_high (97.5%)']]
            predicted[[i]] <- results1[['predictor']]
            tissue[[i]] <- tissuetypes[i]

            OR[[length(tissuetypes)+i]] <- results2[['oddsratio']]
            CI_low[[length(tissuetypes)+i]] <- results2[['CI_low (2.5%)']]
            CI_high[[length(tissuetypes)+i]] <- results2[['CI_high (97.5%)']]
            predicted[[length(tissuetypes)+i]] <- results2[['predictor']]
            tissue[[length(tissuetypes)+i]] <- tissuetypes[i]

            cat("Done \n")
        }

        # Format output
        results = cbind(cbind(cbind(cbind(OR, CI_low), CI_high), predicted), tissue)
        write.table(results, file = paste("/cellar/users/ramarty/Data/hla_ii/generated_data/OR/Both.tissue.thresh_", mutation_threshold, ".", name, ".txt", sep=''))
    }
}
