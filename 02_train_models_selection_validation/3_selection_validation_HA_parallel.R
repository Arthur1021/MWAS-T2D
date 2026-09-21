library(BEDMatrix)
library(data.table)
library(stringr)
library(dplyr)
library(RNOmni)
library(foreach)
library(doParallel)


# PatchUp
PatchUp <- function(M) {
    M <- apply(M, 2, function(x) {
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        return(x)
    })

    return(M)
}

# Standardize
Standardize <- function(M) {
    # Centralize
    M <- M - matrix(rep(colMeans(M), times = nrow(M)), nrow = nrow(M) , ncol = ncol(M), byrow = T)

    # Standardize
    M <- sweep(M, 2, sqrt(apply(M, 2, crossprod) / nrow(M)), "/")

    return(M)
}

# Check if the directory exists
if (!file.exists("6_HA_model_selection_validation")) {
  # Create the directory
  dir.create("6_HA_model_selection_validation")
}

finished_files <- Sys.glob('6_HA_model_selection_validation/*best*')
finished_files_df <- data.frame(finished_files)
finished_files_df$CpG <- sub(".*/(cg\\d+)_.*", "\\1", finished_files_df$finished_files)

#read pheno raw data
pheno_true <- fread('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/6_HA/methylation_m_value.txt', data.table =F)
rownames(pheno_true) <- pheno_true$ID
pheno_true <- pheno_true[,-1]

#pheno_true_selection <- pheno_true[,colnames(pheno_true) %in% rownames(geno_selection)]


for (chr_id in 1:22){
    #read geno for selection
    geno_selection <- BEDMatrix(paste0('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/6_HA/All_selection_chr',chr_id ), simple_names = TRUE)
    #impute missing geno with mean value
    geno_selection <- PatchUp(geno_selection)

    #Standardize
    geno_selection <- Standardize(geno_selection)

    #replace na with 0 
    geno_selection[is.na(geno_selection)] <- 0

    #bim for selection
    bim_selection <- fread(paste0('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/6_HA/All_selection_chr', chr_id, '.bim'), data.table =F)
    names(bim_selection) <- LETTERS[1:6]




    #cov for selection samples
    cov_selection1 <- fread('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/1_covariate/covariate.txt', data.table =F)
    rownames(cov_selection1) <- cov_selection1$NWDID
    cov_selection1 <- cov_selection1[,-1]
    cov_selection1 <- cov_selection1[rownames(cov_selection1) %in% rownames(geno_selection),]
    cov_selection1 <- cov_selection1[,-1]

    cov_selection_pca <- fread('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/6_HA/selection_PCA/combined_indep_PCA.pca.evec', data.table = F, skip = 1)
    rownames(cov_selection_pca) <- cov_selection_pca$V1
    cov_selection_pca <- cov_selection_pca[,-1]
    cov_selection_pca <- cov_selection_pca[,-11]
    colnames(cov_selection_pca) <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')

    cov_selection <- merge(cov_selection1, cov_selection_pca, by = 'row.names')
    rownames(cov_selection) <- cov_selection$Row.names
    cov_selection <- cov_selection[,-1]
    cov_selection <- cov_selection[,c(1:7,9:22)] #select sex age Bcell CD4T CD8T Eos Mono NK cig1c pkyrs1c BMI PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10

    #read geno for validation
    geno_validation <- BEDMatrix(paste0('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/6_HA/All_validation_chr',chr_id), simple_names = TRUE)
    #impute missing geno with mean value
    geno_validation <- PatchUp(geno_validation)

    #Standardize
    geno_validation <- Standardize(geno_validation)

    #replace na with 0 
    geno_validation[is.na(geno_validation)] <- 0

    #bim for validation
    bim_validation <- fread(paste0('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/6_HA/All_validation_chr',chr_id,'.bim'), data.table =F)
    names(bim_validation) <- LETTERS[1:6]

    #cov for validation samples
    cov_validation1 <- fread('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/1_covariate/covariate.txt', data.table =F)
    rownames(cov_validation1) <- cov_validation1$NWDID
    cov_validation1 <- cov_validation1[,-1]
    cov_validation1 <- cov_validation1[rownames(cov_validation1) %in% rownames(geno_validation),]
    cov_validation1 <- cov_validation1[,-1]

    cov_validation_pca <- fread('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/6_HA/validation_PCA/combined_indep_PCA.pca.evec', data.table = F, skip = 1)
    rownames(cov_validation_pca) <- cov_validation_pca$V1
    cov_validation_pca <- cov_validation_pca[,-1]
    cov_validation_pca <- cov_validation_pca[,-11]
    colnames(cov_validation_pca) <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')

    cov_validation <- merge(cov_validation1, cov_validation_pca, by = 'row.names')
    rownames(cov_validation) <- cov_validation$Row.names
    cov_validation <- cov_validation[,-1]
    #cov_validation <- PatchUp(cov_validation)
    cov_validation <- cov_validation[,c(1:7,9:22)] #select sex age Bcell CD4T CD8T Eos Mono NK cig1c pkyrs1c BMI PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10



    files <- list.files(paste0('/mnt/lvm_vol_2/txu/projects/mwas/fredhutch/updated.data.202408/HA.updated.analysis.202410/p12.HA.chr', chr_id, '/'), pattern =glob2rx("clean.metal.cg*.model"))
    #get overlap cpg id
    cpg_id <- str_split(files, '[.]') %>% unlist() %>% matrix(., ncol = 6, byrow = TRUE)
    cpg_id <- cpg_id[,3]
    cpg_id_overlap <- intersect(rownames(pheno_true), cpg_id[!cpg_id %in% finished_files_df$CpG])


    #for (id in cpg_id_overlap){
    selection_validation <- function(id){
        # PatchUp
        PatchUp <- function(M) {
            M <- apply(M, 2, function(x) {
                    x[is.na(x)] <- mean(x, na.rm = TRUE)
                    return(x)
                })

                return(M)
            }

        # Standardize
        Standardize <- function(M) {
            # Centralize
            M <- M - matrix(rep(colMeans(M), times = nrow(M)), nrow = nrow(M) , ncol = ncol(M), byrow = T)

            # Standardize
            M <- sweep(M, 2, sqrt(apply(M, 2, crossprod) / nrow(M)), "/")

            return(M)
        }

        print (paste0('Analysing ', id))
        #read model weight
        model <- read.table(paste0('/mnt/lvm_vol_2/txu/projects/mwas/fredhutch/updated.data.202408/HA.updated.analysis.202410/p12.HA.chr', chr_id, '/clean.metal.', id, '.1.tbl.model'), sep = ' ')
        # model <- model[,-1]
        #read model snps
        snps <- read.table(paste0('/mnt/lvm_vol_2/txu/projects/mwas/fredhutch/updated.data.202408/HA.updated.analysis.202410/p12.HA.chr', chr_id, '/clean.metal.', id, '.1.tbl.snps'), sep = ' ')
        names(snps) <- letters[2:10]
        # tmp <- str_split(snps$b, ':') %>% unlist() %>% matrix(., ncol = 4, byrow = TRUE)
        # tmp[, 1] <- gsub("chr", "", tmp[, 1])
        # tmp <- paste0(tmp[, 1], ':', tmp[, 2])
        # snps$b <- tmp
        if (nrow(snps) == nrow(model)){
            model_all <- cbind(snps, model)
            model_all <- model_all[model_all$b %in% bim_selection$B, ]

            temp <- left_join(model_all, bim_selection, by = c("b" = 'B'))

            temp <- temp[(temp$c == temp$E & temp$d == temp$F) | (temp$c == temp$F & temp$d == temp$E),]
            model_all[temp$c != temp$E, 10:109] <- model_all[temp$c != temp$E, 10:109] * -1


            geno_selection_tmp <- geno_selection[, model_all$b , drop=FALSE]
            true_y <- pheno_true[id,]


            get_r <- function(weight){
                predict_y <- geno_selection_tmp %*% matrix(weight, ncol = 1, byrow =TRUE)
                true_y_selection <- true_y[, rownames(predict_y)]
                
                true_y_selection <- merge(t(true_y_selection), cov_selection, by = 'row.names')
                rownames(true_y_selection) <- true_y_selection$Row.names
                true_y_selection <- true_y_selection[,-1]
                colnames(true_y_selection)[1] <- 'y'
    
                # Remove rows with y equal to -Inf
                true_y_selection <- true_y_selection[true_y_selection$y != -Inf, ]
                # Remove rows with missing y values
                true_y_selection <- true_y_selection[complete.cases(true_y_selection$y), ]

                true_y_selection$y <- RankNorm(true_y_selection$y)
                fit <- lm(y~.,true_y_selection)
                residule <- resid(fit)
                
                compare <- data.frame(residule)
                compare <- merge(compare, predict_y, by = 'row.names')
                colnames(compare)[3] <- 'predict_y'
                reg <- summary(lm(compare$residule ~ compare$predict_y))
                return(reg$adj.r.sq)
         
            }

            answer <- apply(model_all[, 10:109], 2, get_r)

            #best model

            which.max(answer)

            names(which.max(answer))  

            best_r2 <- answer[names(which.max(answer))]



            #validation model
            geno_validation_tmp <- geno_validation[, model_all$b, drop = FALSE]
            predict_y_validation <- geno_validation_tmp %*% matrix(model_all[,names(which.max(answer))], ncol = 1, byrow =TRUE)
            true_y_validation <- true_y[, rownames(predict_y_validation)]
            true_y_validation <- merge(t(true_y_validation), cov_validation, by = 'row.names')
            rownames(true_y_validation) <- true_y_validation$Row.names
            true_y_validation <- true_y_validation[,-1]
            colnames(true_y_validation)[1] <- 'y'


            # Remove rows with y equal to -Inf
            true_y_validation <- true_y_validation[true_y_validation$y != -Inf, ]
            # Remove rows with missing y values
            true_y_validation <- true_y_validation[complete.cases(true_y_validation$y), ]
            
            true_y_validation$y <- RankNorm(true_y_validation$y)
            fit <- lm(y~.,true_y_validation)
            residule <- resid(fit)

            compare <- data.frame(residule)
            compare <- merge(compare, predict_y_validation, by = 'row.names')
            colnames(compare)[3] <- 'predict_y_validation'

            reg_validate <- summary(lm(compare$residule ~ compare$predict_y_validation))
            res_validate <- reg_validate$adj.r.sq

            merge_res <- data.frame(answer)
            colnames(merge_res) <- 'R2_all_models'
            write.table(merge_res, paste0('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/6_HA_model_selection_validation/', id, '_all_models_R2.txt'), row.names = FALSE, sep = '\t')

            vector <- c(id, best_r2, res_validate)
            validate_res <- t(data.frame(vector)) %>% data.frame()
            
            colnames(validate_res) <- c('cpg', 'best_r2', 'validate_r2')
            validate_res$N_snp_in_models <- nrow(snps)
            validate_res$N_snp_in_selection_validation <- ncol(geno_validation_tmp)
            write.table(validate_res, paste0('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/6_HA_model_selection_validation/', id, '_best_models_R2_validation_r2.txt'), row.names = FALSE, quote =FALSE, sep = '\t')

        }else{
            print (paste0('number of snps in models and information not match ', id))
        }
    }

    # Initialize a parallel backend
    cl <- makeCluster(56)
    # Register the parallel backend
    registerDoParallel(cl)
    # Apply a function to each item in parallel
    foreach (i = cpg_id_overlap, .packages = c("BEDMatrix", "data.table", "stringr", "dplyr", "RNOmni", "foreach", "doParallel", "optparse"),
        .export = c("chr_id", "pheno_true", "geno_selection", "bim_selection", "cov_selection", "geno_validation", "bim_validation", "cov_validation")) %dopar% {
        selection_validation(i)
        gc()  # Release memory
    }
    # Stop the parallel backend
    stopCluster(cl)

}
