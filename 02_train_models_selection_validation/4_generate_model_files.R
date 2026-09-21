library(data.table)
library(dplyr)
library(tidyverse)
library(foreach)
library(doParallel)


races <- c('AA', 'CA', 'EA', 'HA')
for (race in races){
    model_performance <- fread(paste0('7_model_summary/model_summary_', race, '.txt'), data.table = F)
    model_performance <- filter(model_performance, best_r2 >= 0.01, validate_r2 >= 0.01)
    Chr_pos_rsID_converter <- fread(paste0('8_SNP_converter/', race, '_SNP_convert.txt'), data.table = F, header = F)
    colnames(Chr_pos_rsID_converter) <- c('chr_pos', 'rsID')
    if (!file.exists(paste0("9_generate_model_files/",race))) {
        dir.create(paste0("9_generate_model_files/",race), recursive = TRUE)
    }


    ##1.3 generate models RDat files
    make_model_RDat_file <- function(i, Chr_pos_rsID_converter){
        if (race == 'AA'){
            prefix <- '3_AA'
        }else if (race == 'CA'){
            prefix <- '4_CA'
        }else if (race == 'EA'){
            prefix <- '5_EA'
        }else if (race == 'HA'){
            prefix <- '6_HA'
        }
        selection_validation_file <- paste0(prefix, '_model_selection_validation/', model_performance[i,]$cpg, '_all_models_R2.txt')
        selection_validation_value <- fread(selection_validation_file, data.table = F)
        best_n <- which.max(selection_validation_value$R2_all_models)

        if (race == 'AA'){
            weight_file <- Sys.glob(paste0('/mnt/lvm_vol_2/txu/projects/mwas/fredhutch/updated.data.202408/AA.updated.analysis.202410/p10.AA.chr*/clean.metal.', model_performance[i,]$cpg, '.1.tbl.model' )) 
        }else if (race == 'CA' ){
            weight_file <- Sys.glob(paste0('/mnt/lvm_vol_2/txu/projects/mwas/fredhutch/updated.data.202408/JA.updated.analysis.202410/p11.JA.chr*/clean.metal.', model_performance[i,]$cpg, '.1.tbl.model' ))
        }else if (race == 'EA' ){
            weight_file <- Sys.glob(paste0('/mnt/lvm_vol_2/txu/projects/mwas/fredhutch/updated.data.202408/EA.updated.analysis.202410/p14.EA.chr*/clean.metal.', model_performance[i,]$cpg, '.1.tbl.model' ))
        }else if (race == 'HA' ){
            weight_file <- Sys.glob(paste0('/mnt/lvm_vol_2/txu/projects/mwas/fredhutch/updated.data.202408/HA.updated.analysis.202410/p12.HA.chr*/clean.metal.', model_performance[i,]$cpg, '.1.tbl.model' ))
        }
        WEIGHT <- read.table(weight_file, sep = ' ')
        if (race == 'AA'){
            snp_file <- Sys.glob(paste0('/mnt/lvm_vol_2/txu/projects/mwas/fredhutch/updated.data.202408/AA.updated.analysis.202410/p10.AA.chr*/clean.metal.', model_performance[i,]$cpg, '.1.tbl.snps' )) 
        }else if(race == 'CA' ){
            snp_file <- Sys.glob(paste0('/mnt/lvm_vol_2/txu/projects/mwas/fredhutch/updated.data.202408/JA.updated.analysis.202410/p11.JA.chr*/clean.metal.', model_performance[i,]$cpg, '.1.tbl.snps' ))
        }else if(race == 'EA' ){
            snp_file <- Sys.glob(paste0('/mnt/lvm_vol_2/txu/projects/mwas/fredhutch/updated.data.202408/EA.updated.analysis.202410/p14.EA.chr*/clean.metal.', model_performance[i,]$cpg, '.1.tbl.snps' ))
        }else if(race == 'HA' ){
            snp_file <- Sys.glob(paste0('/mnt/lvm_vol_2/txu/projects/mwas/fredhutch/updated.data.202408/HA.updated.analysis.202410/p12.HA.chr*/clean.metal.', model_performance[i,]$cpg, '.1.tbl.snps' ))
        }
        SNP <- read.table(snp_file, sep = ' ')
        if (nrow(SNP) == nrow(WEIGHT)){
            SNP_WEIGHT <- cbind(SNP[,1:3], WEIGHT[,best_n])
            colnames(SNP_WEIGHT)[4] <- 'weight'
            # #split
            # SNP_WEIGHT <- SNP_WEIGHT %>% separate(SNP, into = c("Chr", "Pos", "N1", "N2"), sep = ":")

            SNP_WEIGHT$chr_pos <- sub("chr(\\d+):(\\d+):.*", "\\1:\\2", SNP_WEIGHT$SNP)
            SNP_WEIGHT <- filter(SNP_WEIGHT, !weight == 0)
            if (nrow(SNP_WEIGHT) > 0){

                SNP_WEIGHT <- merge(SNP_WEIGHT, Chr_pos_rsID_converter, by = 'chr_pos', all.x = T)
                SNP_WEIGHT <- SNP_WEIGHT %>% separate(chr_pos, into = c("V1", "V4"), sep = ":")
                SNP_WEIGHT$V2 <- SNP_WEIGHT$rsID
                SNP_WEIGHT$V3 <- 0
                SNP_WEIGHT$V5 <- SNP_WEIGHT$A1
                SNP_WEIGHT$V6 <- SNP_WEIGHT$A2

                
                #cv.performance
                SUMMIT <- c(round(model_performance[i,]$validate_r2, 2), '0')
                cv.performance <- data.frame(SUMMIT)
                rownames(cv.performance) <- c('rsq', 'pval')
                #snps
                SNP_WEIGHT <- SNP_WEIGHT %>% distinct(rsID, .keep_all = TRUE) # remove duplicate SNP with same Chr Pos, only retain first one
                wgt.matrix <- data.frame(SNP_WEIGHT$weight)
                row.names(wgt.matrix) <- SNP_WEIGHT$V2
                colnames(wgt.matrix) <- 'SUMMIT'

                #snp
                snps <- select(SNP_WEIGHT, V1, V2, V3, V4, V5, V6)
                
                #save
                save( wgt.matrix , snps , cv.performance , file = paste0( '9_generate_model_files/',race, '/', model_performance[i,]$cpg, ".wgt.RDat" ) )
            }
        }
        
    }

    # Initialize a parallel backend
    cl <- makeCluster(56)
    # Register the parallel backend
    registerDoParallel(cl)
    # Apply a function to each item in parallel
    foreach (i = 1:nrow(model_performance), .packages = c( "data.table",  "dplyr",  "foreach", "doParallel", "tidyverse")) %dopar% {
        make_model_RDat_file(i,Chr_pos_rsID_converter)
    }
    # Stop the parallel backend
    stopCluster(cl)
}

##################################
### generate filtered pos file ###
##################################

# get all CpG site information Chr pos 38 version
CpG_position <- data.frame()
CpG_position_files <- Sys.glob('/mnt/lvm_vol_2/txu/projects/mwas/fredhutch/model.building.AA/chr*.CpG.tab')
for (file in CpG_position_files){
    df_tmp <- fread(file, data.table = F, header = F)
    CpG_position <- rbind(CpG_position, df_tmp)
}


# get all CpG site with models
races <- c('AA', 'CA', 'EA', 'HA')
for (race in races){

    files <- Sys.glob(paste0('9_generate_model_files/', race, '/*.RDat'))
    model_file_list <- data.frame(files)
    model_file_list$WGT <- sapply(strsplit(model_file_list$files, "/"), "[", 3)
    model_file_list$ID <- sub("\\..*", "", model_file_list$WGT)

    #merge WGT and Chr pos
    model_file_list_Chr_pos <- merge(model_file_list, CpG_position, by.x = 'ID', by.y = 'V1', all.x = T)
    model_file_list_Chr_pos$P0 <- model_file_list_Chr_pos$V3 - 1
    model_file_list_Chr_pos$P1 <- model_file_list_Chr_pos$V3 + 1

    print(dim(model_file_list_Chr_pos)) #[1] 227346      7

    ## Filtering probes based on Pidsley et al 2016 documentation##
    # Removing cross-reactive probes
    crossreactives <- read.csv ('/mnt/lvm_vol_2/sliu/project/meQTL_in_5_races/04_filter_cpg_site/eGTEx_mQTLs_eQTLs_GWAS/DNAm_processing/crossreactives_pidsley.csv')
    model_file_list_Chr_pos <- model_file_list_Chr_pos[!(model_file_list_Chr_pos$ID %in% crossreactives$X), ]
    dim(model_file_list_Chr_pos) #[1] 226342      7

    # SNPs within single base extension SBE or in CpG
    snp_cpg <- read.csv ('/mnt/lvm_vol_2/sliu/project/meQTL_in_5_races/04_filter_cpg_site/eGTEx_mQTLs_eQTLs_GWAS/DNAm_processing/snp_cpg_pidsley.csv')
    model_file_list_Chr_pos <- model_file_list_Chr_pos[!(model_file_list_Chr_pos$ID %in% snp_cpg$PROBE),]
    dim(model_file_list_Chr_pos) #[1] 218259      7

    # Remove RS probes
    snp_sbe <- read.csv ('/mnt/lvm_vol_2/sliu/project/meQTL_in_5_races/04_filter_cpg_site/eGTEx_mQTLs_eQTLs_GWAS/DNAm_processing/snp_sbe_pidsley.csv')
    model_file_list_Chr_pos <- model_file_list_Chr_pos[!(model_file_list_Chr_pos$ID %in% snp_sbe$PROBE),]
    print (dim(model_file_list_Chr_pos))

    model_file_list_Chr_pos$PANEL <- race
    model_file_list_Chr_pos$N <- 'NA'
    model_file_list_Chr_pos <- select(model_file_list_Chr_pos, PANEL, WGT, ID, V2, P0, P1, N)
    colnames(model_file_list_Chr_pos)[4] <- 'CHR'
    model_file_list_Chr_pos$CHR <- as.numeric(model_file_list_Chr_pos$CHR)
    model_file_list_Chr_pos <- model_file_list_Chr_pos[order(model_file_list_Chr_pos$CHR, model_file_list_Chr_pos$P0),]
    write.table(model_file_list_Chr_pos, file = paste0("9_generate_model_files/", race, "_methylation_models_hg38.pos"), quote = F, sep = "\t", row.names = F)

}
