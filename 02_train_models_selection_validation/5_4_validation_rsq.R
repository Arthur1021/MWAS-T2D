library(data.table)
library(dplyr)
library(haven)
library(lumi)
library(stringr)
library(RNOmni)
library(doParallel)
library(foreach)

if(!dir.exists('12_validate_models_using_MEC_former_smoking/07_validation_rsq_using_MEC')){
    dir.create('12_validate_models_using_MEC_former_smoking/07_validation_rsq_using_MEC')
}


ID_match <- read_sas('/mnt/lvm_vol_2/sliu/project/MEC_former_smoker_data/raw_data/fs_lung_sentrix_cc_202303_ids.sas7bdat')
ID_match <- data.frame(ID_match)
ID_match <- select(ID_match, methyl_sampleID, ecid)

methy_beta <- fread('/mnt/lvm_vol_2/sliu/project/MEC_former_smoker_data/raw_data/p1_mec_fs_betas_rcp.csv', data.table = F)
row.names(methy_beta) <- methy_beta$V1
methy_beta <- methy_beta[, -1]
#replace 204229100049_R04C01 with ECxxxxxx ID
ID <- data.frame(
    sample = as.character(colnames(methy_beta))
)
ID <- left_join(ID, ID_match, by = c('sample' = 'methyl_sampleID'))
missing_rows <- is.na(ID$ecid)
ID$ecid[missing_rows] <- ID$sample[missing_rows]
colnames(methy_beta) <- ID$ecid

# convert beta to M value

message('Convert beta value to M value')
m_value <- beta2m(methy_beta)
m_value <- t(m_value)
m_value <- data.frame(m_value)


message('merge methylatio M value and covariates')
# AA methylation
cov_AA <- fread('12_validate_models_using_MEC_former_smoking/06_cov/cov_AA.txt', data.table = F)
phe_AA <- m_value[rownames(m_value) %in% cov_AA$sample_ID,]
AA_predicted_methylation <- Sys.glob('12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/AA/*profile')

# CA methylation
cov_CA <- fread('12_validate_models_using_MEC_former_smoking/06_cov/cov_CA.txt', data.table = F)
phe_CA <- m_value[rownames(m_value) %in% cov_CA$sample_ID,]
CA_predicted_methylation <- Sys.glob('12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/CA/*profile')

# EA methylation
cov_EA <- fread('12_validate_models_using_MEC_former_smoking/06_cov/cov_EA.txt', data.table = F)
phe_EA <- m_value[rownames(m_value) %in% cov_EA$sample_ID,]
EA_predicted_methylation <- Sys.glob('12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/EA/*profile')

# HA methylation
cov_HA <- fread('12_validate_models_using_MEC_former_smoking/06_cov/cov_HA.txt', data.table = F)
phe_HA <- m_value[rownames(m_value) %in% cov_HA$sample_ID,]
HA_predicted_methylation <- Sys.glob('12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/HA/*profile')




# Set the number of cores to utilize
num_cores <- 56  # Adjust the number of cores as per your system

# Register parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)





# Create a function to process each file in parallel
process_file <- function(file, pheno, cov) {
  id <- first(str_split(last(str_split(file, '/')[[1]]), '[.]')[[1]])
  if (id %in% colnames(pheno)){
  
    phe_extract <- data.frame(cbind(row.names(pheno), pheno[, id]))
    if (sum(!is.na(phe_extract$X2))/nrow(phe_extract) < 0.5 ){
      return(NULL)
    }

    colnames(phe_extract) <- c('sample', 'pheno')
    
    # Merge phenotype and covariate
    phe_extract <- left_join(phe_extract, cov, by = c('sample' = 'sample_ID'))
    rownames(phe_extract) <- phe_extract$sample
    # Remove rows with NA values
    phe_extract <- phe_extract[complete.cases(phe_extract), ]

    #phe_extract$pheno <- RankNorm(phe_extract$pheno)  # Ranknorm pheno
    phe_extract <- phe_extract[, -1]
    
    # Calculate residuals
    fit <- lm(pheno ~ ., phe_extract)
    residule <- resid(fit)
    residule <- data.frame(residule)
    residule$sample <- row.names(residule)

    # Compare residuals and predicted value
    predict_y_validation <- fread(file, data.table = FALSE)
    row.names(predict_y_validation) <- predict_y_validation$IID
    predict_y_validation <- predict_y_validation[row.names(pheno), ]
    predict_y_validation <- select(predict_y_validation, IID, SCORE)
    #merge
    residule <- left_join(residule, predict_y_validation, by = c('sample'= 'IID'))

    residule <- residule[complete.cases(residule), ]

    
    reg_validate <- summary(lm(residule$residule ~ residule$SCORE))
    rsq <- reg_validate$adj.r.sq
    
    # Return the result
    return(data.frame(CpG = id, rsq = rsq))
  }
}


# AA validation rsq using MEC
message('Validation methylation prediction models in AA ancestry using MEC former smokers')
# Iterate over files in parallel and generate results_df directly
results_df <- foreach(file = AA_predicted_methylation, .packages = c("stringr", "dplyr", "RNOmni", "data.table")) %dopar% {
  process_file(file, phe_AA, cov_AA)
}
# Combine the individual data frames using bind_rows
results_df <- bind_rows(results_df)
fwrite(results_df, '12_validate_models_using_MEC_former_smoking/07_validation_rsq_using_MEC/validation_using_MEC_AA.txt', sep = '\t', row.names = FALSE, quote = FALSE)



# CA validation rsq using MEC
message('Validation methylation prediction models in CA ancestry using MEC former smokers')
# Iterate over files in parallel and generate results_df directly
results_df <- foreach(file = CA_predicted_methylation, .packages = c("stringr", "dplyr", "RNOmni", "data.table")) %dopar% {
  process_file(file, phe_CA, cov_CA)
}
# Combine the individual data frames using bind_rows
results_df <- bind_rows(results_df)
fwrite(results_df, '12_validate_models_using_MEC_former_smoking/07_validation_rsq_using_MEC/validation_using_MEC_CA.txt', sep = '\t', row.names = FALSE, quote = FALSE)




# EA validation rsq using MEC
message('Validation methylation prediction models in EA ancestry using MEC former smokers')
# Iterate over files in parallel and generate results_df directly
results_df <- foreach(file = EA_predicted_methylation, .packages = c("stringr", "dplyr", "RNOmni", "data.table")) %dopar% {
  process_file(file, phe_EA, cov_EA)
}
# Combine the individual data frames using bind_rows
results_df <- bind_rows(results_df)
fwrite(results_df, '12_validate_models_using_MEC_former_smoking/07_validation_rsq_using_MEC/validation_using_MEC_EA.txt', sep = '\t', row.names = FALSE, quote = FALSE)





# HA validation rsq using MEC
message('Validation methylation prediction models in HA ancestry using MEC former smokers')
# Iterate over files in parallel and generate results_df directly
results_df <- foreach(file = HA_predicted_methylation, .packages = c("stringr", "dplyr", "RNOmni", "data.table")) %dopar% {
  process_file(file, phe_HA, cov_HA)
}
# Combine the individual data frames using bind_rows
results_df <- bind_rows(results_df)
fwrite(results_df, '12_validate_models_using_MEC_former_smoking/07_validation_rsq_using_MEC/validation_using_MEC_HA.txt', sep = '\t', row.names = FALSE, quote = FALSE)





# Stop the parallel backend
stopCluster(cl)

