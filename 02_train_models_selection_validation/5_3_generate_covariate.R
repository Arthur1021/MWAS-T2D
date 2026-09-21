
library(data.table)
library(dplyr)
library(haven)

if(!dir.exists('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/06_cov')){
    dir.create('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/06_cov')
}
cov <- read.csv('/mnt/lvm_vol_2/sliu/project/MEC_former_smoker_data/raw_data/covariate_file_former_smoker.csv') 
ID_match <- read_sas('/mnt/lvm_vol_2/sliu/project/MEC_former_smoker_data/raw_data/fs_lung_sentrix_cc_202303_ids.sas7bdat')
ID_match <- data.frame(ID_match)
ID_match <- select(ID_match, methyl_sampleID, ecid)
cov <- left_join(cov, ID_match, by = c('methyl_sampleID' = 'methyl_sampleID'))
rownames(cov) <- cov$ecid
# age at blood draw = age_draw
# sex = Q1_CORR_SEX
# Body Mass Index = Q1_newqi_corr2
# smoking packyears = SPEC_packyrs
# H_BCELL = IDOL_Bcell
# H_CD4T = IDOL_CD4T
# H_CD8T = IDOL_CD8T
# H_MONO = IDOL_Mono
# H_NK = IDOL_NK
# PC1 ~ PC10: filter genotype with imputation r2 and only r2>= 0.99 will be used for generate PCA for each ancestry

cov <- select(cov, age_draw, Q1_CORR_SEX, Q1_newqi_corr2, Q1_packyrs, SPEC_packyrs, IDOL_Bcell, IDOL_CD4T, IDOL_CD8T, IDOL_Mono, IDOL_NK)

cov <- cov %>% mutate(row_name = row.names(.))



# AA
PCA_AA <- fread('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/05_PCA/AA/combined_indep_PCA.pca.evec', data.table = F)
row.names(PCA_AA) <- PCA_AA$V1
PCA_AA <- PCA_AA[1:11]
colnames(PCA_AA) <- c('row_name', 'PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5', 'PCA6', 'PCA7', 'PCA8', 'PCA9', 'PCA10')
cov_AA <- left_join(PCA_AA, cov, by = 'row_name')
row.names(cov_AA) <- cov_AA$row_name
cov_AA <- cov_AA[,-1]
cov_AA <- cbind(row.names(cov_AA), cov_AA)
colnames(cov_AA)[1] <- 'sample_ID'
write.table (cov_AA, '/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/06_cov/cov_AA.txt', sep = '\t', row.names = F, quote = F)

# CA
PCA_CA <- fread('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/05_PCA/CA/combined_indep_PCA.pca.evec', data.table = F)
row.names(PCA_CA) <- PCA_CA$V1
PCA_CA <- PCA_CA[1:11]
colnames(PCA_CA) <- c('row_name', 'PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5', 'PCA6', 'PCA7', 'PCA8', 'PCA9', 'PCA10')
cov_CA <- left_join(PCA_CA, cov, by = 'row_name')
row.names(cov_CA) <- cov_CA$row_name
cov_CA <- cov_CA[,-1]
cov_CA <- cbind(row.names(cov_CA), cov_CA)
colnames(cov_CA)[1] <- 'sample_ID'
write.table (cov_CA, '/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/06_cov/cov_CA.txt', sep = '\t', row.names = F, quote = F)

# EA
PCA_EA <- fread('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/05_PCA/EA/combined_indep_PCA.pca.evec', data.table = F)
row.names(PCA_EA) <- PCA_EA$V1
PCA_EA <- PCA_EA[1:11]
colnames(PCA_EA) <- c('row_name', 'PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5', 'PCA6', 'PCA7', 'PCA8', 'PCA9', 'PCA10')
cov_EA <- left_join(PCA_EA, cov, by = 'row_name')
row.names(cov_EA) <- cov_EA$row_name
cov_EA <- cov_EA[,-1]
cov_EA <- cbind(row.names(cov_EA), cov_EA)
colnames(cov_EA)[1] <- 'sample_ID'
write.table (cov_EA, '/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/06_cov/cov_EA.txt', sep = '\t', row.names = F, quote = F)

# HA
PCA_HA <- fread('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/05_PCA/HA/combined_indep_PCA.pca.evec', data.table = F)
row.names(PCA_HA) <- PCA_HA$V1
PCA_HA <- PCA_HA[1:11]
colnames(PCA_HA) <- c('row_name', 'PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5', 'PCA6', 'PCA7', 'PCA8', 'PCA9', 'PCA10')
cov_HA <- left_join(PCA_HA, cov, by = 'row_name')
row.names(cov_HA) <- cov_HA$row_name
cov_HA <- cov_HA[,-1]
cov_HA <- cbind(row.names(cov_HA), cov_HA)
colnames(cov_HA)[1] <- 'sample_ID'
write.table (cov_HA, '/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/06_cov/cov_HA.txt', sep = '\t', row.names = F, quote = F)

