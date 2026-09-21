library(data.table)
library(dplyr)
library(haven)

if (!dir.exists('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/')) {
  dir.create('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/')
}



methylation_beta <- fread('/mnt/lvm_vol_2/sliu/project/MEC_former_smoker_data/raw_data/p1_mec_fs_betas_rcp.csv', data.table = F)
row.names(methylation_beta) <- methylation_beta$V1
methylation_beta <- methylation_beta[,-1]

ID_match <- read_sas('/mnt/lvm_vol_2/sliu/project/MEC_former_smoker_data/raw_data/fs_lung_sentrix_cc_202303_ids.sas7bdat')
ID_match <- data.frame(ID_match)

############################################
### 1 extract genotype for each ancestry ###
############################################
# AA
AA_df <- fread('12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/AA/cg00000103.profile', data.table = F)
nrow(AA_df) # 112
AA_sample <- filter(ID_match, ID_match$ecid%in% AA_df$IID)
AA_sample_with_methy <- filter(AA_sample, AA_sample$methyl_sampleID %in% colnames(methylation_beta))
nrow(AA_sample_with_methy) # 112

AA_sample <- fread('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/genotype_for_validation_AA.fam', data.table = F, header =F)
AA_sample <- AA_sample[AA_sample$V2 %in% AA_sample_with_methy$ecid,]
AA_sample <- select(AA_sample, V1, V2)
if (!dir.exists('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/AA/')) {
  dir.create('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/AA/')
}
write.table (AA_sample, paste0('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/AA/sample_list.txt'), quote = F, sep = '\t', col.names = F, row.names = F)
# extract sample
for (i in 1:22){
  cmd <- paste0('plink --bfile /mnt/lvm_vol_2/sliu/project/MEC_former_smoker_data/all_SNP_QC_0.99/bfile_merge/chr', i, ' --keep 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/AA/sample_list.txt --make-bed --geno 0.05 --maf 0.05 --out 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/AA/chr', i)
  system(cmd)
}
# merge 22 chromosomes
df <- data.frame(
  file = as.character()
)
for (i in 2:22){
  tmp <- paste0('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/AA/chr', i)
  df[nrow(df)+1,] <- tmp
}
write.table (df, '12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/AA/merge_list.txt', quote = F, sep = '\t', col.names = F, row.names = F)
cmd <- paste0('plink --bfile 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/AA/chr1 --merge-list 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/AA/merge_list.txt --make-bed --out 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/AA/all_chromosomes')
system(cmd)


# CA
CA_df <- fread('12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/CA/cg00000029.profile', data.table = F)
nrow(CA_df) # 283
CA_sample <- filter(ID_match, ID_match$ecid%in% CA_df$IID)
CA_sample_with_methy <- filter(CA_sample, CA_sample$methyl_sampleID %in% colnames(methylation_beta))
nrow(CA_sample_with_methy) # 283


CA_sample <- fread('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/genotype_for_validation_CA.fam', data.table = F, header =F)
CA_sample <- CA_sample[CA_sample$V2 %in% CA_sample_with_methy$ecid,]
CA_sample <- select(CA_sample, V1, V2)
if (!dir.exists('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/CA/')) {
  dir.create('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/CA/')
}
write.table (CA_sample, paste0('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/CA/sample_list.txt'), quote = F, sep = '\t', col.names = F, row.names = F)
# extract sample
for (i in 1:22){
  cmd <- paste0('plink --bfile /mnt/lvm_vol_2/sliu/project/MEC_former_smoker_data/all_SNP_QC_0.99/bfile_merge/chr', i, ' --keep 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/CA/sample_list.txt --make-bed --geno 0.05 --maf 0.05 --out 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/CA/chr', i)
  system(cmd)
}
# merge 22 chromosomes
df <- data.frame(
  file = as.character()
)
for (i in 2:22){
  tmp <- paste0('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/CA/chr', i)
  df[nrow(df)+1,] <- tmp
}
write.table (df, '12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/CA/merge_list.txt', quote = F, sep = '\t', col.names = F, row.names = F)
cmd <- paste0('plink --bfile 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/CA/chr1 --merge-list 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/CA/merge_list.txt --make-bed --out 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/CA/all_chromosomes')
system(cmd)


# EA
EA_df <- fread('12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/EA/cg00000658.profile', data.table = F)
nrow(EA_df) # 267
EA_sample <- filter(ID_match, ID_match$ecid%in% EA_df$IID)
EA_sample_with_methy <- filter(EA_sample, EA_sample$methyl_sampleID %in% colnames(methylation_beta))
nrow(EA_sample_with_methy) # 267


EA_sample <- fread('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/genotype_for_validation_EA.fam', data.table = F, header =F)
EA_sample <- EA_sample[EA_sample$V2 %in% EA_sample_with_methy$ecid,]
EA_sample <- select(EA_sample, V1, V2)
if (!dir.exists('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/EA/')) {
  dir.create('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/EA/')
}
write.table (EA_sample, paste0('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/EA/sample_list.txt'), quote = F, sep = '\t', col.names = F, row.names = F)
# extract sample
for (i in 1:22){
  cmd <- paste0('plink --bfile /mnt/lvm_vol_2/sliu/project/MEC_former_smoker_data/all_SNP_QC_0.99/bfile_merge/chr', i, ' --keep 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/EA/sample_list.txt --make-bed --geno 0.05 --maf 0.05 --out 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/EA/chr', i)
  system(cmd)
}
# merge 22 chromosomes
df <- data.frame(
  file = as.character()
)
for (i in 2:22){
  tmp <- paste0('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/EA/chr', i)
  df[nrow(df)+1,] <- tmp
}
write.table (df, '12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/EA/merge_list.txt', quote = F, sep = '\t', col.names = F, row.names = F)
cmd <- paste0('plink --bfile 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/EA/chr1 --merge-list 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/EA/merge_list.txt --make-bed --out 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/EA/all_chromosomes')
system(cmd)

# HA
HA_df <- fread('12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/HA/cg00000158.profile', data.table = F)
nrow(HA_df) # 138
HA_sample <- filter(ID_match, ID_match$ecid%in% HA_df$IID)
HA_sample_with_methy <- filter(HA_sample, HA_sample$methyl_sampleID %in% colnames(methylation_beta))
nrow(HA_sample_with_methy) # 137


HA_sample <- fread('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/genotype_for_validation_HA.fam', data.table = F, header =F)
HA_sample <- HA_sample[HA_sample$V2 %in% HA_sample_with_methy$ecid,]
HA_sample <- select(HA_sample, V1, V2)
if (!dir.exists('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/HA/')) {
  dir.create('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/HA/')
}
write.table (HA_sample, paste0('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/HA/sample_list.txt'), quote = F, sep = '\t', col.names = F, row.names = F)
# extract sample
for (i in 1:22){
  cmd <- paste0('plink --bfile /mnt/lvm_vol_2/sliu/project/MEC_former_smoker_data/all_SNP_QC_0.99/bfile_merge/chr', i, ' --keep 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/HA/sample_list.txt --make-bed --geno 0.05 --maf 0.05 --out 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/HA/chr', i)
  system(cmd)
}
# merge 22 chromosomes
df <- data.frame(
  file = as.character()
)
for (i in 2:22){
  tmp <- paste0('12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/HA/chr', i)
  df[nrow(df)+1,] <- tmp
}
write.table (df, '12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/HA/merge_list.txt', quote = F, sep = '\t', col.names = F, row.names = F)
cmd <- paste0('plink --bfile 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/HA/chr1 --merge-list 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/HA/merge_list.txt --make-bed --out 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/HA/all_chromosomes')
system(cmd)


###############################################
### 2 generate PCA matrix for each ancestry ###
###############################################
if (!dir.exists('12_validate_models_using_MEC_former_smoking/05_PCA/')) {
  dir.create('12_validate_models_using_MEC_former_smoking/05_PCA/')
}

system('python /mnt/lvm_vol_1/sliu/work/pipeline/TWAS_fusion/PCA_using_EIGENSOFT_update.py -in_gen 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/AA/all_chromosomes -Eur N -odir 12_validate_models_using_MEC_former_smoking/05_PCA/AA')
system('python /mnt/lvm_vol_1/sliu/work/pipeline/TWAS_fusion/PCA_using_EIGENSOFT_update.py -in_gen 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/CA/all_chromosomes -Eur N -odir 12_validate_models_using_MEC_former_smoking/05_PCA/CA')
system('python /mnt/lvm_vol_1/sliu/work/pipeline/TWAS_fusion/PCA_using_EIGENSOFT_update.py -in_gen 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/EA/all_chromosomes -Eur Y -odir 12_validate_models_using_MEC_former_smoking/05_PCA/EA')
system('python /mnt/lvm_vol_1/sliu/work/pipeline/TWAS_fusion/PCA_using_EIGENSOFT_update.py -in_gen 12_validate_models_using_MEC_former_smoking/04_SNP_rsq_0.99/HA/all_chromosomes -Eur N -odir 12_validate_models_using_MEC_former_smoking/05_PCA/HA')
