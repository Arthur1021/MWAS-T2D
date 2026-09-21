# 1. get all SNP in models

library(data.table)
library(doParallel)
cl <- makeCluster(56)
registerDoParallel(cl)
# AA
dir.create('12_validate_models_using_MEC_former_smoking/01_SNP_from_models/')
dir.create('12_validate_models_using_MEC_former_smoking/01_SNP_from_models/AA/')
AA <- fread('9_generate_model_files/AA_methylation_models_hg38.pos', data.table = F)
foreach(w = 1:nrow(AA), .packages = c("data.table")) %dopar% {
  # Load the data
  load(paste0('9_generate_model_files/AA/', AA[w, 2]))
  
  # Write the table to file
  write.table(snps, paste0('12_validate_models_using_MEC_former_smoking/01_SNP_from_models/AA/', AA[w, 3], '_SNP.txt'), sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
}

# CA
dir.create('12_validate_models_using_MEC_former_smoking/01_SNP_from_models/')
dir.create('12_validate_models_using_MEC_former_smoking/01_SNP_from_models/CA/')
CA <- fread('9_generate_model_files/CA_methylation_models_hg38.pos', data.table = F)
foreach(w = 1:nrow(CA), .packages = c("data.table")) %dopar% {
  # Load the data
  load(paste0('9_generate_model_files/CA/', CA[w, 2]))
  
  # Write the table to file
  write.table(snps, paste0('12_validate_models_using_MEC_former_smoking/01_SNP_from_models/CA/', CA[w, 3], '_SNP.txt'), sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
}

# EA
dir.create('12_validate_models_using_MEC_former_smoking/01_SNP_from_models/')
dir.create('12_validate_models_using_MEC_former_smoking/01_SNP_from_models/EA/')
EA <- fread('9_generate_model_files/EA_methylation_models_hg38.pos', data.table = F)
foreach(w = 1:nrow(EA), .packages = c("data.table")) %dopar% {
  # Load the data
  load(paste0('9_generate_model_files/EA/', EA[w, 2]))
  
  # Write the table to file
  write.table(snps, paste0('12_validate_models_using_MEC_former_smoking/01_SNP_from_models/EA/', EA[w, 3], '_SNP.txt'), sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
}

# HA
dir.create('12_validate_models_using_MEC_former_smoking/01_SNP_from_models/')
dir.create('12_validate_models_using_MEC_former_smoking/01_SNP_from_models/HA/')
HA <- fread('9_generate_model_files/HA_methylation_models_hg38.pos', data.table = F)
foreach(w = 1:nrow(HA), .packages = c("data.table")) %dopar% {
  # Load the data
  load(paste0('9_generate_model_files/HA/', HA[w, 2]))
  
  # Write the table to file
  write.table(snps, paste0('12_validate_models_using_MEC_former_smoking/01_SNP_from_models/HA/', HA[w, 3], '_SNP.txt'), sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
}

# Stop the parallel processing
stopCluster(cl)


# 2. combine SNP used in models and get unique SNPs
# 2.1 AA
folder_path <- "12_validate_models_using_MEC_former_smoking/01_SNP_from_models/AA/"
file_names <- list.files(folder_path, full.names = TRUE)
AA_SNP <- rbindlist(lapply(file_names, fread), fill = TRUE)
AA_SNP <- AA_SNP[order(AA_SNP$V1, AA_SNP$V4), ]
AA_SNP <- AA_SNP[!duplicated(AA_SNP), ] 
nrow(unique(AA_SNP[,c(1,4)])) #[1] 4636886
length(unique(AA_SNP$V2)) # 4636886
nrow(AA_SNP) # 4638348
write.table(AA_SNP, '12_validate_models_using_MEC_former_smoking/01_SNP_from_models/SNPs_in_models_AA.txt', sep ='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

# 2.1 CA
folder_path <- "12_validate_models_using_MEC_former_smoking/01_SNP_from_models/CA/"
file_names <- list.files(folder_path, full.names = TRUE)
CA_SNP <- rbindlist(lapply(file_names, fread), fill = TRUE)
CA_SNP <- CA_SNP[order(CA_SNP$V1, CA_SNP$V4), ]
CA_SNP <- CA_SNP[!duplicated(CA_SNP), ] 
nrow(unique(CA_SNP[,c(1,4)])) #[1] 966581
length(unique(CA_SNP$V2)) # 966581
nrow(CA_SNP) # 966791
write.table(CA_SNP, '12_validate_models_using_MEC_former_smoking/01_SNP_from_models/SNPs_in_models_CA.txt', sep ='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

# 2.3 EA
folder_path <- "12_validate_models_using_MEC_former_smoking/01_SNP_from_models/EA/"
file_names <- list.files(folder_path, full.names = TRUE)
EA_SNP <- rbindlist(lapply(file_names, fread), fill = TRUE)
EA_SNP <- EA_SNP[order(EA_SNP$V1, EA_SNP$V4), ]
EA_SNP <- EA_SNP[!duplicated(EA_SNP), ] 
nrow(unique(EA_SNP[,c(1,4)])) #[1] 1692705
length(unique(EA_SNP$V2)) # 1692705
nrow(EA_SNP) # 1692830
write.table(EA_SNP, '12_validate_models_using_MEC_former_smoking/01_SNP_from_models/SNPs_in_models_EA.txt', sep ='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

# 2.4 HA
folder_path <- "12_validate_models_using_MEC_former_smoking/01_SNP_from_models/HA/"
file_names <- list.files(folder_path, full.names = TRUE)
HA_SNP <- rbindlist(lapply(file_names, fread), fill = TRUE)
HA_SNP <- HA_SNP[order(HA_SNP$V1, HA_SNP$V4), ]
HA_SNP <- HA_SNP[!duplicated(HA_SNP), ] 
nrow(unique(HA_SNP[,c(1,4)])) #[1] 2282065
length(unique(HA_SNP$V2)) # 2282065
nrow(HA_SNP) # 2282223
write.table(HA_SNP, '12_validate_models_using_MEC_former_smoking/01_SNP_from_models/SNPs_in_models_HA.txt', sep ='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)




# 3. get sample ID in each ancestry
library(haven)
library(dplyr)
sample_ID <- read_sas('/mnt/lvm_vol_2/sliu/project/MEC_former_smoker_data/raw_data/fs_lung_sentrix_cc_202303_ids.sas7bdat')
sample_ID <- data.frame(sample_ID)
race <- read.csv('/mnt/lvm_vol_2/sliu/project/MEC_former_smoker_data/raw_data/covariate_file_former_smoker.csv')
race <- select(race, methyl_sampleID, Q1_eth)
sample_ID_race <- merge(sample_ID, race, by = 'methyl_sampleID')
CA_sample <- filter(sample_ID_race, Q1_eth == 'J')
LA_sample <- filter(sample_ID_race, Q1_eth == 'L')
BA_sample <- filter(sample_ID_race, Q1_eth == 'B')
WA_sample <- filter(sample_ID_race, Q1_eth == 'W')
HA_sample <- filter(sample_ID_race, Q1_eth == 'H')


# 4. extract SNP for each ancestry
dir.create('12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/')
genotype <- data.frame()
for (chr in 1:22){
    file <- paste0('/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/12_validate_models_using_MEC_former_smoking/01_QC/chr', chr, '.vcf')
    message(paste0('Analysing ', file))
    df <- fread(file, data.table = F)
    genotype <- rbind(genotype, df)
}

# 4.1 AA
AA_SNP$ID_1 <- paste0('chr',AA_SNP$V1, ":", AA_SNP$V4, ":", AA_SNP$V5, ":", AA_SNP$V6)
AA_SNP$ID_2 <- paste0('chr',AA_SNP$V1, ":", AA_SNP$V4, ":", AA_SNP$V6, ":", AA_SNP$V5)

AA_genotype_filter <- genotype[ (genotype$ID %in% AA_SNP$ID_1 | genotype$ID %in% AA_SNP$ID_2),]
dim(AA_genotype_filter) #4638213     4638213/4638348 = 0.9999709

# find which columns belong to AA sample
matching_columns <- colnames(AA_genotype_filter)[sapply(strsplit(colnames(AA_genotype_filter), "_"), function(x) x[2]) %in% BA_sample$ecid]
AA_genotype_filter_sample <- select(AA_genotype_filter,all_of(matching_columns) )

# count genotype rate for each sample
genotype_rate <- data.frame(
  Sample = as.character(),
  rate = as.numeric()
)
for (i in 1:ncol(AA_genotype_filter_sample)){
  sample <- colnames(AA_genotype_filter_sample)[i]
  rate <- sum(!AA_genotype_filter_sample[i] == './.')/nrow(AA_genotype_filter_sample[i])
  genotype_rate[i,1] <- sample
  genotype_rate[i,2] <- rate
}
# Split Sample column
genotype_rate$Sample_Name <- sapply(strsplit(genotype_rate$Sample, "_"), function(x) x[2])
# Group by Sample_Name and select the sample with the highest rate
highest_rate_samples <- genotype_rate %>%
  group_by(Sample_Name) %>%
  filter(rate == max(rate)) %>%
  ungroup()

# extract sample based on highest rate
AA_genotype_filter_sample <- select(AA_genotype_filter_sample,all_of(highest_rate_samples$Sample) )

# change chr:pos to rsID
SNP_information <- AA_genotype_filter[,1:9]
SNP_information <- SNP_information %>% mutate(CHR_POS = paste0(`#CHROM`, ':', POS))
AA_SNP <- AA_SNP %>% mutate(CHR_POS = paste0('chr', V1, ':', V4))
AA_SNP <- select(AA_SNP, CHR_POS, V2)
AA_SNP <- AA_SNP[!duplicated(AA_SNP), ] 
SNP_information <- left_join(SNP_information, AA_SNP, by = 'CHR_POS')
SNP_information <- select(SNP_information, '#CHROM', POS, V2, REF, ALT, QUAL, FILTER, INFO, FORMAT)
colnames(SNP_information)[3] <- 'ID'
AA_genotype_filter_sample <- cbind(SNP_information, AA_genotype_filter_sample)

# write out result
write.table(AA_genotype_filter_sample, '12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/genotype_for_validation_AA.vcf', sep = '\t', quote = F, row.names = F)




# 4.2 CA
CA_SNP <- CA_SNP %>%
  mutate(ID_1 = paste0('chr', V1, ":", V4, ":", V5, ":", V6),
         ID_2 = paste0('chr', V1, ":", V4, ":", V6, ":", V5))



CA_genotype_filter <- genotype[genotype$ID %in% c(CA_SNP$ID_1, CA_SNP$ID_2), ]

dim(CA_genotype_filter) #[1] 883197   1137 # 883197/966791=0.9135346

# find which columns belong to CA sample
matching_columns <- colnames(CA_genotype_filter)[sapply(strsplit(colnames(CA_genotype_filter), "_"), function(x) x[2]) %in% CA_sample$ecid]
CA_genotype_filter_sample <- select(CA_genotype_filter,all_of(matching_columns) )

# count genotype rate for each sample
genotype_rate <- data.frame(
  Sample = as.character(),
  rate = as.numeric()
)
for (i in 1:ncol(CA_genotype_filter_sample)){
  sample <- colnames(CA_genotype_filter_sample)[i]
  rate <- sum(!CA_genotype_filter_sample[i] == './.')/nrow(CA_genotype_filter_sample[i])
  genotype_rate[i,1] <- sample
  genotype_rate[i,2] <- rate
}
# Split Sample column
genotype_rate$Sample_Name <- sapply(strsplit(genotype_rate$Sample, "_"), function(x) x[2])
# Group by Sample_Name and select the sample with the highest rate
highest_rate_samples <- genotype_rate %>%
  group_by(Sample_Name) %>%
  filter(rate == max(rate)) %>%
  ungroup()

# extract sample based on highest rate
CA_genotype_filter_sample <- select(CA_genotype_filter_sample,all_of(highest_rate_samples$Sample) )

# change chr:pos to rsID
SNP_information <- CA_genotype_filter[,1:9]
SNP_information <- SNP_information %>% mutate(CHR_POS = paste0(`#CHROM`, ':', POS))
CA_SNP <- CA_SNP %>% mutate(CHR_POS = paste0('chr', V1, ':', V4))
CA_SNP <- select(CA_SNP, CHR_POS, V2)
CA_SNP <- CA_SNP[!duplicated(CA_SNP), ] 
SNP_information <- left_join(SNP_information, CA_SNP, by = 'CHR_POS')
SNP_information <- select(SNP_information, '#CHROM', POS, V2, REF, ALT, QUAL, FILTER, INFO, FORMAT)
colnames(SNP_information)[3] <- 'ID'
CA_genotype_filter_sample <- cbind(SNP_information, CA_genotype_filter_sample)

# write out result
write.table(CA_genotype_filter_sample, '12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/genotype_for_validation_CA.vcf', sep = '\t', quote = F, row.names = F)




# 4.3 EA
EA_SNP <- EA_SNP %>%
  mutate(ID_1 = paste0('chr', V1, ":", V4, ":", V5, ":", V6),
         ID_2 = paste0('chr', V1, ":", V4, ":", V6, ":", V5))



EA_genotype_filter <- genotype[genotype$ID %in% c(EA_SNP$ID_1, EA_SNP$ID_2), ]

dim(EA_genotype_filter) #[1] 1692783    1137     #1692783/1692830=0.9999722

# find which columns belong to EA sample
matching_columns <- colnames(EA_genotype_filter)[sapply(strsplit(colnames(EA_genotype_filter), "_"), function(x) x[2]) %in% WA_sample$ecid]
EA_genotype_filter_sample <- select(EA_genotype_filter,all_of(matching_columns) )

# count genotype rate for each sample
genotype_rate <- data.frame(
  Sample = as.character(),
  rate = as.numeric()
)
for (i in 1:ncol(EA_genotype_filter_sample)){
  sample <- colnames(EA_genotype_filter_sample)[i]
  rate <- sum(!EA_genotype_filter_sample[i] == './.')/nrow(EA_genotype_filter_sample[i])
  genotype_rate[i,1] <- sample
  genotype_rate[i,2] <- rate
}
# Split Sample column
genotype_rate$Sample_Name <- sapply(strsplit(genotype_rate$Sample, "_"), function(x) x[2])
# Group by Sample_Name and select the sample with the highest rate
highest_rate_samples <- genotype_rate %>%
  group_by(Sample_Name) %>%
  filter(rate == max(rate)) %>%
  ungroup()

# extract sample based on highest rate
EA_genotype_filter_sample <- select(EA_genotype_filter_sample,all_of(highest_rate_samples$Sample) )

# change chr:pos to rsID
SNP_information <- EA_genotype_filter[,1:9]
SNP_information <- SNP_information %>% mutate(CHR_POS = paste0(`#CHROM`, ':', POS))
EA_SNP <- EA_SNP %>% mutate(CHR_POS = paste0('chr', V1, ':', V4))
EA_SNP <- select(EA_SNP, CHR_POS, V2)
EA_SNP <- EA_SNP[!duplicated(EA_SNP), ] 
SNP_information <- left_join(SNP_information, EA_SNP, by = 'CHR_POS')
SNP_information <- select(SNP_information, '#CHROM', POS, V2, REF, ALT, QUAL, FILTER, INFO, FORMAT)
colnames(SNP_information)[3] <- 'ID'
EA_genotype_filter_sample <- cbind(SNP_information, EA_genotype_filter_sample)

# write out result
write.table(EA_genotype_filter_sample, '12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/genotype_for_validation_EA.vcf', sep = '\t', quote = F, row.names = F)



# 4.4 HA
HA_SNP <- HA_SNP %>%
  mutate(ID_1 = paste0('chr', V1, ":", V4, ":", V5, ":", V6),
         ID_2 = paste0('chr', V1, ":", V4, ":", V6, ":", V5))



HA_genotype_filter <- genotype[genotype$ID %in% c(HA_SNP$ID_1, HA_SNP$ID_2), ]

dim(HA_genotype_filter) #[1] 2282028    1137 # 2282028/2282223=0.9999146

# find which columns belong to HA sample
matching_columns <- colnames(HA_genotype_filter)[sapply(strsplit(colnames(HA_genotype_filter), "_"), function(x) x[2]) %in% LA_sample$ecid]
HA_genotype_filter_sample <- select(HA_genotype_filter,all_of(matching_columns) )

# count genotype rate for each sample
genotype_rate <- data.frame(
  Sample = as.character(),
  rate = as.numeric()
)
for (i in 1:ncol(HA_genotype_filter_sample)){
  sample <- colnames(HA_genotype_filter_sample)[i]
  rate <- sum(!HA_genotype_filter_sample[i] == './.')/nrow(HA_genotype_filter_sample[i])
  genotype_rate[i,1] <- sample
  genotype_rate[i,2] <- rate
}
# Split Sample column
genotype_rate$Sample_Name <- sapply(strsplit(genotype_rate$Sample, "_"), function(x) x[2])
# Group by Sample_Name and select the sample with the highest rate
highest_rate_samples <- genotype_rate %>%
  group_by(Sample_Name) %>%
  filter(rate == max(rate)) %>%
  ungroup()

# extract sample based on highest rate
HA_genotype_filter_sample <- select(HA_genotype_filter_sample,all_of(highest_rate_samples$Sample) )

# change chr:pos to rsID
SNP_information <- HA_genotype_filter[,1:9]
SNP_information <- SNP_information %>% mutate(CHR_POS = paste0(`#CHROM`, ':', POS))
HA_SNP <- HA_SNP %>% mutate(CHR_POS = paste0('chr', V1, ':', V4))
HA_SNP <- select(HA_SNP, CHR_POS, V2)
HA_SNP <- HA_SNP[!duplicated(HA_SNP), ] 
SNP_information <- left_join(SNP_information, HA_SNP, by = 'CHR_POS')
SNP_information <- select(SNP_information, '#CHROM', POS, V2, REF, ALT, QUAL, FILTER, INFO, FORMAT)
colnames(SNP_information)[3] <- 'ID'
HA_genotype_filter_sample <- cbind(SNP_information, HA_genotype_filter_sample) # 


# write out result
write.table(HA_genotype_filter_sample, '12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/genotype_for_validation_HA.vcf', sep = '\t', quote = F, row.names = F)


# 5. convert vcf to plink bfile format
races <- c('AA', 'CA', 'EA', 'HA')
for (race in races){
  cmd <- paste0('plink2 --vcf 12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/genotype_for_validation_', race, '.vcf --id-delim --rm-dup force-first --make-bed --out 12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/genotype_for_validation_', race)
  system (cmd)
}



# 6. generates a score file for individual level prediction and predict into individual level data
dir.create('12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/')
dir.create('12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/AA')
dir.create('12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/CA')
dir.create('12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/EA')
dir.create('12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/HA')
library(data.table)
library(doParallel)
num_cores <- 56
# Initialize the parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Read AA_models data
AA_models <- fread('9_generate_model_files/AA_methylation_models_hg38.pos', data.table = FALSE)

# Define the function to process each iteration
process_iteration <- function(row) {
  cmd <- paste0('/mnt/lvm_vol_1/hzhong/R-4.1.3/bin/Rscript 12_validate_models_using_MEC_former_smoking/make_score.R 9_generate_model_files/AA/', row[2], ' > 12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/AA/', row[3])
  system(cmd)
  
  cmd <- paste0('plink --bfile 12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/genotype_for_validation_AA --score ', ' 12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/AA/', row[3], ' 1 2 4 --out 12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/AA/', row[3])
  system(cmd)
}

# Use foreach to iterate in parallel
foreach(w = 1:nrow(AA_models)) %dopar% {
  process_iteration(AA_models[w, ])
}

# Read CA_models data
CA_models <- fread('9_generate_model_files/CA_methylation_models_hg38.pos', data.table = FALSE)

# Define the function to process each iteration
process_iteration <- function(row) {
  cmd <- paste0('/mnt/lvm_vol_1/hzhong/R-4.1.3/bin/Rscript 12_validate_models_using_MEC_former_smoking/make_score.R 9_generate_model_files/CA/', row[2], ' > 12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/CA/', row[3])
  system(cmd)
  
  cmd <- paste0('plink --bfile 12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/genotype_for_validation_CA --score ', ' 12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/CA/', row[3], ' 1 2 4 --out 12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/CA/', row[3])
  system(cmd)
}

# Use foreach to iterate in parallel
foreach(w = 1:nrow(CA_models)) %dopar% {
  process_iteration(CA_models[w, ])
}

# Read EA_models data
EA_models <- fread('9_generate_model_files/EA_methylation_models_hg38.pos', data.table = FALSE)

# Define the function to process each iteration
process_iteration <- function(row) {
  cmd <- paste0('/mnt/lvm_vol_1/hzhong/R-4.1.3/bin/Rscript 12_validate_models_using_MEC_former_smoking/make_score.R 9_generate_model_files/EA/', row[2], ' > 12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/EA/', row[3])
  system(cmd)
  
  cmd <- paste0('plink --bfile 12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/genotype_for_validation_EA --score ', ' 12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/EA/', row[3], ' 1 2 4 --out 12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/EA/', row[3])
  system(cmd)
}

# Use foreach to iterate in parallel
foreach(w = 1:nrow(EA_models)) %dopar% {
  process_iteration(EA_models[w, ])
}


# Read HA_models data
HA_models <- fread('9_generate_model_files/HA_methylation_models_hg38.pos', data.table = FALSE)

# Define the function to process each iteration
process_iteration <- function(row) {
  cmd <- paste0('/mnt/lvm_vol_1/hzhong/R-4.1.3/bin/Rscript 12_validate_models_using_MEC_former_smoking/make_score.R 9_generate_model_files/HA/', row[2], ' > 12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/HA/', row[3])
  system(cmd)
  
  cmd <- paste0('plink --bfile 12_validate_models_using_MEC_former_smoking/02_SNP_for_validation/genotype_for_validation_HA --score ', ' 12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/HA/', row[3], ' 1 2 4 --out 12_validate_models_using_MEC_former_smoking/03_predicted_phenotype/HA/', row[3])
  system(cmd)
}

# Use foreach to iterate in parallel
foreach(w = 1:nrow(HA_models)) %dopar% {
  process_iteration(HA_models[w, ])
}







# Stop the parallel backend
stopCluster(cl)

