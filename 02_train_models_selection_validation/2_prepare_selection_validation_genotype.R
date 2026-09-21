library(data.table)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)
library(lumi)

#####################################################################
###### 1. get samples with methylation data (exam1), covariate ######
#####################################################################
if(!dir.exists('1_covariate')){
    dir.create('1_covariate')
}

#sample with methylation
load('/home/sliu/work/project/MESA_methylation_WGS_data/phs001416_TOPMed_Methylomics_MESA/TOPMed_MESA_normalized_methylation_files/norm.beta.pc11.Robj')
methylation_norm_sample <- colnames(norm.beta)
#remove dup sample
duplicate_sample <- fread('/home/sliu/work/project/MESA_methylation_WGS_data/phs001416_TOPMed_Methylomics_MESA/TOPMed_MESA_normalized_methylation_files/remove_duplicates.txt', data.table = F , header =FALSE)
methylation_norm_sample_no_dup <- methylation_norm_sample[!(methylation_norm_sample %in% duplicate_sample$V1)]

#1 CAU
#2 CHI
#3 AFR
#4 HISP


#only select exam1 sample #935 samples
methylation_norm_sample_no_dup_visit1 <- methylation_norm_sample_no_dup[stringr::str_detect(methylation_norm_sample_no_dup, 'exam1')]
map <- fread('/home/sliu/work/project/MESA_methylation_WGS_data/phs001416_TOPMed_Methylomics_MESA/TOPMed_MESA_normalized_methylation_files/TOPMed_MESA.methylomics.samplesheet_with_feno.mixup_fix.only_BIS.pass_QC.txt', data.table = F)
map_filter <- map[map$sidno_exam.ID %in% methylation_norm_sample_no_dup_visit1,]

map_filter <- select(map_filter, 'sidno', 'Sample_Name', 'NWDID', 'race', 'sex' , 'age')

#select sample with cell counts
cell_type <- fread('/home/sliu/work/project/MESA_methylation_WGS_data/phs001416_TOPMed_Methylomics_MESA/TOPMed_MESA_normalized_methylation_files/estimated_cellcounts_houseman.txt', data.table = F)
colnames(cell_type)[1] <- 'Sample_Name'
cell_type <- select(cell_type, 'Sample_Name', 'Bcell', 'CD4T', 'CD8T', 'Eos', 'Mono', 'Neu', 'NK')

map_filter_cell_type <- merge(map_filter, cell_type, by = 'Sample_Name')

#select sample with covariate
covariate <- fread('/mnt/lvm_vol_1/langwu/MESA_methylome_proteome_WGS_data/covariates/MESAe1and5_CoVarLangWu_20210630/MESAe1and5_CoVarLangWu_20210630.txt', data.table = F)
covariate <- select(covariate, 'sidno', 'age1c', 'gender1', 'htcm1', 'wtlb1', 'cig1c', 'pkyrs1c')
map_filter_cell_type_covariate <- merge(map_filter_cell_type, covariate, by = 'sidno')


#check age1c = age, gender1 = sex, calculate BMI
sum(map_filter_cell_type_covariate$age1c == map_filter_cell_type_covariate$age) == nrow(map_filter_cell_type_covariate)
sum(map_filter_cell_type_covariate$gender1 == map_filter_cell_type_covariate$sex) == nrow(map_filter_cell_type_covariate)
map_filter_cell_type_covariate$BMI <- map_filter_cell_type_covariate$wtlb1/map_filter_cell_type_covariate$htcm1/map_filter_cell_type_covariate$htcm1*703

res <- map_filter_cell_type_covariate[c(3:13,18,19,20)]

# not run
# polymorphic_df <- fread('/home/sliu/work/project/MESA_methylation_WGS_data/phs001416_TOPMed_Methylomics_MESA/TOPMed_MESA_normalized_methylation_files/EPIC.hg38.manifest.polymorphic_in_MESA.txt.gz', data.table =F)


#write out result
write.table(res, '1_covariate/covariate.txt', sep ='\t', quote = F, row.names =F)

write.table(res$NWDID, '1_covariate/sample_with_covariate.txt', sep ='\t', quote = F, row.names =F, col.names = F)

# #AA
# write.table(map_filter_cell_type_covariate[map_filter_cell_type_covariate$race == 3,]$NWDID, '1_covariate/AA_sample.txt', sep ='\t', quote = F, row.names =F, col.names = F)
# #CA
# write.table(map_filter_cell_type_covariate[map_filter_cell_type_covariate$race == 2,]$NWDID, '1_covariate/CA_sample.txt', sep ='\t', quote = F, row.names =F, col.names = F)
# #EA
# write.table(map_filter_cell_type_covariate[map_filter_cell_type_covariate$race == 1,]$NWDID, '1_covariate/EA_sample.txt', sep ='\t', quote = F, row.names =F, col.names = F)
# #HA
# write.table(map_filter_cell_type_covariate[map_filter_cell_type_covariate$race == 4,]$NWDID, '1_covariate/HA_sample.txt', sep ='\t', quote = F, row.names =F, col.names = F)


#check samples with genotype, all samples have genotype 
sample_genotype <- fread('/mnt/lvm_vol_1/langwu/MESA_methylome_proteome_WGS_data/wgs/minDP10/freeze.7a.chr1.pass_and_fail.gtonly.minDP10.fam', data.table = F)
nrow(res) == sum(res$NWDID %in% sample_genotype$V2)



###########################################################
###### 2. identify unrelated samples for each ethnic ######
###########################################################

if(!dir.exists('2_genotype')){
    dir.create('2_genotype')
}

#filter genotype with for 935 samples
geno_files <- Sys.glob('/mnt/lvm_vol_1/langwu/MESA_methylome_proteome_WGS_data/wgs/minDP10/*.vcf')
geno_files <- geno_files[1:22]



genotype_select_sample <- function (infile){
    tmp <- last(unlist(str_split(infile, '/'))) %>% str_split('[.]') %>% unlist()
    outfile <- tmp[3]
    cmd <- paste0('plink2 --vcf ', infile, ' --var-filter --snps-only --keep 1_covariate/sample_with_covariate.txt --make-bed --out 2_genotype/', outfile)
    system(cmd)
}

registerDoParallel(50)
foreach (x = iter(geno_files)) %dopar% {
  genotype_select_sample(x)
}


#filter genotype with maf 0.01 for each ethnic AA CA EA HA


# race <- unlist(list('AA', 'CA', 'EA', 'HA'))
# genotype_maf <- function(infile, race){
#     sample_file <- paste0('1_covariate/', race , '_sample.txt')
#     cmd <- paste0('plink2 --bfile 2_genotype/', infile, ' --maf 0.01 --keep ', sample_file, ' --make-bed --out 2_genotype/', race, '_', infile, '_maf_0.01')
#     print (cmd)
#     system (cmd)
# }
# for (r in race){
#     foreach (i = 1:22) %dopar% {
#     genotype_maf(paste0('chr', i), r)
#     }
# }
# do not use above code, need to seperate into 2 steps, first sample then maf filter

# Define the ethnic groups
race <- c('AA', 'CA', 'EA', 'HA')

# Define the function to process in two steps
genotype_keep_maf <- function(infile, race) {
    sample_file <- paste0('1_covariate/', race, '_sample.txt')
    
    # Step 1: Filter by sample file (--keep)
    step1_outfile <- paste0('2_genotype/', race, '_', infile, '_keep_tmp')
    cmd1 <- paste0('plink2 --bfile 2_genotype/', infile, 
                   ' --keep ', sample_file, 
                   ' --make-bed --out ', step1_outfile)
    print(cmd1)
    system(cmd1)
    
    # Step 2: Filter by MAF (--maf 0.01)
    step2_outfile <- paste0('2_genotype/', race, '_', infile, '_maf_0.01')
    cmd2 <- paste0('plink2 --bfile ', step1_outfile, 
                   ' --maf 0.01 --make-bed --out ', step2_outfile)
    print(cmd2)
    system(cmd2)
}

# Parallel processing for all chromosomes and ethnic groups
registerDoParallel(cores = 50)
foreach(r = race, .combine = c) %dopar% {
    foreach(i = 1:22) %do% {
        genotype_keep_maf(paste0('chr', i), r)
    }
}

#filter genotype with geno 0.05, add snp id with chr:pos, remove duplicate id for each ethnic AA CA EA HA
genotype_geno <- function(infile, race){
    cmd <- paste0('plink2 --bfile 2_genotype/',race , '_' ,infile, '_maf_0.01 --geno 0.05 --set-missing-var-ids @:# --rm-dup exclude-all --make-bed --out 2_genotype/', race ,'_',infile, '_maf_0.01_geno_0.05')
    system (cmd)
}
for (r in race){
    foreach (i = 1:22) %dopar% {
    genotype_geno(paste0('chr', i),r)
    }
}

#merge chromosomes for each ethnic AA CA EA HA

for (r in race){
    merge_list <- list()
    for (i in 2:22){
        merge_list <- append(merge_list, paste0('2_genotype/',r,'_chr', i, '_maf_0.01_geno_0.05'))
    }
    write.table(unlist(merge_list), paste0('2_genotype/',r, '_merge_list.txt'), sep ='\t', quote = F, row.names =F, col.names = F)
    rm (merge_list)
}

#merge for each ethnic AA CA EA HA
for (r in race){
    cmd <- paste0('plink --bfile 2_genotype/', r, '_chr1_maf_0.01_geno_0.05 --merge-list 2_genotype/', r, '_merge_list.txt --make-bed --out 2_genotype/', r, '_All_chr_maf_0.01_geno_0.05')
    system (cmd)

}

#identify unrelated samples
for (r in race){
    cmd <- paste0('python /mnt/lvm_vol_1/sliu/work/pipeline/TWAS_fusion/relatedness_analysis_plink.py -in_gen 2_genotype/', r, '_All_chr_maf_0.01_geno_0.05 -odir 2_genotype/independent_', r)
    system (cmd)
}




#######
# get snps used in models
# nohup python3 get_snps_in_model.py > get_snps_in_model.out


########################################################################
###### 3. generate PCA for selection and validation data ######
########################################################################

race_analysis <- c('3_AA', '4_CA', '5_EA', '6_HA')
for (r in race_analysis){
    if (!dir.exists(r)){
        dir.create(r)
    }
}

# filter maf with 0.01, add snp id with chr:pos:allele:allele (allele ASCII-sort order)
command <- vector()
for (r in race_analysis){
    race <- r %>% str_split('_') %>% unlist() %>% last()
    id_add <- paste('@:#:', '$1:', '$2', sep = '\\')
    for (i in 1:22){
        cmd <- paste0('plink --bfile 2_genotype/chr',i , ' --keep 2_genotype/independent_', race, '/chr1-22_independent.fam --maf 0.05 --set-missing-var-ids @:#:\\$1:\\$2 --make-bed  --out ', r, '/chr',i,'_maf_0.01')
        command <- append(command, cmd)
    }
}

foreach (x = iter(command)) %dopar% {
    system(x)
}

# filter geno with 0.01
command <- vector()
for (r in race_analysis){
    for (i in 1:22){
        cmd <- paste0('plink2 --bfile ', r, '/chr', i, '_maf_0.01 --geno 0.05  --make-bed --out ', r, '/chr', i, '_maf_0.01_geno_0.05')
        command <- append(command, cmd)
    }
}

foreach (x = iter(command)) %dopar% {
    system(x)
}

#merge chromosomes for each ethnic AA CA EA HA

for (r in race_analysis){
    merge_list <- list()
    for (i in 2:22){
        merge_list <- append(merge_list, paste0(r,'/chr', i, '_maf_0.01_geno_0.05'))
    }
    write.table(unlist(merge_list), paste0(r, '/merge_list.txt'), sep ='\t', quote = F, row.names =F, col.names = F)

    cmd <- paste0('plink --bfile ', r, '/chr1_maf_0.01_geno_0.05 --merge-list ', r, '/merge_list.txt --make-bed --out ', r, '/All_chr_maf_0.01_geno_0.05')
    system(cmd)
    rm (merge_list)
}




#randomly extract half samples as selection cohort and the remaining as validation cohort
rm(r)
#for ( r in race_analysis){
split_and_PCA <- function(r){    
    vector <- fread(paste0(r, '/All_chr_maf_0.01_geno_0.05.fam'), data.table = F)
    vector <- vector$V2
    vector_length <- length(vector)
    set.seed(12345)
    random_index <- sample(1:vector_length)
    split_point <- vector_length %/% 2
    selection_cohort <- vector[random_index[1:split_point]]
    selection_cohort <- data.frame(selection_cohort)
    colnames(selection_cohort) <- 'IID'
    selection_cohort$FID <- 0
    selection_cohort <- select(selection_cohort, 'FID', 'IID')
    write.table(selection_cohort, paste0(r,'/selection_sample.txt'), sep ='\t', quote = F, row.names =F, col.names = TRUE)
    
    #extract selection sample geno
    cmd <- paste0('plink2 --bfile ', r, '/All_chr_maf_0.01_geno_0.05 --keep ',r ,'/selection_sample.txt --make-bed --out ', r, '/All_selection_filter')
    system(cmd)


    

    if (r == '4_EA'){
        race_para <- 'Y'
    }else{
        race_para <- 'N'
    }

    #PCA matrix for selection samples
    cmd <- paste0('python3 /mnt/lvm_vol_1/sliu/work/pipeline/TWAS_fusion/PCA_using_EIGENSOFT_update.py -in_gen ',r ,'/All_selection_filter -Eur ',race_para, ' -odir ', r, '/selection_PCA'  )
    system(cmd)


    validation_cohort <- vector[random_index[(split_point + 1):vector_length]]
    validation_cohort <- data.frame(validation_cohort)
    colnames(validation_cohort) <- 'IID'
    validation_cohort$FID <- 0
    validation_cohort <- select(validation_cohort, 'FID', 'IID')
    write.table(validation_cohort, paste0(r,'/validation_sample.txt'), sep ='\t', quote = F, row.names =F, col.names = TRUE)
    cmd <- paste0('plink2 --bfile ', r, '/All_chr_maf_0.01_geno_0.05 --keep ',r ,'/validation_sample.txt --make-bed --out ', r, '/All_validation_filter')
    system(cmd)
    #PCA matrix for validation samples
    cmd <- paste0('python /mnt/lvm_vol_1/sliu/work/pipeline/TWAS_fusion/PCA_using_EIGENSOFT_update.py -in_gen ',r ,'/All_validation_filter -Eur ',race_para, ' -odir ', r, '/validation_PCA'  )
    system(cmd)

}

registerDoParallel(50)
foreach (r = race_analysis) %dopar% {
    split_and_PCA(r)
}

##############################################################################################
#### 4. split selection and validation sample geno into chromosome from raw genotype data ####
##############################################################################################

race_analysis <- c('3_AA', '4_CA', '5_EA', '6_HA')
# Define function to process a single chromosome for a specific race
process_chromosome <- function(chr_n, r) {
    r_abbre <- str_split(r, '_') %>% unlist()
    if (r_abbre[2] == 'CA'){
       r_abbre[2] <- 'JA'
    }
    
    # Selection
    cmd <- paste0('plink --bfile 2_genotype/chr', chr_n, 
                  ' --keep ', r, '/selection_sample.txt',
                  ' --set-missing-var-ids @:#:\\$1:\\$2',
                  ' --update-name 1_QC/SNP_chr_pos_', r_abbre[2], '.txt',
                  ' --extract 1_QC/SNP_used_in_', r_abbre[2], '.txt',
                  ' --make-bed --out ', r, '/All_selection_chr', chr_n)
    system(cmd)
    
    # Validation
    cmd <- paste0('plink --bfile 2_genotype/chr', chr_n, 
                  ' --keep ', r, '/validation_sample.txt',
                  ' --set-missing-var-ids @:#:\\$1:\\$2',
                  ' --update-name 1_QC/SNP_chr_pos_', r_abbre[2], '.txt',
                  ' --extract 1_QC/SNP_used_in_', r_abbre[2], '.txt',
                  ' --make-bed --out ', r, '/All_validation_chr', chr_n)
    system(cmd)
}

# Parallel processing for each race
for (r in race_analysis) {
    mclapply(1:22, function(chr_n) process_chromosome(chr_n, r), mc.cores = 50)
}


###################################################################
###### 4. extract methylation and covariate for each ethnic  ######
###################################################################


#methylation QC
load('/home/sliu/work/project/MESA_methylation_WGS_data/phs001416_TOPMed_Methylomics_MESA/TOPMed_MESA_normalized_methylation_files/norm.beta.pc11.Robj')
load('/home/sliu/work/project/MESA_methylation_WGS_data/phs001416_TOPMed_Methylomics_MESA/TOPMed_MESA_normalized_methylation_files/detection_p_value_matrix.RData')
norm.beta[det > 0.01] <- NA
norm.beta <- as.data.frame(norm.beta)
EPIC_hg38_manifest <- fread('EPIC.hg38.manifest.tsv', data.table =F)
probes_exclude <- filter(EPIC_hg38_manifest, EPIC_hg38_manifest$MASK_mapping == 'TRUE' | EPIC_hg38_manifest$MASK_sub30_copy == 'TRUE' | EPIC_hg38_manifest$MASK_typeINextBaseSwitch == 'TRUE' | EPIC_hg38_manifest$MASK_extBase == 'TRUE' | EPIC_hg38_manifest$probeType != 'cg')

race <- c('3_AA', '4_CA', '5_EA', '6_HA')

M_value_regression_out_coverate <- function(r){
    print (r)
    sample_in_race <- fread(paste0(r, '/All_chr_maf_0.01_geno_0.05.fam'), header = F)
    match_tmp <- filter(map_filter_cell_type_covariate, map_filter_cell_type_covariate$NWDID %in% sample_in_race$V2)
    methylation_beta_select <- norm.beta %>% select(match_tmp$Sample_Name)
    
    # calculate probe CR
    callrate <- apply(methylation_beta_select, 1, function(x){sum(!is.na(x))})
    callrate <- callrate/ncol(methylation_beta_select)

    # keep probes with CR > 90%
    methylation_beta_select <- methylation_beta_select[callrate > 0.9, ]


    methylation_beta_select_good <- methylation_beta_select[!(rownames(methylation_beta_select) %in% probes_exclude$probeID),]
    m_value <- beta2m(methylation_beta_select_good)

    headerline <- as.data.frame(match_tmp$NWDID)
    rownames(headerline) <- match_tmp$Sample_Name
    colnames(headerline) <- 'NWDID'

    if (identical(colnames(m_value), rownames(headerline))){
        colnames(m_value) <- c(headerline$NWDID)
    }else{
        break
    }

    #write result 
    m_value <- cbind(rownames(m_value), m_value)
    colnames(m_value)[1] <- 'ID'

    write.table (m_value, paste0(r, '/methylation_m_value.txt'),  sep ='\t', quote = F, row.names = F)

}
registerDoParallel(50)
foreach (r = race) %dopar% {
    M_value_regression_out_coverate(r)
}

