
# AA 156738
# CA 427504
# EA 1812017
# HA 88743
# prepare input files
library(data.table)
library(dplyr)
AA <- fread('../10_mwas_on_T2D_rsq_0.01/summary_mwas_on_T2D_AA.txt', data.table = F)
AA <- select(AA, ID, TWAS.Z, TWAS.P)

AA$WEIGHT <- 156738
CA <- fread('../10_mwas_on_T2D_rsq_0.01/summary_mwas_on_T2D_CA.txt', data.table = F)
CA <- select(CA, ID, TWAS.Z, TWAS.P)
CA$WEIGHT <- 427504
EA <- fread('../10_mwas_on_T2D_rsq_0.01/summary_mwas_on_T2D_EA.txt', data.table = F)
EA <- select(EA, ID, TWAS.Z, TWAS.P)
EA$WEIGHT <- 1812017
HA <- fread('../10_mwas_on_T2D_rsq_0.01/summary_mwas_on_T2D_HA.txt', data.table = F)
HA <- select(HA, ID, TWAS.Z, TWAS.P)
HA$WEIGHT <- 88743

fwrite(AA, 'MWAS_on_T2D_AA.txt', sep = '\t')
fwrite(CA, 'MWAS_on_T2D_CA.txt', sep = '\t')
fwrite(EA, 'MWAS_on_T2D_EA.txt', sep = '\t')
fwrite(HA, 'MWAS_on_T2D_HA.txt', sep = '\t')

# perform meta analysis
system('/mnt/lvm_vol_2/sliu/software/generic-metal/metal metal_T2D_MWAS.txt')


# merge all results
meta <- fread('meta_T2D_MWAS_result1.tbl', data.table = F)
meta <- select(meta, MarkerName, Zscore, 'P-value', Direction, HetPVal)
colnames(meta)[1] <- 'ID'
meta <- meta[order(meta$'P-value'),]
meta$meta_FDR <- p.adjust(meta$'P-value', method = 'fdr')
meta$meta_Bonferroni <- p.adjust(meta$'P-value', method = 'bonferroni')

library(dplyr)

# Assuming 'common_column' is the column you want to merge on
merged_data <- left_join(meta, AA[1:3], by = "ID")
merged_data <- left_join(merged_data, CA[1:3], by = "ID")
merged_data <- left_join(merged_data, EA[1:3], by = "ID")
merged_data <- left_join(merged_data, HA[1:3], by = "ID")

colnames(merged_data) <- c('CpG', 'meta_Z', 'meta_P', 'meta_Direction', 'meta_HetPVal', 'meta_FDR', 'meta_Bonferroni', 'AA_Z', 'AA_P', 'CA_Z', 'CA_P', 'EA_Z', 'EA_P', 'HA_Z', 'HA_P')


hg38 <- fread('../EPIC.hg38.manifest.tsv', data.table = F)
hg38 <- select(hg38, probeID, CpG_chrm, CpG_beg, CpG_end)
colnames(hg38) <- c('ID', 'CHR_38', 'P0_38', 'P1_38')
hg19 <- fread('../EPIC.hg19.manifest.txt', data.table = F)
hg19 <- select(hg19, V4, V1, V2, V3)
colnames(hg19) <- c('ID', 'CHR_37', 'P0_37', 'P1_37')

merged_data <- left_join(merged_data, hg38, by = c('CpG' = 'ID'))
merged_data <- left_join(merged_data, hg19, by = c('CpG' = 'ID'))
merged_data <- select(merged_data, CpG, CHR_38, P0_38, P1_38, CHR_37, P0_37, P1_37, meta_Z, meta_P, meta_Direction, meta_HetPVal, meta_FDR, meta_Bonferroni, AA_Z, AA_P, CA_Z, CA_P, EA_Z, EA_P, HA_Z, HA_P)


write.table(merged_data, 'all_CpG_meta_results.txt', sep = '\t', quote = F, row.names =F)

merged_data <- filter(merged_data, meta_FDR < 0.05)
write.table(merged_data, 'sig_CpG_meta_results.txt', sep = '\t', quote = F, row.names =F)

merged_data_2 <- filter(merged_data, meta_Bonferroni < 0.05)

AA$Bonferroni_p <- p.adjust(AA$TWAS.P , method = 'bonferroni')
CA$Bonferroni_p <- p.adjust(CA$TWAS.P , method = 'bonferroni')
EA$Bonferroni_p <- p.adjust(EA$TWAS.P , method = 'bonferroni')
HA$Bonferroni_p <- p.adjust(HA$TWAS.P , method = 'bonferroni')

merged_data_2$novel <- ifelse(merged_data_2$CpG %in% c(AA[AA$Bonferroni_p < 0.05,]$ID, CA[CA$Bonferroni_p < 0.05,]$ID, EA[EA$Bonferroni_p < 0.05,]$ID, HA[HA$Bonferroni_p < 0.05,]$ID), 'none', 'novel')
write.table(merged_data_2, 'sig_CpG_meta_results_bonferroni_novel.txt', sep = '\t', quote = F, row.names =F)