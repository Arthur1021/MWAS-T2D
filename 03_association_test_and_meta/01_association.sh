
################################
### 6. mwas on T2D Rsq 0.01 ###
################################
mkdir 10_mwas_on_T2D_rsq_0.01
#AA
for i in $(seq 1 22); do {
    nohup /mnt/lvm_vol_1/hzhong/R-4.1.3/bin/Rscript /mnt/lvm_vol_2/sliu/pipeline/TWAS_fusion/bin/FUSION.assoc_test.R --sumstats /mnt/lvm_vol_2/sliu/database/GWAS_summary/T2D_2023/T2D_2023_GWAS_summary_AFR.txt --weights 9_generate_model_files/AA_methylation_models_hg38.pos --weights_dir 9_generate_model_files/AA/ --ref_ld_chr /mnt/lvm_vol_2/sliu/database/LD-reference/AFR_chr_hg38/1000G.AFR.ALLSNP.QC. --chr ${i} --out 10_mwas_on_T2D_rsq_0.01/mwas_AA_T2D_${i}.txt > 10_mwas_on_T2D_rsq_0.01/mwas_AA_T2D_${i}.out&
} done

#CA
for i in $(seq 1 22); do {
    nohup /mnt/lvm_vol_1/hzhong/R-4.1.3/bin/Rscript /mnt/lvm_vol_2/sliu/pipeline/TWAS_fusion/bin/FUSION.assoc_test.R --sumstats /mnt/lvm_vol_2/sliu/database/GWAS_summary/T2D_2023/T2D_2023_GWAS_summary_EAS.txt --weights 9_generate_model_files/CA_methylation_models_hg38.pos --weights_dir 9_generate_model_files/CA/ --ref_ld_chr /mnt/lvm_vol_2/sliu/database/LD-reference/EAS_chr_hg38/1000G.EAS.ALLSNP.QC. --chr ${i} --out 10_mwas_on_T2D_rsq_0.01/mwas_CA_T2D_${i}.txt > 10_mwas_on_T2D_rsq_0.01/mwas_CA_T2D_${i}.out&
} done

#EA
for i in $(seq 1 22); do {
    nohup /mnt/lvm_vol_1/hzhong/R-4.1.3/bin/Rscript /mnt/lvm_vol_2/sliu/pipeline/TWAS_fusion/bin/FUSION.assoc_test.R --sumstats /mnt/lvm_vol_2/sliu/database/GWAS_summary/T2D_2023/T2D_2023_GWAS_summary_EUR.txt --weights 9_generate_model_files/EA_methylation_models_hg38.pos --weights_dir 9_generate_model_files/EA/ --ref_ld_chr /mnt/lvm_vol_2/sliu/database/LD-reference/EUR_chr_hg38/1000G.EUR.ALLSNP.QC. --chr ${i} --out 10_mwas_on_T2D_rsq_0.01/mwas_EA_T2D_${i}.txt > 10_mwas_on_T2D_rsq_0.01/mwas_EA_T2D_${i}.out&
} done

#HA 
for i in $(seq 1 22); do {
    nohup /mnt/lvm_vol_1/hzhong/R-4.1.3/bin/Rscript /mnt/lvm_vol_2/sliu/pipeline/TWAS_fusion/bin/FUSION.assoc_test.R --sumstats /mnt/lvm_vol_2/sliu/database/GWAS_summary/T2D_2023/T2D_2023_GWAS_summary_AMR.txt --weights 9_generate_model_files/HA_methylation_models_hg38.pos --weights_dir 9_generate_model_files/HA/ --ref_ld_chr /mnt/lvm_vol_2/sliu/database/LD-reference/AMR_chr_hg38/1000G.AMR.ALLSNP.QC. --chr ${i} --out 10_mwas_on_T2D_rsq_0.01/mwas_HA_T2D_${i}.txt > 10_mwas_on_T2D_rsq_0.01/mwas_HA_T2D_${i}.out&
} done

