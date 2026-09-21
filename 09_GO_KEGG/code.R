library(missMethyl)
library(data.table)
library(dplyr)
# files contain one CpG ID per line, no header required

AA <- fread('/mnt/Data1/sliu7/data/sliu/project/MESA_methylation_WGS_data/10_mwas_on_T2D_rsq_0.01/summary_mwas_on_T2D_AA.txt')
AA_sig <- AA[AA$TWAS.P < 0.05/nrow(AA), ]$ID

CA <- fread('/mnt/Data1/sliu7/data/sliu/project/MESA_methylation_WGS_data/10_mwas_on_T2D_rsq_0.01/summary_mwas_on_T2D_CA.txt')
CA_sig <- CA[CA$TWAS.P < 0.05/nrow(CA), ]$ID

EA <- fread('/mnt/Data1/sliu7/data/sliu/project/MESA_methylation_WGS_data/10_mwas_on_T2D_rsq_0.01/summary_mwas_on_T2D_EA.txt')
EA_sig <- EA[EA$TWAS.P < 0.05/nrow(EA), ]$ID

HA <- fread('/mnt/Data1/sliu7/data/sliu/project/MESA_methylation_WGS_data/10_mwas_on_T2D_rsq_0.01/summary_mwas_on_T2D_HA.txt')
HA_sig <- HA[HA$TWAS.P < 0.05/nrow(HA), ]$ID

META <- fread('/mnt/Data1/sliu7/data/sliu/project/MESA_methylation_WGS_data/15_Meta_analysis/all_CpG_meta_results.txt')
META_sig <- META[META$meta_P < 0.05/nrow(META), ]$CpG

sig.cpg <- unique(c(AA_sig, CA_sig, EA_sig, HA_sig, META_sig))

all.cpg <- unique(c(AA$ID, CA$ID, EA$ID, HA$ID, META$CpG))



# sig.cpg <- scan("query_cpg_list.txt", what = character(), quiet = TRUE)
# all.cpg <- scan("background_cpg_list.txt", what = character(), quiet = TRUE)

sig.cpg <- unique(sig.cpg)
all.cpg <- unique(all.cpg)

# keep query CpGs that are in background
sig.cpg <- intersect(sig.cpg, all.cpg)

# GO enrichment
go.res <- gometh(
  sig.cpg = sig.cpg,
  all.cpg = all.cpg,
  collection = "GO",
  array.type = "EPIC",
  prior.prob = TRUE
)

go.res <- go.res[order(go.res$FDR), ]
write.table(go.res, "missMethyl_GO_results.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# KEGG enrichment
kegg.res <- gometh(
  sig.cpg = sig.cpg,
  all.cpg = all.cpg,
  collection = "KEGG",
  array.type = "EPIC",
  prior.prob = TRUE
)

kegg.res <- kegg.res[order(kegg.res$FDR), ]
write.table(kegg.res, "missMethyl_KEGG_results.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# FDR < 0.05 results
write.table(subset(go.res, FDR < 0.05),
            "missMethyl_GO_FDR0.05.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(subset(kegg.res, FDR < 0.05),
            "missMethyl_KEGG_FDR0.05.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)