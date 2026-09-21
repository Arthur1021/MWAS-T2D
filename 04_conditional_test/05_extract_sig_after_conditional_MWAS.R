library(data.table)
library(dplyr)
## ===============================
## 0. Define ancestry mapping
## ===============================

ancestry_map <- data.table(
  model_anc = c("AA", "HA", "CA", "EA"),
  gwas_anc  = c("AFR", "AMR", "EAS", "EUR")
)

anno <- fread('/mnt/Data1/sliu7/data/sliu/project/MESA_methylation_WGS_data/18_Tables/table_anno.hg38_multianno.txt')
anno <- select(anno, Chr, Start, Gene.refGene)
## ===============================
## 1. Define paths
## ===============================

base_project <- "/mnt/Data1/sliu7/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data"

conditional_dir <- "04_MWAS_using_conditional_GWAS"
out_dir <- "05_extract_sig_after_conditional_MWAS"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ===============================
## 2. Loop over ancestries
## ===============================

for (i in seq_len(nrow(ancestry_map))) {

  model_anc <- ancestry_map$model_anc[i]   # AA, HA, CA, EA
  gwas_anc  <- ancestry_map$gwas_anc[i]    # AFR, AMR, EAS, EUR

  cat("\n==============================\n")
  cat("Processing:", model_anc, "using", gwas_anc, "conditional GWAS\n")
  cat("==============================\n")

  ## ===============================
  ## Load original MWAS results
  ## Used for Bonferroni denominator
  ## ===============================

  mwas_file <- file.path(
    base_project,
    "10_mwas_on_T2D_rsq_0.01",
    paste0("summary_mwas_on_T2D_", model_anc, ".txt")
  )

  if (!file.exists(mwas_file)) {
    cat("Original MWAS file not found:", mwas_file, "\n")
    next
  }

  df <- fread(mwas_file)

  bonf_threshold <- 0.05 / nrow(df)
  df_sig <- df[df$TWAS.P < bonf_threshold, ]
  cat("Bonferroni threshold:", bonf_threshold, "\n")

  ## ===============================
  ## Read conditional MWAS results
  ## ===============================

  files <- Sys.glob(
    file.path(
      conditional_dir,
      paste0("mwas_", model_anc, "_T2D_cond_", gwas_anc, "_chr*.txt")
    )
  )

  if (length(files) == 0) {
    cat("No conditional MWAS files found for", model_anc, gwas_anc, "\n")
    next
  }

  cat("Number of conditional MWAS files:", length(files), "\n")

  all <- rbindlist(
    lapply(files, function(file) {
      cat("Reading:", file, "\n")
      fread(file)
    }),
    fill = TRUE
  )

  cat("Total conditional MWAS rows:", nrow(all), "\n")

  ## ===============================
  ## Filter significant results
  ## ===============================

  if (!"TWAS.P" %in% colnames(all)) {
    cat("Column TWAS.P not found for", model_anc, ". Skipping.\n")
    next
  }

  all_filter <- all[TWAS.P < bonf_threshold]
  all_filter <- all_filter[all_filter$ID %in% df_sig$ID, ]
  all_filter <- select(all_filter, ID, CHR, P0, NSNP, NWGT, MODELCV.R2, TWAS.Z, TWAS.P)
  all_filter <- left_join(all_filter, anno, by = c('CHR' = 'Chr', 'P0' = 'Start'))
  all_filter <- left_join(all_filter, df_sig[,c(2, 13,14)], by = c('ID' = 'ID'))
  all_filter <- select(all_filter, ID, CHR, P0, NSNP, NWGT, MODELCV.R2, TWAS.Z.y, TWAS.P.y, TWAS.Z.x, TWAS.P.x, Gene.refGene)
  cat("Significant after conditional MWAS:", nrow(all_filter), "\n")

  ## ===============================
  ## Save output
  ## ===============================

  out_file <- file.path(
    out_dir,
    paste0("conditional_significant_T2D_", model_anc, "_", gwas_anc, ".txt")
  )

  fwrite(all_filter, out_file, sep = ",")

  cat("Saved:", out_file, "\n")
}

cat("\nAll ancestries finished!\n")