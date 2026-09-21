library(data.table)

## ===============================
## 0. Define ancestry mapping
## ===============================

ancestry_map <- data.table(
  model_anc = c("AA", "HA", "CA", "EA"),
  gwas_anc  = c("AFR", "AMR", "EAS", "EUR")
)

## ===============================
## 1. Common paths
## ===============================

base_project <- "/mnt/Data1/sliu7/mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data"

fusion_script <- "/mnt/Data1/sliu7/mnt/lvm_vol_2/sliu/pipeline/TWAS_fusion/bin/FUSION.assoc_test.R"

out_dir <- "04_MWAS_using_conditional_GWAS/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ===============================
## 2. Loop over ancestries
## ===============================

for (i in seq_len(nrow(ancestry_map))) {

  model_anc <- ancestry_map$model_anc[i]   # AA, HA, CA, EA
  gwas_anc  <- ancestry_map$gwas_anc[i]    # AFR, AMR, EAS, EUR

  cat("\n===============================\n")
  cat("Processing model ancestry:", model_anc, "\n")
  cat("Using conditional GWAS ancestry:", gwas_anc, "\n")
  cat("===============================\n")

  ## ===============================
  ## 2.1 Load original MWAS results
  ## ===============================

  mwas_file <- file.path(
    base_project,
    "10_mwas_on_T2D_rsq_0.01",
    paste0("summary_mwas_on_T2D_", model_anc, ".txt")
  )

  if (!file.exists(mwas_file)) {
    cat("MWAS file not found:", mwas_file, "\n")
    next
  }

  df <- fread(mwas_file)

  ## Select significant CpGs using Bonferroni correction
  df_sig <- df[df$FDR < 0.05,]

  cat("Significant CpGs for", model_anc, ":", nrow(df_sig), "\n")

  if (nrow(df_sig) == 0) {
    cat("No significant CpGs for", model_anc, ". Skipping.\n")
    next
  }

  ## ===============================
  ## 2.2 Extract CpG IDs and chromosomes
  ## ===============================

  sig_ids <- df_sig$ID
  sig_chr <- sort(unique(df_sig$CHR_38))

  cat("Chromosomes used for", model_anc, ":", paste(sig_chr, collapse = ", "), "\n")

  ## ===============================
  ## 2.3 Filter .pos file
  ## ===============================

  pos_file <- file.path(
    base_project,
    "9_generate_model_files",
    paste0(model_anc, "_methylation_models_hg38.pos")
  )

  if (!file.exists(pos_file)) {
    cat("POS file not found:", pos_file, "\n")
    next
  }

  pos <- fread(pos_file)

  pos_sig <- pos[pos$ID %in% sig_ids]

  cat("CpGs retained in pos file for", model_anc, ":", nrow(pos_sig), "\n")

  if (nrow(pos_sig) == 0) {
    cat("No CpGs matched between significant MWAS results and pos file for", model_anc, ". Skipping.\n")
    next
  }

  out_pos <- file.path(
    out_dir,
    paste0(model_anc, "_methylation_models_sig.pos")
  )

  write.table(
    pos_sig,
    out_pos,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("Filtered pos file saved:", out_pos, "\n")

  ## ===============================
  ## 2.4 Prepare FUSION input paths
  ## ===============================

  sumstats <- file.path(
    "03_prepare_conditional_GWAS",
    paste0(gwas_anc, "_TWAS_input.txt")
  )

  weights_dir <- file.path(
    base_project,
    "9_generate_model_files",
    model_anc
  )

  ld_prefix <- paste0(
    "/mnt/Data1/sliu7/data/sliu/database/LD-reference/",
    gwas_anc,
    "_chr_hg38/1000G.",
    gwas_anc,
    ".ALLSNP.QC."
  )

  if (!file.exists(sumstats)) {
    cat("Conditional GWAS sumstats not found:", sumstats, "\n")
    next
  }

  if (!dir.exists(weights_dir)) {
    cat("Weights directory not found:", weights_dir, "\n")
    next
  }

  ## ===============================
  ## 2.5 Run FUSION for selected chromosomes
  ## ===============================

  for (chr in sig_chr) {

    out_file <- file.path(
      out_dir,
      paste0("mwas_", model_anc, "_T2D_cond_", gwas_anc, "_chr", chr, ".txt")
    )

    log_file <- file.path(
      out_dir,
      paste0("mwas_", model_anc, "_T2D_cond_", gwas_anc, "_chr", chr, ".out")
    )

    ## Skip if result already exists and is not empty
    if (file.exists(out_file) && file.info(out_file)$size > 0) {
      cat("Result already exists. Skipping:", out_file, "\n")
      next
    }

    cat("Submitting", model_anc, "using", gwas_anc, "conditional GWAS, chr", chr, "\n")

    cmd <- paste(
      "nohup Rscript", fusion_script,
      "--sumstats", sumstats,
      "--weights", out_pos,
      "--weights_dir", weights_dir,
      "--ref_ld_chr", ld_prefix,
      "--chr", chr,
      "--out", out_file,
      ">",
      log_file,
      "2>&1 &"
    )

    system(cmd)
  }

  cat("Finished submitting jobs for", model_anc, "with", gwas_anc, "conditional GWAS.\n")
}

cat("\nAll ancestry jobs checked/submitted!\n")