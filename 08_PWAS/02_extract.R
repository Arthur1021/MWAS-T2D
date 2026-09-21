library(data.table)
library(dplyr)

# ============================================================
# 0. SETTINGS
# ============================================================

out_dir <- "/mnt/Data1/sliu7/project/MWAS_T2D/PWAS_MESA_1K/PWAS_overlap_with_MWAS_sig_CpGs_table"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ann_file <- "/mnt/Data1/sliu7/data/sliu/project/MESA_SomaScan_7K_PWAS/01_preprocess_data/annotation/SomaScan_annotation.txt"

mwas_base <- "/mnt/Data1/sliu7/project/MWAS_T2D/PWAS_MESA_1K/nearest_gene_lists_from_sig_CpGs"

pwas_files <- data.table(
  pop = c("AA", "CA", "EA", "HA"),
  pwas_file = c(
    "/mnt/Data1/sliu7/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/03_association/association_T2D_MVP_AA/all_association.out",
    "/mnt/Data1/sliu7/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/03_association/association_T2D_MVP_CA/all_association.out",
    "/mnt/Data1/sliu7/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/03_association/association_T2D_MVP_EA/all_association.out",
    "/mnt/Data1/sliu7/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/03_association/association_T2D_MVP_HA/all_association.out"
  )
)

# ============================================================
# 1. HELPER FUNCTIONS
# ============================================================

clean_gene_vector <- function(x) {
  x <- as.character(x)
  x <- unlist(strsplit(x, ";|,|/"))
  x <- trimws(x)
  x <- x[!is.na(x)]
  x <- x[x != ""]
  x <- x[x != "."]
  unique(x)
}

collapse_unique <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x)]
  x <- x[x != ""]
  x <- x[x != "."]
  paste(sort(x), collapse = ";")
}

get_mwas_cpg_file <- function(pop) {
  file.path(
    mwas_base,
    pop,
    paste0(pop, "_Bonferroni_sig_CpGs_with_nearest_gene.txt")
  )
}

process_one_pop <- function(pop, pwas_file, ann) {
  
  cat("\n============================================================\n")
  cat("Processing:", pop, "\n")
  cat("============================================================\n")
  
  mwas_cpg_file <- get_mwas_cpg_file(pop)
  
  if (!file.exists(pwas_file)) {
    warning("Missing PWAS file: ", pwas_file)
    return(list(result = data.table(), summary = data.table(pop = pop, status = "missing_pwas_file")))
  }
  
  if (!file.exists(mwas_cpg_file)) {
    warning("Missing MWAS CpG file: ", mwas_cpg_file)
    return(list(result = data.table(), summary = data.table(pop = pop, status = "missing_mwas_cpg_file")))
  }
  
  # ----------------------------------------------------------
  # Read PWAS and annotation
  # ----------------------------------------------------------
  
  pwas <- fread(pwas_file)
  
  required_pwas_cols <- c("ID", "PWAS.Z", "PWAS.P")
  missing_pwas_cols <- setdiff(required_pwas_cols, colnames(pwas))
  
  if (length(missing_pwas_cols) > 0) {
    stop("Missing PWAS columns in ", pwas_file, ": ", paste(missing_pwas_cols, collapse = ", "))
  }
  
  pwas_anno <- left_join(
    as.data.frame(pwas),
    as.data.frame(ann),
    by = c("ID" = "SomaId")
  )
  
  pwas_anno <- as.data.table(pwas_anno)
  
  keep_cols <- c(
    "ID",
    "TargetFullName",
    "Target",
    "EntrezGeneSymbol",
    "PWAS.Z",
    "PWAS.P"
  )
  
  keep_cols <- intersect(keep_cols, colnames(pwas_anno))
  pwas_anno <- pwas_anno[, ..keep_cols]
  
  pwas_anno[, pop := pop]
  pwas_anno[, pwas_row_id := .I]
  
  # Expand EntrezGeneSymbol in case one SomaScan target maps to multiple genes
  pwas_gene_long <- pwas_anno[
    ,
    .(MatchedGene = clean_gene_vector(EntrezGeneSymbol)),
    by = pwas_row_id
  ]
  
  pwas_long <- merge(
    pwas_anno,
    pwas_gene_long,
    by = "pwas_row_id",
    all.x = TRUE
  )
  
  # ----------------------------------------------------------
  # Read MWAS significant CpGs with nearest gene
  # ----------------------------------------------------------
  
  mwas_cpg <- fread(mwas_cpg_file)
  
  required_mwas_cols <- c("ID", "CHR_38", "P0_38", "P1_38", "TWAS.Z", "TWAS.P", "FDR", "Nearest Gene")
  missing_mwas_cols <- setdiff(required_mwas_cols, colnames(mwas_cpg))
  
  if (length(missing_mwas_cols) > 0) {
    stop("Missing MWAS columns in ", mwas_cpg_file, ": ", paste(missing_mwas_cols, collapse = ", "))
  }
  
  setnames(
    mwas_cpg,
    old = c("ID", "TWAS.Z", "TWAS.P", "FDR"),
    new = c("CpG_ID", "MWAS.Z", "MWAS.P", "MWAS.FDR")
  )
  
  mwas_cpg[, cpg_row_id := .I]
  
  # Expand Nearest Gene in case one CpG maps to multiple genes
  mwas_gene_long <- mwas_cpg[
    ,
    .(MatchedGene = clean_gene_vector(`Nearest Gene`)),
    by = cpg_row_id
  ]
  
  mwas_long <- merge(
    mwas_cpg,
    mwas_gene_long,
    by = "cpg_row_id",
    all.x = TRUE
  )
  
  # ----------------------------------------------------------
  # Match PWAS proteins and MWAS CpGs by gene
  # ----------------------------------------------------------
  
  overlap_long <- merge(
    pwas_long,
    mwas_long,
    by = "MatchedGene",
    allow.cartesian = TRUE
  )
  
  if (nrow(overlap_long) == 0) {
    
    pop_out_file <- file.path(
      out_dir,
      paste0(pop, "_PWAS_table_with_matched_CpGs.txt")
    )
    
    fwrite(data.table(), pop_out_file, sep = "\t")
    
    summary_dt <- data.table(
      pop = pop,
      status = "no_overlap",
      n_pwas_total = nrow(pwas),
      n_mwas_sig_cpg = uniqueN(mwas_cpg$CpG_ID),
      n_output_rows = 0,
      output_file = pop_out_file
    )
    
    return(list(result = data.table(), summary = summary_dt))
  }
  
  # ----------------------------------------------------------
  # Collapse multiple CpGs into one row per PWAS protein/gene
  # This creates the table style you showed.
  # ----------------------------------------------------------
  
  final_table <- overlap_long[
    ,
    .(
      CpG_ID = collapse_unique(CpG_ID),
      n_matched_CpG = uniqueN(CpG_ID),
      
      CpG_chr38 = collapse_unique(CHR_38),
      CpG_pos38 = collapse_unique(P0_38),
      
      MWAS_min_P = min(MWAS.P, na.rm = TRUE),
      MWAS_min_FDR = min(MWAS.FDR, na.rm = TRUE),
      MWAS_best_CpG = CpG_ID[which.min(MWAS.P)][1],
      MWAS_best_Z = MWAS.Z[which.min(MWAS.P)][1],
      MWAS_best_P = MWAS.P[which.min(MWAS.P)][1],
      MWAS_best_FDR = MWAS.FDR[which.min(MWAS.P)][1]
    ),
    by = .(
      pop,
      ID,
      TargetFullName,
      Target,
      EntrezGeneSymbol,
      PWAS.Z,
      PWAS.P,
      MatchedGene
    )
  ]
  
  final_table <- final_table[
    order(pop, PWAS.P, MWAS_min_P)
  ]
  
  # Put columns in the exact style you want
  final_table <- final_table[
    ,
    .(
      pop,
      ID,
      TargetFullName,
      Target,
      EntrezGeneSymbol,
      PWAS.Z,
      PWAS.P,
      MatchedGene,
      CpG_ID,
      n_matched_CpG,
      MWAS_best_CpG,
      MWAS_best_Z,
      MWAS_best_P,
      MWAS_best_FDR,
      MWAS_min_P,
      MWAS_min_FDR,
      CpG_chr38,
      CpG_pos38
    )
  ]
  
  # ----------------------------------------------------------
  # Write per-population table
  # ----------------------------------------------------------
  
  pop_out_file <- file.path(
    out_dir,
    paste0(pop, "_PWAS_table_with_matched_CpGs.txt")
  )
  
  fwrite(
    final_table,
    pop_out_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  
  # ----------------------------------------------------------
  # Summary
  # ----------------------------------------------------------
  
  summary_dt <- data.table(
    pop = pop,
    status = "ok",
    n_pwas_total = nrow(pwas),
    n_pwas_with_gene_symbol = sum(
      !is.na(pwas_anno$EntrezGeneSymbol) &
        pwas_anno$EntrezGeneSymbol != "" &
        pwas_anno$EntrezGeneSymbol != "."
    ),
    n_mwas_sig_cpg = uniqueN(mwas_cpg$CpG_ID),
    n_mwas_unique_nearest_gene = uniqueN(mwas_long$MatchedGene),
    n_output_rows = nrow(final_table),
    n_unique_matched_protein = uniqueN(final_table$ID),
    n_unique_matched_gene = uniqueN(final_table$MatchedGene),
    n_unique_matched_cpg = uniqueN(overlap_long$CpG_ID),
    n_pwas_p_lt_0_05 = sum(final_table$PWAS.P < 0.05, na.rm = TRUE),
    min_pwas_p = min(final_table$PWAS.P, na.rm = TRUE),
    min_mwas_p = min(final_table$MWAS_min_P, na.rm = TRUE),
    output_file = pop_out_file
  )
  
  cat("Output rows:", nrow(final_table), "\n")
  cat("Matched proteins:", uniqueN(final_table$ID), "\n")
  cat("Matched genes:", uniqueN(final_table$MatchedGene), "\n")
  cat("Matched CpGs:", uniqueN(overlap_long$CpG_ID), "\n")
  cat("Output:", pop_out_file, "\n")
  
  return(list(result = final_table, summary = summary_dt))
}

# ============================================================
# 2. LOAD SOMASCAN ANNOTATION
# ============================================================

ann <- fread(ann_file)

required_ann_cols <- c("SomaId", "TargetFullName", "Target", "EntrezGeneSymbol")
missing_ann_cols <- setdiff(required_ann_cols, colnames(ann))

if (length(missing_ann_cols) > 0) {
  stop("Missing annotation columns: ", paste(missing_ann_cols, collapse = ", "))
}

ann <- ann[, ..required_ann_cols]
ann <- unique(ann, by = "SomaId")

# ============================================================
# 3. RUN ALL POPULATIONS
# ============================================================

res_list <- list()
summary_list <- list()

for (i in seq_len(nrow(pwas_files))) {
  
  tmp <- process_one_pop(
    pop = pwas_files$pop[i],
    pwas_file = pwas_files$pwas_file[i],
    ann = ann
  )
  
  res_list[[pwas_files$pop[i]]] <- tmp$result
  summary_list[[pwas_files$pop[i]]] <- tmp$summary
}

# ============================================================
# 4. COMBINED TABLE
# ============================================================

combined_table <- rbindlist(res_list, fill = TRUE)

combined_file <- file.path(
  out_dir,
  "ALL_PWAS_table_with_matched_CpGs.txt"
)

fwrite(
  combined_table,
  combined_file,
  sep = ",",
  quote = FALSE,
  na = "NA"
)

# ============================================================
# 5. SUMMARY TABLE
# ============================================================

summary_dt <- rbindlist(summary_list, fill = TRUE)

summary_file <- file.path(
  out_dir,
  "summary_PWAS_table_with_matched_CpGs.txt"
)

fwrite(
  summary_dt,
  summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\n============================================================\n")
cat("Done.\n")
cat("Main table:\n")
cat(combined_file, "\n\n")
cat("Summary table:\n")
cat(summary_file, "\n")
cat("Output directory:\n")
cat(out_dir, "\n")
cat("============================================================\n\n")

print(summary_dt)