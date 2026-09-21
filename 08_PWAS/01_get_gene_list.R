library(data.table)
library(dplyr)

# ============================================================
# 0. SETTINGS
# ============================================================

work_dir <- "/mnt/Data1/sliu7/project/MWAS_T2D/PWAS_MESA_1K"
setwd(work_dir)

mwas_dir <- "/mnt/Data1/sliu7/data/sliu/project/MESA_methylation_WGS_data/10_mwas_on_T2D_rsq_0.01"

anno_file <- "/mnt/Data1/sliu7/data/sliu/project/MESA_methylation_WGS_data/18_Tables/table_anno.hg38_multianno.txt"

out_base <- file.path(work_dir, "nearest_gene_lists_from_sig_CpGs")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

ancestry_map <- data.table(
  pop = c("AA", "CA", "EA", "HA")
)

# ============================================================
# 1. LOAD CpG ANNOTATION
# ============================================================

ann <- fread(anno_file)

ann <- ann %>%
  select(Chr, Start, Gene.refGene)

colnames(ann)[3] <- "Nearest Gene"

ann <- as.data.table(ann)

# Remove duplicated genomic annotation rows if present
ann <- unique(ann, by = c("Chr", "Start"))

# Make sure chromosome format matches MWAS CHR_38
ann[, Chr := as.character(Chr)]
ann[, Chr := gsub("^chr", "", Chr, ignore.case = TRUE)]

# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

get_mwas_file <- function(pop) {
  file.path(
    mwas_dir,
    paste0("summary_mwas_on_T2D_", pop, ".txt")
  )
}

clean_gene_list <- function(x) {
  genes <- unique(unlist(strsplit(as.character(x), ";")))
  genes <- trimws(genes)
  genes <- genes[!is.na(genes)]
  genes <- genes[genes != ""]
  genes <- genes[genes != "."]
  return(sort(unique(genes)))
}

# ============================================================
# 3. MAIN LOOP: GENERATE GENE LIST ONLY
# ============================================================

summary_list <- list()
all_gene_list <- list()

for (k in seq_len(nrow(ancestry_map))) {
  
  pop <- ancestry_map$pop[k]
  
  cat("\n============================================================\n")
  cat("Processing:", pop, "\n")
  cat("============================================================\n")
  
  out_dir <- file.path(out_base, pop)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  mwas_file <- get_mwas_file(pop)
  
  if (!file.exists(mwas_file)) {
    warning("Missing MWAS file: ", mwas_file)
    
    summary_list[[length(summary_list) + 1]] <- data.table(
      pop = pop,
      mwas_file = mwas_file,
      n_total_cpg = NA_integer_,
      bonf_threshold = NA_real_,
      n_sig_cpg = NA_integer_,
      n_sig_cpg_with_gene = NA_integer_,
      n_unique_nearest_gene = NA_integer_,
      sig_cpg_anno_file = NA_character_,
      gene_list_file = NA_character_,
      status = "missing_mwas_file"
    )
    
    next
  }
  
  # ----------------------------------------------------------
  # Read MWAS
  # ----------------------------------------------------------
  
  mwas <- fread(mwas_file)
  
  required_mwas_cols <- c("ID", "CHR_38", "P0_38", "TWAS.P")
  missing_mwas_cols <- setdiff(required_mwas_cols, colnames(mwas))
  
  if (length(missing_mwas_cols) > 0) {
    warning(
      "Missing columns in MWAS file ",
      mwas_file,
      ": ",
      paste(missing_mwas_cols, collapse = ", ")
    )
    
    summary_list[[length(summary_list) + 1]] <- data.table(
      pop = pop,
      mwas_file = mwas_file,
      n_total_cpg = nrow(mwas),
      bonf_threshold = NA_real_,
      n_sig_cpg = NA_integer_,
      n_sig_cpg_with_gene = NA_integer_,
      n_unique_nearest_gene = NA_integer_,
      sig_cpg_anno_file = NA_character_,
      gene_list_file = NA_character_,
      status = paste0("missing_mwas_cols:", paste(missing_mwas_cols, collapse = ","))
    )
    
    next
  }
  
  # Make chromosome format consistent
  mwas[, CHR_38 := as.character(CHR_38)]
  mwas[, CHR_38 := gsub("^chr", "", CHR_38, ignore.case = TRUE)]
  
  # ----------------------------------------------------------
  # Select Bonferroni-significant CpGs
  # ----------------------------------------------------------
  
  bonf_threshold <- 0.05 / nrow(mwas)
  
  sig_mwas <- mwas[
    !is.na(TWAS.P) & TWAS.P < bonf_threshold
  ]
  
  sig_mwas <- unique(sig_mwas, by = "ID")
  
  cat("Total CpGs:", nrow(mwas), "\n")
  cat("Bonferroni threshold:", bonf_threshold, "\n")
  cat("Significant CpGs:", nrow(sig_mwas), "\n")
  
  if (nrow(sig_mwas) == 0) {
    
    summary_list[[length(summary_list) + 1]] <- data.table(
      pop = pop,
      mwas_file = mwas_file,
      n_total_cpg = nrow(mwas),
      bonf_threshold = bonf_threshold,
      n_sig_cpg = 0,
      n_sig_cpg_with_gene = 0,
      n_unique_nearest_gene = 0,
      sig_cpg_anno_file = NA_character_,
      gene_list_file = NA_character_,
      status = "no_bonferroni_sig_cpg"
    )
    
    next
  }
  
  # ----------------------------------------------------------
  # Add nearest gene annotation
  # ----------------------------------------------------------
  
  sig_mwas_anno <- left_join(
    as.data.frame(sig_mwas),
    as.data.frame(ann),
    by = c("CHR_38" = "Chr", "P0_38" = "Start")
  )
  
  sig_mwas_anno <- as.data.table(sig_mwas_anno)
  
  sig_cpg_anno_file <- file.path(
    out_dir,
    paste0(pop, "_Bonferroni_sig_CpGs_with_nearest_gene.txt")
  )
  
  fwrite(
    sig_mwas_anno,
    sig_cpg_anno_file,
    sep = "\t",
    quote = FALSE
  )
  
  # ----------------------------------------------------------
  # Extract nearest genes only
  # ----------------------------------------------------------
  
  gene_list <- clean_gene_list(sig_mwas_anno$`Nearest Gene`)
  
  cat("Unique nearest genes:", length(gene_list), "\n")
  
  gene_list_file <- file.path(
    out_dir,
    paste0(pop, "_nearest_gene_list_from_Bonferroni_sig_CpGs.txt")
  )
  
  fwrite(
    data.table(Gene = gene_list),
    gene_list_file,
    sep = "\t",
    quote = FALSE,
    col.names = TRUE
  )
  
  all_gene_list[[pop]] <- data.table(
    pop = pop,
    Gene = gene_list
  )
  
  summary_list[[length(summary_list) + 1]] <- data.table(
    pop = pop,
    mwas_file = mwas_file,
    n_total_cpg = nrow(mwas),
    bonf_threshold = bonf_threshold,
    n_sig_cpg = nrow(sig_mwas),
    n_sig_cpg_with_gene = sum(
      !is.na(sig_mwas_anno$`Nearest Gene`) &
        sig_mwas_anno$`Nearest Gene` != "." &
        sig_mwas_anno$`Nearest Gene` != ""
    ),
    n_unique_nearest_gene = length(gene_list),
    sig_cpg_anno_file = sig_cpg_anno_file,
    gene_list_file = gene_list_file,
    status = "ok"
  )
}

# ============================================================
# 4. WRITE COMBINED OUTPUTS
# ============================================================

summary_dt <- rbindlist(summary_list, fill = TRUE)

summary_file <- file.path(
  out_base,
  "summary_generate_nearest_gene_lists_from_sig_CpGs.txt"
)

fwrite(
  summary_dt,
  summary_file,
  sep = "\t",
  quote = FALSE
)

if (length(all_gene_list) > 0) {
  
  all_gene_dt <- rbindlist(all_gene_list, fill = TRUE)
  
  all_gene_file <- file.path(
    out_base,
    "all_pop_nearest_gene_list_from_Bonferroni_sig_CpGs.txt"
  )
  
  fwrite(
    all_gene_dt,
    all_gene_file,
    sep = "\t",
    quote = FALSE
  )
  
  shared_gene_file <- file.path(
    out_base,
    "nearest_gene_overlap_by_pop.txt"
  )
  
  gene_overlap <- all_gene_dt[
    ,
    .(
      pops = paste(sort(unique(pop)), collapse = ";"),
      n_pop = uniqueN(pop)
    ),
    by = Gene
  ][order(-n_pop, Gene)]
  
  fwrite(
    gene_overlap,
    shared_gene_file,
    sep = "\t",
    quote = FALSE
  )
}

cat("\n============================================================\n")
cat("Done. Gene lists generated only. No TWAS/FUSION was run.\n")
cat("Summary written to:\n")
cat(summary_file, "\n")
cat("Output base directory:\n")
cat(out_base, "\n")
cat("============================================================\n\n")

print(summary_dt[, .N, by = .(pop, status)])