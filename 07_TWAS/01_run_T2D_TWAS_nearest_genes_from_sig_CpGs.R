library(data.table)
library(dplyr)
library(parallel)

# ============================================================
# 0. SETTINGS
# ============================================================

work_dir <- "/mnt/Data1/sliu7/project/MWAS_T2D/MESA_TWAS"
setwd(work_dir)

fusion_script <- "/mnt/Data1/sliu7/mnt/lvm_vol_2/sliu/pipeline/TWAS_fusion/bin/FUSION.assoc_test.R"

mwas_dir <- "/mnt/Data1/sliu7/data/sliu/project/MESA_methylation_WGS_data/10_mwas_on_T2D_rsq_0.01"

anno_file <- "/mnt/Data1/sliu7/data/sliu/project/MESA_methylation_WGS_data/18_Tables/table_anno.hg38_multianno.txt"

mrna_model_base <- "/mnt/Data1/sliu7/data/sliu/database/MESA_mrna_models/mesa.mrna.fusion"

gwas_base <- "/mnt/Data1/sliu7/data/sliu/database/GWAS_Summary/T2D_2023"

ld_base <- "/mnt/Data1/sliu7/data/sliu/database/LD-reference"

out_base <- "TWAS_for_nearest_genes_from_sig_CpGs"

dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# Number of parallel chromosome jobs per ancestry
# Start with 4 or 6; increase if LD/GWAS I/O is stable.
n_cores <- 20

# ============================================================
# 1. ANCESTRY MAP
# ============================================================

ancestry_map <- data.table(
  pop = c("AA", "EA", "HA"),
  model_pop = c("AFR", "EUR", "HIS"),
  gwas_pop = c("AFR", "EUR", "AMR"),
  ld_pop = c("AFR", "EUR", "AMR")
)

# ============================================================
# 2. LOAD CpG ANNOTATION
# ============================================================

ann <- fread(anno_file)

ann <- ann %>%
  select(Chr, Start, Gene.refGene)

colnames(ann)[3] <- "Nearest Gene"

ann <- as.data.table(ann)

# Remove duplicated genomic annotation rows if present
ann <- unique(ann, by = c("Chr", "Start"))

# ============================================================
# 3. HELPER FUNCTIONS
# ============================================================

get_mwas_file <- function(pop) {
  file.path(
    mwas_dir,
    paste0("summary_mwas_on_T2D_", pop, ".txt")
  )
}

get_model_pos_file <- function(model_pop) {
  file.path(
    mrna_model_base,
    paste0(model_pop, "_MESA_rna_models.pos")
  )
}

get_weights_dir <- function(model_pop) {
  file.path(
    mrna_model_base,
    model_pop,
    "rdata"
  )
}

get_gwas_file <- function(gwas_pop) {
  file.path(
    gwas_base,
    paste0("T2D_2023_GWAS_summary_", gwas_pop, ".txt")
  )
}

get_ld_ref <- function(ld_pop) {
  file.path(
    ld_base,
    paste0(ld_pop, "_chr_hg38"),
    paste0("1000G.", ld_pop, ".ALLSNP.QC.")
  )
}

clean_gene_list <- function(x) {
  genes <- unique(unlist(strsplit(x, ";")))
  genes <- trimws(genes)
  genes <- genes[!is.na(genes)]
  genes <- genes[genes != ""]
  genes <- genes[genes != "."]
  return(unique(genes))
}

# ============================================================
# 4. FUNCTION TO RUN ONE CHROMOSOME
# ============================================================

run_one_chr <- function(chr,
                        sig_pos,
                        pop,
                        model_pop,
                        gwas_pop,
                        ld_pop,
                        fusion_script,
                        gwas_file,
                        weights_dir,
                        ref_ld_chr,
                        pos_dir,
                        result_dir,
                        log_dir) {
  
  chr_pos <- sig_pos[CHR == chr]
  
  chr_pos_file <- file.path(
    pos_dir,
    paste0(pop, "_nearest_gene_sig_models_chr", chr, ".pos")
  )
  
  fwrite(
    chr_pos[, !"Gene"],
    chr_pos_file,
    sep = "\t",
    quote = FALSE
  )
  
  out_file <- file.path(
    result_dir,
    paste0("T2D_", gwas_pop, "_", pop, "_nearest_gene_chr", chr, ".txt")
  )
  
  log_file <- file.path(
    log_dir,
    paste0("T2D_", gwas_pop, "_", pop, "_nearest_gene_chr", chr, ".log")
  )
  
  # Skip if output already exists
  if (file.exists(out_file) && file.info(out_file)$size > 0) {
    message("Skip existing: ", out_file)
    
    return(data.table(
      pop = pop,
      model_pop = model_pop,
      gwas_pop = gwas_pop,
      ld_pop = ld_pop,
      chr = chr,
      n_models = nrow(chr_pos),
      pos_file = chr_pos_file,
      out_file = out_file,
      log_file = log_file,
      status = "skipped_existing_output",
      exit_status = NA_integer_
    ))
  }
  
  cmd <- paste(
    "Rscript", shQuote(fusion_script),
    "--sumstats", shQuote(gwas_file),
    "--weights", shQuote(chr_pos_file),
    "--weights_dir", shQuote(weights_dir),
    "--ref_ld_chr", shQuote(ref_ld_chr),
    "--chr", chr,
    "--out", shQuote(out_file),
    ">", shQuote(log_file),
    "2>&1"
  )
  
  message("Running chromosome ", chr, " for ", pop)
  message(cmd)
  
  exit_status <- system(cmd)
  
  if (exit_status == 0 && file.exists(out_file) && file.info(out_file)$size > 0) {
    status <- "ok"
  } else {
    status <- paste0("failed_exit_", exit_status)
  }
  
  return(data.table(
    pop = pop,
    model_pop = model_pop,
    gwas_pop = gwas_pop,
    ld_pop = ld_pop,
    chr = chr,
    n_models = nrow(chr_pos),
    pos_file = chr_pos_file,
    out_file = out_file,
    log_file = log_file,
    status = status,
    exit_status = exit_status
  ))
}

# ============================================================
# 5. MAIN LOOP
# ============================================================

summary_list <- list()

for (k in seq_len(nrow(ancestry_map))) {
  
  pop <- ancestry_map$pop[k]
  model_pop <- ancestry_map$model_pop[k]
  gwas_pop <- ancestry_map$gwas_pop[k]
  ld_pop <- ancestry_map$ld_pop[k]
  
  cat("\n============================================================\n")
  cat("Processing:", pop, "\n")
  cat("RNA model:", model_pop, "\n")
  cat("GWAS:", gwas_pop, "\n")
  cat("LD:", ld_pop, "\n")
  cat("============================================================\n")
  
  # ----------------------------------------------------------
  # Input/output paths
  # ----------------------------------------------------------
  
  mwas_file <- get_mwas_file(pop)
  model_pos_file <- get_model_pos_file(model_pop)
  weights_dir <- get_weights_dir(model_pop)
  gwas_file <- get_gwas_file(gwas_pop)
  ref_ld_chr <- get_ld_ref(ld_pop)
  
  out_dir <- file.path(out_base, pop)
  pos_dir <- file.path(out_dir, "pos")
  log_dir <- file.path(out_dir, "logs")
  result_dir <- file.path(out_dir, "results")
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pos_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ----------------------------------------------------------
  # Check input files/directories
  # ----------------------------------------------------------
  
  if (!file.exists(mwas_file)) {
    warning("Missing MWAS file: ", mwas_file)
    summary_list[[length(summary_list) + 1]] <- data.table(
      pop = pop,
      model_pop = model_pop,
      gwas_pop = gwas_pop,
      ld_pop = ld_pop,
      chr = NA_integer_,
      n_models = NA_integer_,
      pos_file = NA_character_,
      out_file = NA_character_,
      log_file = NA_character_,
      status = "missing_mwas_file",
      exit_status = NA_integer_
    )
    next
  }
  
  if (!file.exists(model_pos_file)) {
    warning("Missing model pos file: ", model_pos_file)
    summary_list[[length(summary_list) + 1]] <- data.table(
      pop = pop,
      model_pop = model_pop,
      gwas_pop = gwas_pop,
      ld_pop = ld_pop,
      chr = NA_integer_,
      n_models = NA_integer_,
      pos_file = model_pos_file,
      out_file = NA_character_,
      log_file = NA_character_,
      status = "missing_model_pos_file",
      exit_status = NA_integer_
    )
    next
  }
  
  if (!dir.exists(weights_dir)) {
    warning("Missing weights dir: ", weights_dir)
    summary_list[[length(summary_list) + 1]] <- data.table(
      pop = pop,
      model_pop = model_pop,
      gwas_pop = gwas_pop,
      ld_pop = ld_pop,
      chr = NA_integer_,
      n_models = NA_integer_,
      pos_file = NA_character_,
      out_file = NA_character_,
      log_file = NA_character_,
      status = "missing_weights_dir",
      exit_status = NA_integer_
    )
    next
  }
  
  if (!file.exists(gwas_file)) {
    warning("Missing GWAS file: ", gwas_file)
    summary_list[[length(summary_list) + 1]] <- data.table(
      pop = pop,
      model_pop = model_pop,
      gwas_pop = gwas_pop,
      ld_pop = ld_pop,
      chr = NA_integer_,
      n_models = NA_integer_,
      pos_file = NA_character_,
      out_file = gwas_file,
      log_file = NA_character_,
      status = "missing_gwas_file",
      exit_status = NA_integer_
    )
    next
  }
  
  # ----------------------------------------------------------
  # Read MWAS and select significant CpGs
  # ----------------------------------------------------------
  
  mwas <- fread(mwas_file)
  
  required_mwas_cols <- c("ID", "CHR_38", "P0_38", "TWAS.P")
  missing_mwas_cols <- setdiff(required_mwas_cols, colnames(mwas))
  
  if (length(missing_mwas_cols) > 0) {
    warning(
      "Missing columns in MWAS file ",
      mwas_file,
      ": ",
      paste(missing_mwas_cols, collapse = ",")
    )
    
    summary_list[[length(summary_list) + 1]] <- data.table(
      pop = pop,
      model_pop = model_pop,
      gwas_pop = gwas_pop,
      ld_pop = ld_pop,
      chr = NA_integer_,
      n_models = NA_integer_,
      pos_file = NA_character_,
      out_file = NA_character_,
      log_file = NA_character_,
      status = paste0("missing_mwas_cols:", paste(missing_mwas_cols, collapse = ",")),
      exit_status = NA_integer_
    )
    next
  }
  
  bonf_threshold <- 0.05 / nrow(mwas)
  
  sig_mwas <- mwas[
    !is.na(TWAS.P) & TWAS.P < bonf_threshold
  ]
  
  sig_mwas <- unique(sig_mwas, by = "ID")
  
  cat("Bonferroni threshold:", bonf_threshold, "\n")
  cat("Significant CpGs:", nrow(sig_mwas), "\n")
  
  if (nrow(sig_mwas) == 0) {
    summary_list[[length(summary_list) + 1]] <- data.table(
      pop = pop,
      model_pop = model_pop,
      gwas_pop = gwas_pop,
      ld_pop = ld_pop,
      chr = NA_integer_,
      n_models = 0,
      pos_file = NA_character_,
      out_file = NA_character_,
      log_file = NA_character_,
      status = "no_bonferroni_sig_cpg",
      exit_status = NA_integer_
    )
    next
  }
  
  # ----------------------------------------------------------
  # Add nearest gene annotation
  # ----------------------------------------------------------
  
  sig_mwas <- left_join(
    as.data.frame(sig_mwas),
    as.data.frame(ann),
    by = c("CHR_38" = "Chr", "P0_38" = "Start")
  )
  
  sig_mwas <- as.data.table(sig_mwas)
  
  sig_mwas_anno_file <- file.path(
    out_dir,
    paste0(pop, "_Bonferroni_sig_CpGs_with_nearest_gene.txt")
  )
  
  fwrite(
    sig_mwas,
    sig_mwas_anno_file,
    sep = "\t",
    quote = FALSE
  )
  
  # ----------------------------------------------------------
  # Extract nearest genes
  # ----------------------------------------------------------
  
  gene_list <- clean_gene_list(sig_mwas$`Nearest Gene`)
  
  cat("Unique nearest genes:", length(gene_list), "\n")
  
  if (length(gene_list) == 0) {
    warning("No nearest genes found for ", pop)
    
    summary_list[[length(summary_list) + 1]] <- data.table(
      pop = pop,
      model_pop = model_pop,
      gwas_pop = gwas_pop,
      ld_pop = ld_pop,
      chr = NA_integer_,
      n_models = 0,
      pos_file = NA_character_,
      out_file = NA_character_,
      log_file = NA_character_,
      status = "no_nearest_gene",
      exit_status = NA_integer_
    )
    next
  }
  
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
  
  # ----------------------------------------------------------
  # Read RNA model position file and filter target genes
  # ----------------------------------------------------------
  
  model_pos <- fread(model_pos_file)
  
  required_pos_cols <- c("PANEL", "WGT", "ID", "CHR", "P0", "P1", "N")
  missing_pos_cols <- setdiff(required_pos_cols, colnames(model_pos))
  
  if (length(missing_pos_cols) > 0) {
    warning(
      "Missing columns in model pos file ",
      model_pos_file,
      ": ",
      paste(missing_pos_cols, collapse = ",")
    )
    
    summary_list[[length(summary_list) + 1]] <- data.table(
      pop = pop,
      model_pop = model_pop,
      gwas_pop = gwas_pop,
      ld_pop = ld_pop,
      chr = NA_integer_,
      n_models = NA_integer_,
      pos_file = model_pos_file,
      out_file = NA_character_,
      log_file = NA_character_,
      status = paste0("missing_pos_cols:", paste(missing_pos_cols, collapse = ",")),
      exit_status = NA_integer_
    )
    next
  }
  
  # Original ID may look like tissue_gene; keep gene symbol after last "_"
  model_pos[, Gene := sub(".*_", "", ID)]
  
  sig_pos <- model_pos[Gene %in% gene_list]
  
  cat("Matched RNA models:", nrow(sig_pos), "\n")
  
  if (nrow(sig_pos) == 0) {
    warning("No RNA model matched nearest genes for ", pop)
    
    summary_list[[length(summary_list) + 1]] <- data.table(
      pop = pop,
      model_pop = model_pop,
      gwas_pop = gwas_pop,
      ld_pop = ld_pop,
      chr = NA_integer_,
      n_models = 0,
      pos_file = NA_character_,
      out_file = NA_character_,
      log_file = NA_character_,
      status = "no_matched_rna_model",
      exit_status = NA_integer_
    )
    next
  }
  
  all_sig_pos_file <- file.path(
    pos_dir,
    paste0(pop, "_nearest_gene_sig_models_all_chr.pos")
  )
  
  fwrite(
    sig_pos[, !"Gene"],
    all_sig_pos_file,
    sep = "\t",
    quote = FALSE
  )
  
  # ----------------------------------------------------------
  # Parallel run FUSION by chromosome
  # ----------------------------------------------------------
  
  chrs <- sort(unique(sig_pos$CHR))
  
  cat("Chromosomes to run:", paste(chrs, collapse = ", "), "\n")
  
  use_cores <- min(n_cores, length(chrs))
  
  cat("Running chromosomes in parallel with", use_cores, "cores\n")
  
  res_list <- mclapply(
    X = chrs,
    FUN = run_one_chr,
    sig_pos = sig_pos,
    pop = pop,
    model_pop = model_pop,
    gwas_pop = gwas_pop,
    ld_pop = ld_pop,
    fusion_script = fusion_script,
    gwas_file = gwas_file,
    weights_dir = weights_dir,
    ref_ld_chr = ref_ld_chr,
    pos_dir = pos_dir,
    result_dir = result_dir,
    log_dir = log_dir,
    mc.cores = use_cores
  )
  
  ancestry_summary <- rbindlist(res_list, fill = TRUE)
  
  summary_list <- c(summary_list, res_list)
  
  ancestry_summary_file <- file.path(
    out_dir,
    paste0("run_T2D_TWAS_nearest_genes_summary_", pop, ".txt")
  )
  
  fwrite(
    ancestry_summary,
    ancestry_summary_file,
    sep = "\t",
    quote = FALSE
  )
  
  cat("Ancestry summary written to:\n")
  cat(ancestry_summary_file, "\n")
  
  print(ancestry_summary[, .N, by = status])
}

# ============================================================
# 6. WRITE FINAL SUMMARY
# ============================================================

summary_dt <- rbindlist(summary_list, fill = TRUE)

summary_file <- file.path(
  out_base,
  "run_T2D_TWAS_for_nearest_genes_from_sig_CpGs_summary.txt"
)

fwrite(
  summary_dt,
  summary_file,
  sep = "\t",
  quote = FALSE
)

cat("\n============================================================\n")
cat("Done.\n")
cat("Summary written to:\n")
cat(summary_file, "\n")
cat("============================================================\n\n")

print(summary_dt[, .N, by = .(pop, status)])