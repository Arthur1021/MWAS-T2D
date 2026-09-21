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
  sig_independent <- all_filter[all_filter$ID %in% df_sig$ID, ]$ID
  all_filter <- df_sig[!df_sig$ID %in% all_filter$ID, ] #   non independent
  all_filter <- select(all_filter, ID, CHR_38, P0_38, NSNP, NWGT, MODELCV.R2, TWAS.Z, TWAS.P)
  all_filter <- left_join(all_filter, anno, by = c('CHR_38' = 'Chr', 'P0_38' = 'Start'))
  cat("non Significant after conditional MWAS:", nrow(all_filter), "\n")

  ## ===============================
  ## Save output
  ## ===============================

  out_file <- file.path(
    out_dir,
    paste0("conditional_non_significant_T2D_", model_anc, "_", gwas_anc, ".txt")
  )

  fwrite(all_filter, out_file, sep = ",")

  cat("Saved:", out_file, "\n")
}

cat("\nAll ancestries finished!\n")


## ============================================================
## 3. Cross-population analysis:
##    established/non-independent methylation associations
##    detected only through non-European populations
## ============================================================

cat("\n==============================\n")
cat("Cross-population established-locus analysis\n")
cat("==============================\n")

## Read all non-independent / established-locus CpG files
nonind_files <- list.files(
  out_dir,
  pattern = "^conditional_non_significant_T2D_.*\\.txt$",
  full.names = TRUE
)

if (length(nonind_files) == 0) {
  stop("No non-independent conditional output files found.")
}

nonind_all <- rbindlist(
  lapply(nonind_files, function(f) {
    x <- fread(f)
    
    ## If model_anc was not saved inside previous loop, recover it from file name
    if (!"model_anc" %in% colnames(x)) {
      fname <- basename(f)
      x$model_anc <- sub("^conditional_non_significant_T2D_([^_]+)_.*$", "\\1", fname)
    }
    
    x
  }),
  fill = TRUE
)

## Standardize ancestry label
nonind_all[, ancestry := fifelse(model_anc == "AA", "AA",
                          fifelse(model_anc == "CA", "AsA",
                          fifelse(model_anc == "EA", "EA",
                          fifelse(model_anc == "HA", "HA", model_anc))))]

cat("Total non-independent CpG rows:", nrow(nonind_all), "\n")
cat("Unique non-independent CpGs:", uniqueN(nonind_all$ID), "\n")

## ============================================================
## 4. CpG-level cross-population sharing
## ============================================================

cpg_summary <- nonind_all[
  ,
  .(
    Significant_Ancestries = paste(sort(unique(ancestry)), collapse = ","),
    Any_AA  = any(ancestry == "AA"),
    Any_AsA = any(ancestry == "AsA"),
    Any_EA  = any(ancestry == "EA"),
    Any_HA  = any(ancestry == "HA"),
    CHR_38 = CHR_38[1],
    P0_38 = P0_38[1],
    Gene.refGene = Gene.refGene[1]
  ),
  by = ID
]

cpg_summary[
  ,
  Any_non_EA := Any_AA | Any_AsA | Any_HA
]

cpg_summary[
  ,
  Non_European_only_established_CpG := Any_non_EA & !Any_EA
]

non_ea_only_cpgs <- cpg_summary[
  Non_European_only_established_CpG == TRUE
]

cat("\nNumber of established-locus CpG associations detected only through non-European analyses:\n")
cat(nrow(non_ea_only_cpgs), "\n")

## ============================================================
## 5. Locus / nearest-gene-level summary
## ============================================================

## Clean gene/locus labels
non_ea_only_cpgs[
  ,
  locus_label := Gene.refGene
]

non_ea_only_cpgs[
  is.na(locus_label) | locus_label == "",
  locus_label := paste0("chr", CHR_38, ":", P0_38)
]

## Collapse by nearest-gene locus label
locus_summary <- non_ea_only_cpgs[
  ,
  .(
    Number_of_CpGs = uniqueN(ID),
    CpGs = paste(sort(unique(ID)), collapse = ","),
    Chr = paste(sort(unique(CHR_38)), collapse = ","),
    Position_range = paste0(min(P0_38, na.rm = TRUE), "-", max(P0_38, na.rm = TRUE)),
    Significant_Ancestries = paste(sort(unique(unlist(strsplit(Significant_Ancestries, ",")))), collapse = ",")
  ),
  by = locus_label
]

setorder(locus_summary, -Number_of_CpGs, locus_label)

cat("\nNumber of established T2D locus labels detected only through non-European analyses:\n")
cat(nrow(locus_summary), "\n")

## ============================================================
## 6. Optional: split semicolon-separated genes
##    This gives unique gene symbols instead of locus labels
## ============================================================

gene_split <- non_ea_only_cpgs[
  ,
  .(gene_symbol = unlist(strsplit(as.character(locus_label), ";"))),
  by = .(ID, Significant_Ancestries)
]

gene_split[, gene_symbol := trimws(gene_symbol)]
gene_split <- gene_split[!is.na(gene_symbol) & gene_symbol != ""]

gene_summary <- gene_split[
  ,
  .(
    Number_of_CpGs = uniqueN(ID),
    CpGs = paste(sort(unique(ID)), collapse = ","),
    Significant_Ancestries = paste(sort(unique(Significant_Ancestries)), collapse = ";")
  ),
  by = gene_symbol
]

setorder(gene_summary, -Number_of_CpGs, gene_symbol)

cat("\nNumber of unique gene symbols after splitting semicolon-separated gene labels:\n")
cat(nrow(gene_summary), "\n")

## ============================================================
## 7. Save output tables
## ============================================================

fwrite(
  cpg_summary,
  file.path(out_dir, "Established_locus_CpG_cross_population_summary.txt"),
  sep = "\t"
)

fwrite(
  non_ea_only_cpgs,
  file.path(out_dir, "Non_European_only_established_locus_CpGs.txt"),
  sep = "\t"
)

fwrite(
  locus_summary,
  file.path(out_dir, "Non_European_only_established_T2D_locus_labels.txt"),
  sep = "\t"
)

fwrite(
  gene_summary,
  file.path(out_dir, "Non_European_only_established_T2D_gene_symbols_split.txt"),
  sep = "\t"
)

## Also save Excel file
library(openxlsx)

xlsx_out <- file.path(
  out_dir,
  "Non_European_only_established_T2D_loci_summary.xlsx"
)

write.xlsx(
  list(
    CpG_cross_population_summary = cpg_summary,
    Non_EA_only_CpGs = non_ea_only_cpgs,
    Non_EA_only_locus_labels = locus_summary,
    Non_EA_only_gene_symbols_split = gene_summary
  ),
  file = xlsx_out,
  overwrite = TRUE
)

cat("\nSaved cross-population established-locus analysis to:\n")
cat(xlsx_out, "\n")