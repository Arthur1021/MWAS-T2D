library(TwoSampleMR)
library(data.table)
library(future.apply)

plan(multisession, workers = 20)

# ==============================
# ✅ PATH SETUP
# ==============================

BASE <- "/mnt/Data1/sliu7/project/MWAS_T2D"

XQTL_BASE <- file.path(BASE, "coloc_analysis/01_extract_sig_cis_region")
GWAS_BASE <- file.path(BASE, "conditional_test/01_format_gwas_summary")
OUT_BASE  <- file.path(BASE, "Mendelian_randomization/results")

LD_REF_DIR <- "/mnt/Data1/sliu7/database/1kg.v3"
PLINK_BIN  <- "plink"

dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

# ==============================
# ✅ ANCESTRY MAP
# ==============================

gwas_map <- list(
  AA = "AFR",
  CA = "EAS",
  EA = "EUR",
  HA = "AMR"
)

ld_map <- list(
  AA = "AFR",
  CA = "EAS",
  EA = "EUR",
  HA = "AMR"
)

# ==============================
# ✅ BUILD TASK LIST
# ==============================

tasks <- list()

for (anc in names(gwas_map)) {
  
  xqtl_dir <- file.path(XQTL_BASE, anc)
  gwas_file <- file.path(
    GWAS_BASE,
    paste0("T2D_2023_GWAS_summary_", gwas_map[[anc]], ".txt")
  )
  
  files <- list.files(xqtl_dir, pattern = "\\.txt$", full.names = TRUE)
  
  for (f in files) {
    tasks[[length(tasks) + 1]] <- list(
      anc = anc,
      xqtl = f,
      gwas = gwas_file
    )
  }
}

cat("✅ Total tasks:", length(tasks), "\n")

# ==============================
# ✅ MR FUNCTION
# ==============================

run_mr <- function(task) {
  
  anc <- task$anc
  xqtl_file <- task$xqtl
  gwas_file <- task$gwas
  
  cpg <- tools::file_path_sans_ext(basename(xqtl_file))
  
  message("Running: ", anc, " - ", cpg)
  
  out_dir <- file.path(OUT_BASE, anc)
  dir.create(out_dir, showWarnings = FALSE)
  
  outfile <- file.path(out_dir, paste0(cpg, "_mr.txt"))
  if (file.exists(outfile)) return(NULL)
  
  # ==============================
  # ✅ LOAD DATA
  # ==============================
  exp_raw <- fread(xqtl_file, data.table = FALSE)
  gwas_raw <- fread(gwas_file, data.table = FALSE)
  
  # ==============================
  # ✅ FORMAT EXPOSURE
  # ==============================
  exp_dat <- tryCatch({
    
    format_data(
      exp_raw,
      type = "exposure",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "pval",
      samplesize_col = "samplesize"
    )
    
  }, error = function(e) {
    message("❌ format exposure failed: ", cpg)
    return(NULL)
  })
  
  if (is.null(exp_dat) || nrow(exp_dat) == 0) return(NULL)
  
  # ✅ loose threshold
  exp_dat <- subset(exp_dat, pval.exposure < 1e-3)
  if (nrow(exp_dat) < 1) return(NULL)
  
  # ==============================
  # ✅ LOCAL LD CLUMPING
  # ==============================
  ref_prefix <- file.path(LD_REF_DIR, ld_map[[anc]])
  
  exp_dat <- tryCatch({
    
    clump_data(
      exp_dat,
      clump_kb = 10000,
      clump_r2 = 0.001,
      bfile = ref_prefix,
      plink_bin = PLINK_BIN
    )
    
  }, error = function(e) {
    message("⚠️ Clumping failed → skip: ", cpg)
    return(NULL)
  })
  
  if (is.null(exp_dat) || nrow(exp_dat) < 1) return(NULL)
  
  # ==============================
  # ✅ FORMAT OUTCOME
  # ==============================
  out_dat_full <- tryCatch({
    
    format_data(
      gwas_raw,
      type = "outcome",
      snp_col = "SNP",
      beta_col = "b",
      se_col = "se",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      eaf_col = "freq",
      pval_col = "p",
      samplesize_col = "N"
    )
    
  }, error = function(e) {
    message("❌ format outcome failed: ", cpg)
    return(NULL)
  })
  
  if (is.null(out_dat_full)) return(NULL)
  
  # ✅ overlap SNP
  out_dat <- out_dat_full[out_dat_full$SNP %in% exp_dat$SNP, ]
  if (nrow(out_dat) < 1) return(NULL)
  
  # ==============================
  # ✅ HARMONISE
  # ==============================
  dat_h <- tryCatch({
    harmonise_data(exp_dat, out_dat)
  }, error = function(e) {
    message("❌ harmonisation failed: ", cpg)
    return(NULL)
  })
  
  if (is.null(dat_h) || nrow(dat_h) < 1) return(NULL)
  
  # ==============================
  # ✅ MR ( allow 1 SNP)
  # ==============================
  res <- tryCatch({
    
    mr(
      dat_h,
      method_list = c(
        "mr_ivw",
        "mr_wald_ratio",
        "mr_egger_regression"
      )
    )
    
  }, error = function(e) {
    message("❌ MR failed: ", cpg)
    return(NULL)
  })
  
  if (is.null(res)) return(NULL)
  
  # ==============================
  # ✅ SENSITIVITY
  # ==============================
  het <- tryCatch(mr_heterogeneity(dat_h), error = function(e) NULL)
  pleio <- tryCatch(mr_pleiotropy_test(dat_h), error = function(e) NULL)
  loo <- tryCatch(mr_leaveoneout(dat_h), error = function(e) NULL)
  
  # ==============================
  # ✅ SAVE
  # ==============================
  fwrite(res, file.path(out_dir, paste0(cpg, "_mr.txt")))
  
  if (!is.null(het))
    fwrite(het, file.path(out_dir, paste0(cpg, "_het.txt")))
  
  if (!is.null(pleio))
    fwrite(pleio, file.path(out_dir, paste0(cpg, "_pleio.txt")))
  
  if (!is.null(loo))
    fwrite(loo, file.path(out_dir, paste0(cpg, "_loo.txt")))
  
  return(NULL)
}

# ==============================
# ✅ RUN
# ==============================

future_lapply(tasks, run_mr)