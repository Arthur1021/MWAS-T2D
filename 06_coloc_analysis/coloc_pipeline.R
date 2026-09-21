suppressMessages(library(data.table))
suppressMessages(library(valr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(coloc))
suppressMessages(library(optparse))

## ============================================================
## 0. Options
## ============================================================

option_list = list(
  make_option("--xQTL", action = "store", default = NA, type = "character",
              help = "Path to xQTL summary file"),

  make_option("--xQTL_type", action = "store", default = "quant", type = "character",
              help = "xQTL trait type for coloc: quant or cc. Default: quant"),

  make_option("--GWAS", action = "store", default = NA, type = "character",
              help = "Path to disease GWAS summary file"),

  make_option("--GWAS_type", action = "store", default = "cc", type = "character",
              help = "GWAS trait type for coloc: quant or cc. Default: cc"),

  make_option("--outcome_file", action = "store", default = NA, type = "character",
              help = "Output prefix / target name"),

  make_option("--p12", action = "store", default = 1e-5, type = "numeric",
              help = "Prior probability that a SNP is associated with both traits. Default: 1e-5"),

  make_option("--output_dir", action = "store", default = NA, type = "character",
              help = "Path to output directory"),

  make_option("--window", action = "store", default = 250000, type = "integer",
              help = "Window size around significant xQTL SNPs. Default: 250000"),

  make_option("--xqtl_p_threshold", action = "store", default = 1e-5, type = "numeric",
              help = "P-value threshold to define significant xQTL loci. Default: 1e-5"),

  make_option("--min_snps", action = "store", default = 50, type = "integer",
              help = "Minimum number of overlapping SNPs required for coloc. Default: 50"),

  make_option("--gwas_case_fraction", action = "store", default = NA, type = "numeric",
              help = "Case fraction for case-control GWAS, required if GWAS_type=cc and no s column exists")
)

opt = parse_args(OptionParser(option_list = option_list))

## ============================================================
## 1. Basic input checks
## ============================================================

if (is.na(opt$xQTL) || !file.exists(opt$xQTL)) {
  stop("xQTL file not found: ", opt$xQTL)
}

if (is.na(opt$GWAS) || !file.exists(opt$GWAS)) {
  stop("GWAS file not found: ", opt$GWAS)
}

if (is.na(opt$output_dir)) {
  stop("Please provide --output_dir")
}

if (is.na(opt$outcome_file)) {
  stop("Please provide --outcome_file")
}

if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

## ============================================================
## 2. Helper functions
## ============================================================

is_palindromic <- function(a1, a2) {
  pair <- paste0(a1, a2)
  pair %in% c("AT", "TA", "CG", "GC")
}

standardize_chr <- function(chr) {
  chr <- as.character(chr)
  chr <- gsub("^chr", "", chr, ignore.case = TRUE)
  return(chr)
}

## ============================================================
## 3. Load disease GWAS
## Required columns:
## SNP, chr, pos, effect_allele, other_allele, beta, se
## Optional:
## samplesize, eaf, ncase, ncontrol, s
## ============================================================

cat("Loading GWAS file...\n")

gwas <- fread(opt$GWAS, data.table = FALSE)

required_gwas_cols <- c("SNP", "chr", "pos", "effect_allele", "other_allele", "beta", "se")
missing_gwas_cols <- setdiff(required_gwas_cols, colnames(gwas))

if (length(missing_gwas_cols) > 0) {
  stop("GWAS file is missing required columns: ", paste(missing_gwas_cols, collapse = ", "))
}

gwas <- gwas %>%
  filter(
    effect_allele %in% c("A", "T", "C", "G"),
    other_allele %in% c("A", "T", "C", "G"),
    str_starts(SNP, "rs")
  )

gwas$chr <- standardize_chr(gwas$chr)
gwas$pos <- as.numeric(gwas$pos)
gwas$beta <- as.numeric(gwas$beta)
gwas$se <- as.numeric(gwas$se)

if (!"samplesize" %in% colnames(gwas)) {
  gwas$samplesize <- NA_real_
} else {
  gwas$samplesize <- as.numeric(gwas$samplesize)
}

## If case/control counts exist, calculate sample size and case fraction
if (all(c("ncase", "ncontrol") %in% colnames(gwas))) {
  gwas$ncase <- as.numeric(gwas$ncase)
  gwas$ncontrol <- as.numeric(gwas$ncontrol)
  gwas$samplesize <- gwas$ncase + gwas$ncontrol
  gwas$s <- gwas$ncase / (gwas$ncase + gwas$ncontrol)
}

## If no s column and user provides case fraction, use it
if (!"s" %in% colnames(gwas)) {
  gwas$s <- NA_real_
}

if (opt$GWAS_type == "cc") {
  if (all(is.na(gwas$s)) && !is.na(opt$gwas_case_fraction)) {
    gwas$s <- opt$gwas_case_fraction
  }
}

cat("GWAS SNPs retained:", nrow(gwas), "\n")

## ============================================================
## 4. Load xQTL
## Required columns:
## SNP, chr, pos, effect_allele, other_allele, beta, se, pval, samplesize
## Optional:
## eaf
## ============================================================

cat("Loading xQTL file...\n")

xqtl <- fread(opt$xQTL, data.table = FALSE)

required_xqtl_cols <- c("SNP", "chr", "pos", "effect_allele", "other_allele", "beta", "se", "pval", "samplesize")
missing_xqtl_cols <- setdiff(required_xqtl_cols, colnames(xqtl))

if (length(missing_xqtl_cols) > 0) {
  stop("xQTL file is missing required columns: ", paste(missing_xqtl_cols, collapse = ", "))
}

xqtl <- xqtl %>%
  filter(
    effect_allele %in% c("A", "T", "C", "G"),
    other_allele %in% c("A", "T", "C", "G"),
    str_starts(SNP, "rs")
  )

xqtl$chr <- standardize_chr(xqtl$chr)
xqtl$pos <- as.numeric(xqtl$pos)
xqtl$beta <- as.numeric(xqtl$beta)
xqtl$se <- as.numeric(xqtl$se)
xqtl$pval <- as.numeric(xqtl$pval)
xqtl$samplesize <- as.numeric(xqtl$samplesize)

cat("xQTL SNPs retained:", nrow(xqtl), "\n")

## ============================================================
## 5. Define significant xQTL regions
## ============================================================

xqtl_sig <- xqtl %>%
  filter(pval < opt$xqtl_p_threshold)

cat("Significant xQTL SNPs:", nrow(xqtl_sig), "\n")

if (nrow(xqtl_sig) == 0) {
  stop("No significant xQTL SNPs found using threshold: ", opt$xqtl_p_threshold)
}

sig_region <- data.frame(
  chrom = as.numeric(xqtl_sig$chr),
  start = as.numeric(xqtl_sig$pos - opt$window),
  end = as.numeric(xqtl_sig$pos + opt$window)
)

sig_region$start[sig_region$start < 0] <- 0

unique_region <- data.frame(bed_merge(sig_region))

cat("Merged coloc regions:", nrow(unique_region), "\n")

## ============================================================
## 6. Run coloc region by region
## ============================================================

for (l in seq_len(nrow(unique_region))) {

  chr_l <- as.character(unique_region[l, ]$chrom)
  start_l <- as.numeric(unique_region[l, ]$start)
  end_l <- as.numeric(unique_region[l, ]$end)

  index <- paste0(chr_l, "_", start_l, "_", end_l)

  cat("\n==============================\n")
  cat("Processing region:", index, "\n")
  cat("==============================\n")

  out_snp_file <- file.path(opt$output_dir, paste0(opt$outcome_file, "_", index, "_SNP.txt"))
  out_coloc_file <- file.path(opt$output_dir, paste0(opt$outcome_file, "_", index, "_coloc.txt"))

  ## Skip if output already exists and is non-empty
  if (file.exists(out_coloc_file) && file.info(out_coloc_file)$size > 0) {
    cat("Result already exists. Skipping:", out_coloc_file, "\n")
    next
  }

  tmp_xqtl <- xqtl %>%
    filter(chr == chr_l, pos >= start_l, pos <= end_l)

  if (nrow(tmp_xqtl) == 0) {
    cat("No xQTL SNPs in region. Skipping.\n")
    next
  }

  tmp_xqtl_gwas <- inner_join(
    tmp_xqtl,
    gwas,
    by = "SNP",
    suffix = c("_xqtl", "_gwas")
  )

  cat("Overlapping SNPs before allele filtering:", nrow(tmp_xqtl_gwas), "\n")

  if (nrow(tmp_xqtl_gwas) < opt$min_snps) {
    cat("Too few overlapping SNPs before filtering. Skipping.\n")
    next
  }

  ## Remove palindromic SNPs
  tmp_xqtl_gwas <- tmp_xqtl_gwas[
    !is_palindromic(
      tmp_xqtl_gwas$effect_allele_xqtl,
      tmp_xqtl_gwas$other_allele_xqtl
    ),
  ]

  cat("SNPs after removing palindromic SNPs:", nrow(tmp_xqtl_gwas), "\n")

  if (nrow(tmp_xqtl_gwas) < opt$min_snps) {
    cat("Too few SNPs after removing palindromic SNPs. Skipping.\n")
    next
  }

  ## Keep only allele-matched SNPs
  tmp_xqtl_gwas <- tmp_xqtl_gwas[
    (
      tmp_xqtl_gwas$effect_allele_xqtl == tmp_xqtl_gwas$effect_allele_gwas &
        tmp_xqtl_gwas$other_allele_xqtl == tmp_xqtl_gwas$other_allele_gwas
    ) |
      (
        tmp_xqtl_gwas$effect_allele_xqtl == tmp_xqtl_gwas$other_allele_gwas &
          tmp_xqtl_gwas$other_allele_xqtl == tmp_xqtl_gwas$effect_allele_gwas
      ),
  ]

  cat("SNPs after allele matching:", nrow(tmp_xqtl_gwas), "\n")

  if (nrow(tmp_xqtl_gwas) < opt$min_snps) {
    cat("Too few SNPs after allele matching. Skipping.\n")
    next
  }

  ## Flip GWAS beta if GWAS effect allele does not match xQTL effect allele
  tmp_xqtl_gwas$beta_gwas <- ifelse(
    tmp_xqtl_gwas$effect_allele_xqtl == tmp_xqtl_gwas$effect_allele_gwas,
    tmp_xqtl_gwas$beta_gwas,
    -1 * tmp_xqtl_gwas$beta_gwas
  )

  ## Compute varbeta
  tmp_xqtl_gwas$varbeta_xqtl <- tmp_xqtl_gwas$se_xqtl^2
  tmp_xqtl_gwas$varbeta_gwas <- tmp_xqtl_gwas$se_gwas^2

  ## Remove duplicated SNPs
  tmp_xqtl_gwas <- tmp_xqtl_gwas[!duplicated(tmp_xqtl_gwas$SNP), ]

  ## ==========================================================
  ## 6.1 Prepare xQTL MAF
  ## Prefer xQTL EAF. If absent, use GWAS EAF if available.
  ## ==========================================================

  if ("eaf_xqtl" %in% colnames(tmp_xqtl_gwas)) {
    tmp_xqtl_gwas$eaf_for_maf <- as.numeric(tmp_xqtl_gwas$eaf_xqtl)
  } else if ("eaf" %in% colnames(tmp_xqtl_gwas)) {
    tmp_xqtl_gwas$eaf_for_maf <- as.numeric(tmp_xqtl_gwas$eaf)
  } else if ("eaf_gwas" %in% colnames(tmp_xqtl_gwas)) {
    tmp_xqtl_gwas$eaf_for_maf <- as.numeric(tmp_xqtl_gwas$eaf_gwas)
  } else {
    stop("No eaf, eaf_xqtl, or eaf_gwas column found. coloc.abf with beta/varbeta needs MAF.")
  }

  tmp_xqtl_gwas$MAF <- pmin(
    tmp_xqtl_gwas$eaf_for_maf,
    1 - tmp_xqtl_gwas$eaf_for_maf
  )

  ## ==========================================================
  ## 6.2 Remove invalid rows
  ## ==========================================================

  keep <- complete.cases(
    tmp_xqtl_gwas$beta_xqtl,
    tmp_xqtl_gwas$varbeta_xqtl,
    tmp_xqtl_gwas$beta_gwas,
    tmp_xqtl_gwas$varbeta_gwas,
    tmp_xqtl_gwas$MAF,
    tmp_xqtl_gwas$samplesize_xqtl
  )

  tmp_xqtl_gwas <- tmp_xqtl_gwas[keep, ]

  tmp_xqtl_gwas <- tmp_xqtl_gwas[
    tmp_xqtl_gwas$varbeta_xqtl > 0 &
      tmp_xqtl_gwas$varbeta_gwas > 0 &
      tmp_xqtl_gwas$MAF > 0 &
      tmp_xqtl_gwas$MAF < 0.5 &
      tmp_xqtl_gwas$samplesize_xqtl > 0,
  ]

  cat("SNPs after missing/invalid filtering:", nrow(tmp_xqtl_gwas), "\n")

  if (nrow(tmp_xqtl_gwas) < opt$min_snps) {
    cat("Too few valid SNPs for coloc. Skipping.\n")
    next
  }

  ## ==========================================================
  ## 6.3 Prepare coloc datasets
  ## ==========================================================

  data_xqtl <- data.frame(
    beta = tmp_xqtl_gwas$beta_xqtl,
    varbeta = tmp_xqtl_gwas$varbeta_xqtl,
    snp = tmp_xqtl_gwas$SNP,
    MAF = tmp_xqtl_gwas$MAF,
    N = tmp_xqtl_gwas$samplesize_xqtl
  )

  data_gwas <- data.frame(
    beta = tmp_xqtl_gwas$beta_gwas,
    varbeta = tmp_xqtl_gwas$varbeta_gwas,
    snp = tmp_xqtl_gwas$SNP
  )

  ## Add GWAS N if available
  if ("samplesize_gwas" %in% colnames(tmp_xqtl_gwas)) {
    if (!all(is.na(tmp_xqtl_gwas$samplesize_gwas))) {
      data_gwas$N <- tmp_xqtl_gwas$samplesize_gwas
    }
  } else if ("samplesize" %in% colnames(tmp_xqtl_gwas)) {
    if (!all(is.na(tmp_xqtl_gwas$samplesize))) {
      data_gwas$N <- tmp_xqtl_gwas$samplesize
    }
  }

  ## Add case fraction for case-control GWAS
  if (opt$GWAS_type == "cc") {
    if ("s_gwas" %in% colnames(tmp_xqtl_gwas)) {
      s_val <- unique(na.omit(tmp_xqtl_gwas$s_gwas))
      if (length(s_val) > 0) {
        data_gwas$s <- s_val[1]
      }
    } else if ("s" %in% colnames(tmp_xqtl_gwas)) {
      s_val <- unique(na.omit(tmp_xqtl_gwas$s))
      if (length(s_val) > 0) {
        data_gwas$s <- s_val[1]
      }
    } else if (!is.na(opt$gwas_case_fraction)) {
      data_gwas$s <- opt$gwas_case_fraction
    }
  }

  ## Convert to list for coloc
  data_xqtl <- as.list(data_xqtl)
  data_gwas <- as.list(data_gwas)

  data_xqtl$type <- opt$xQTL_type
  data_gwas$type <- opt$GWAS_type

  ## For case-control GWAS, coloc strongly prefers N and s
  if (opt$GWAS_type == "cc") {
    if (is.null(data_gwas$N)) {
      cat("Warning: GWAS N is missing for case-control coloc.\n")
    }
    if (is.null(data_gwas$s)) {
      cat("Warning: GWAS case fraction s is missing for case-control coloc.\n")
    }
  }

  ## ==========================================================
  ## 6.4 Run coloc
  ## ==========================================================

  cat("Running coloc.abf for region:", index, "\n")

  res <- tryCatch(
    {
      coloc.abf(
        dataset1 = data_xqtl,
        dataset2 = data_gwas,
        p12 = opt$p12
      )
    },
    error = function(e) {
      cat("coloc failed for region:", index, "\n")
      cat("Error:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(res)) {
    next
  }

  res_result <- as.data.frame(res$result)

  res_summary_df <- as.data.frame(
    matrix(res$summary, ncol = length(res$summary), byrow = TRUE)
  )
  colnames(res_summary_df) <- names(res$summary)

  target_id <- first(str_split(opt$outcome_file, "_")[[1]])

  res_summary_df$target <- target_id
  res_summary_df$outcome_file <- opt$outcome_file
  res_summary_df$index <- index
  res_summary_df$chr <- chr_l
  res_summary_df$start <- start_l
  res_summary_df$end <- end_l
  res_summary_df$n_snps_used <- nrow(tmp_xqtl_gwas)
  res_summary_df$p12 <- opt$p12

  ## Save SNP-level result and region-level coloc summary
  write.table(
    res_result,
    out_snp_file,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )

  write.table(
    res_summary_df,
    out_coloc_file,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )

  cat("Saved SNP-level result:", out_snp_file, "\n")
  cat("Saved coloc summary:", out_coloc_file, "\n")
}

cat("\nAll coloc regions finished!\n")