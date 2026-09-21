library(data.table)

## ===============================
## 0. Parameters
## ===============================

cis_window <- 1000000   # +/- 1Mb around CpG site

rsid_file <- "/mnt/Data1/sliu7/data/sliu/database/dbSNP/GRCh38_to_rsID.txt"

base_mwas_dir <- "/mnt/Data1/sliu7/data/sliu/project/MESA_methylation_WGS_data/10_mwas_on_T2D_rsq_0.01"

base_raw_dir <- "/mnt/Data1/sliu7/project/MWAS_T2D/coloc_analysis/raw_GWAS_meta"

base_out_dir <- "/mnt/Data1/sliu7/project/MWAS_T2D/coloc_analysis/01_extract_sig_cis_region"

## MWAS ancestry -> raw meQTL folder mapping
## CA MWAS corresponds to JA raw meQTL folder
ancestry_map <- data.table(
  mwas_anc = c("AA", "CA", "EA", "HA"),
  raw_anc  = c("AA", "JA", "EA", "HA")
)
ancestry_map <- data.table(
  mwas_anc = c( "EA", "HA", "AA"),
  raw_anc  = c( "EA", "HA", "AA")
)
ancestry_map <- data.table(
  mwas_anc = c("AA"),
  raw_anc  = c("AA")
)
## ===============================
## 1. Load rsID mapping
## ===============================

cat("Loading rsID mapping...\n")

rsid <- fread(
  rsid_file,
  header = FALSE,
  col.names = c("SNP_chrpos", "rsID")
)

## Use index instead of setkey; lighter than sorting full table
setindex(rsid, SNP_chrpos)

cat("rsID mapping loaded:", nrow(rsid), "rows\n")

## ===============================
## 2. Loop through ancestries
## ===============================

for (a in seq_len(nrow(ancestry_map))) {

  mwas_anc <- ancestry_map[a, mwas_anc]
  raw_anc  <- ancestry_map[a, raw_anc]

  cat("\n====================================\n")
  cat("Processing ancestry:", mwas_anc, "\n")
  cat("Raw meQTL folder:", raw_anc, "\n")
  cat("====================================\n")

  mwas_file <- file.path(
    base_mwas_dir,
    paste0("summary_mwas_on_T2D_", mwas_anc, ".txt")
  )

  raw_dir <- file.path(base_raw_dir, raw_anc)

  out_dir <- file.path(base_out_dir, mwas_anc)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (!file.exists(mwas_file)) {
    cat("MWAS file missing:", mwas_file, "\n")
    next
  }

  if (!dir.exists(raw_dir)) {
    cat("Raw meQTL directory missing:", raw_dir, "\n")
    next
  }

  ## ===============================
  ## 2.1 Load significant CpGs
  ## ===============================

  cat("Loading MWAS file:", mwas_file, "\n")

  mwas <- fread(mwas_file)

  if (!all(c("TWAS.P", "CHR_38", "P0_38", "ID") %in% colnames(mwas))) {
    cat("Required columns missing in MWAS file:", mwas_file, "\n")
    cat("Available columns:\n")
    print(colnames(mwas))
    next
  }

  sig_mwas <- mwas[TWAS.P < 0.05 / nrow(mwas)]

  cat("Total CpGs in MWAS:", nrow(mwas), "\n")
  cat("Significant CpGs:", nrow(sig_mwas), "\n")

  if (nrow(sig_mwas) == 0) {
    cat("No significant CpGs for", mwas_anc, ". Skipping.\n")
    next
  }

  ## ===============================
  ## 2.2 Loop through significant CpGs
  ## ===============================

  for (r in seq_len(nrow(sig_mwas))) {

    chr <- sig_mwas[r, CHR_38]
    cpg <- sig_mwas[r, ID]
    cpg_pos <- sig_mwas[r, P0_38] + 1

    start_pos <- max(0, cpg_pos - cis_window)
    end_pos <- cpg_pos + cis_window

    ## ===============================
    ## Ancestry-specific meQTL file path
    ##
    ## CA:
    ## raw_GWAS_meta/JA/chr<chr>/clean.metal.<CpG>.1.tbl
    ##
    ## AA, EA, HA:
    ## raw_GWAS_meta/<ANC>/<chr>/metal.<CpG>.1.tbl
    ## ===============================

    if (mwas_anc == "CA") {
      meQTL_file <- file.path(
        raw_dir,
        paste0("chr", chr),
        paste0("clean.metal.", cpg, ".1.tbl")
      )
    } else {
      meQTL_file <- file.path(
        raw_dir,
        paste0('metal.chr', chr),
        paste0("metal.", cpg, ".1.tbl")
      )
    }

    out_file <- file.path(
      out_dir,
      paste0(cpg, ".txt")
    )

    ## Skip if output already exists
    if (file.exists(out_file) && file.info(out_file)$size > 0) {
      cat("Already exists, skipping:", mwas_anc, cpg, "\n")
      next
    }

    if (!file.exists(meQTL_file)) {
      cat("Missing file:", meQTL_file, "\n")
      next
    }

    cat("\nProcessing:", mwas_anc, cpg, "chr", chr, "\n")
    cat("  meQTL file:", meQTL_file, "\n")
    tmp <- fread(meQTL_file)
    if (ancestry_map[a,1] == 'EA'){
          tmp$N <- 4370
    }else if(ancestry_map[a,1] == 'HA'){
          tmp$N <- 1200
    }else if(ancestry_map[a,1] == 'AA'){
          tmp$N <- 6805
    }


    required_cols <- c(
      "MarkerName",
      "Allele1",
      "Allele2",
      "Effect",
      "StdErr",
      "P-value",
      "N"
    )

    if (!all(required_cols %in% colnames(tmp))) {
      cat("  Required columns missing in meQTL file:", meQTL_file, "\n")
      cat("  Available columns:\n")
      print(colnames(tmp))
      next
    }

    ## ===============================
    ## Extract chr:pos from MarkerName
    ## Example:
    ## chr11:1987706:C:T -> 11:1987706
    ## 11:1987706:C:T    -> 11:1987706
    ## ===============================

    tmp[, SNP_chrpos := sub("^chr", "", MarkerName)]
    tmp[, SNP_chrpos := sub("^([^:]+:[^:]+):.*$", "\\1", SNP_chrpos)]

    ## Split chr and pos
    tmp[, snp_chr := sub(":.*$", "", SNP_chrpos)]
    tmp[, pos := as.integer(sub("^.*:", "", SNP_chrpos))]

    ## ===============================
    ## Keep only cis-region SNPs
    ## ===============================

    tmp <- tmp[
      snp_chr == as.character(chr) &
        pos >= start_pos &
        pos <= end_pos
    ]

    cat("  SNPs in cis-region:", nrow(tmp), "\n")

    if (nrow(tmp) == 0) {
      cat("  No SNPs in cis-region. Skipping.\n")
      next
    }

    ## ===============================
    ## Add rsID
    ## ===============================

    tmp[, rsID := rsid[tmp, on = "SNP_chrpos", rsID]]

    n_total <- nrow(tmp)
    n_match <- sum(!is.na(tmp$rsID))

    cat(
      "  Matched rsID:", n_match, "/", n_total,
      "(", round(n_match / n_total * 100, 2), "% )\n"
    )

    ## Keep only SNPs with rsID
    tmp <- tmp[!is.na(rsID)]

    if (nrow(tmp) == 0) {
      cat("  No matched rsID. Skipping.\n")
      next
    }

    ## ===============================
    ## Convert to coloc input format
    ## ===============================

    tmp_coloc <- tmp[, .(
      SNP = rsID,
      chr = snp_chr,
      pos = pos,
      effect_allele = toupper(Allele1),
      other_allele = toupper(Allele2),
      beta = Effect,
      se = StdErr,
      pval = `P-value`,
      samplesize = N
    )]

    ## ===============================
    ## Remove invalid rows
    ## ===============================

    tmp_coloc <- tmp_coloc[
      !is.na(SNP) &
        !is.na(chr) &
        !is.na(pos) &
        effect_allele %in% c("A", "T", "C", "G") &
        other_allele %in% c("A", "T", "C", "G") &
        !is.na(beta) &
        !is.na(se) &
        !is.na(pval) &
        !is.na(samplesize) &
        se > 0
    ]

    ## Remove duplicated rsIDs
    tmp_coloc <- tmp_coloc[!duplicated(SNP)]

    cat("  Final SNPs for coloc:", nrow(tmp_coloc), "\n")

    if (nrow(tmp_coloc) == 0) {
      cat("  No valid SNPs after formatting. Skipping.\n")
      next
    }

    ## ===============================
    ## Save coloc-ready xQTL file
    ## ===============================

    fwrite(tmp_coloc, out_file, sep = "\t")

    cat("  Saved:", out_file, "\n")
  }

  cat("\nFinished ancestry:", mwas_anc, "\n")
}

cat("\nAll ancestries finished!\n")