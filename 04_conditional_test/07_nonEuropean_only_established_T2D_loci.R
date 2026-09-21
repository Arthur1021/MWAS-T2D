# ============================================================
# Identify known T2D loci represented only through
# non-European MWAS analyses
#
# Inputs:
# 1) Supplementary Table 6 population-sharing results
# 2) EPIC hg19 manifest
# 3) Suzuki et al. Nature 2024 supplementary workbook
#
# Main outputs:
# - 399 non-European-only CpGs
# - Number of those CpGs overlapping known T2D loci
# - Number of unique known T2D loci represented
# - Population-specific locus counts
# ============================================================


# ============================================================
# 0. LOAD PACKAGES
# ============================================================

library(data.table)
library(readxl)


# ============================================================
# 1. FILE PATHS
# ============================================================

sig_file <- "../TableS6/Supplementary_Table_6_cross_population_CpG_sharing.txt"

manifest_file <- "../../../mnt/lvm_vol_2/sliu/project/MESA_methylation_WGS_data/EPIC.hg19.manifest.txt"

st4_file <- "41586_2024_7019_MOESM3_ESM_index_ST4.xlsx"

out_prefix <- "nonEuropean_only_known_T2D_locus_overlap"


# ============================================================
# 2. READ SUPPLEMENTARY TABLE 6
# ============================================================

sig <- fread(
  sig_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

cat("\n============================================\n")
cat("SUPPLEMENTARY TABLE 6\n")
cat("============================================\n")

print(names(sig))

cat(
  "Total CpGs in Supplementary Table 6:",
  nrow(sig),
  "\n"
)


# Rename important columns
setnames(
  sig,
  old = c(
    "CpG ID",
    "Significant Ancestries"
  ),
  new = c(
    "CpG",
    "SigPop"
  )
)


# ============================================================
# 3. FUNCTION TO CHECK POPULATION MEMBERSHIP
# ============================================================

has_pop <- function(x, pop) {

  x[is.na(x)] <- ""

  grepl(
    paste0(
      "(^|,)",
      pop,
      "(,|$)"
    ),
    x
  )
}


# ============================================================
# 4. CREATE POPULATION-SPECIFIC SIGNIFICANCE FLAGS
# ============================================================

sig[, AA_sig  := has_pop(SigPop, "AA")]
sig[, AsA_sig := has_pop(SigPop, "AsA")]
sig[, EA_sig  := has_pop(SigPop, "EA")]
sig[, HA_sig  := has_pop(SigPop, "HA")]


# ============================================================
# 5. IDENTIFY NON-EUROPEAN-ONLY CpGs
# ============================================================

nonEA <- sig[
  EA_sig == FALSE &
    (
      AA_sig == TRUE |
      AsA_sig == TRUE |
      HA_sig == TRUE
    )
]


cat("\n============================================\n")
cat("NON-EUROPEAN-ONLY CpGs\n")
cat("============================================\n")

cat(
  "Total non-European-only CpGs:",
  nrow(nonEA),
  "\n"
)

# Expected:
# 399


# ============================================================
# 6. VERIFY POPULATION BREAKDOWN
# ============================================================

cat("\nExact significance patterns:\n")

pattern_table <- nonEA[
  ,
  .N,
  by = SigPop
][order(-N)]

print(pattern_table)


n_AA_only <- nonEA[
  SigPop == "AA",
  .N
]

n_AsA_only <- nonEA[
  SigPop == "AsA",
  .N
]

n_HA_only <- nonEA[
  SigPop == "HA",
  .N
]

n_shared_nonEA <- nonEA[
  grepl(",", SigPop),
  .N
]


cat("\nAA only:", n_AA_only, "\n")
cat("AsA only:", n_AsA_only, "\n")
cat("HA only:", n_HA_only, "\n")
cat(
  "Shared across >=2 non-European populations:",
  n_shared_nonEA,
  "\n"
)

cat(
  "Total check:",
  n_AA_only +
    n_AsA_only +
    n_HA_only +
    n_shared_nonEA,
  "\n"
)


# ============================================================
# 7. READ EPIC hg19 MANIFEST
# ============================================================

manifest <- fread(
  manifest_file,
  header = FALSE
)

setnames(
  manifest,
  c(
    "CHR",
    "CPG_START",
    "CPG_END",
    "CpG"
  )
)


manifest[
  ,
  CHR := gsub(
    "^chr",
    "",
    CHR
  )
]

manifest[, CPG_START := as.integer(CPG_START)]
manifest[, CPG_END   := as.integer(CPG_END)]


cat("\n============================================\n")
cat("EPIC HG19 MANIFEST\n")
cat("============================================\n")

cat(
  "Manifest CpGs:",
  nrow(manifest),
  "\n"
)


# ============================================================
# 8. MERGE CpGs WITH hg19 COORDINATES
# ============================================================

nonEA_coord <- merge(
  nonEA,
  manifest,
  by = "CpG",
  all.x = TRUE
)


cat("\n============================================\n")
cat("COORDINATE MATCHING\n")
cat("============================================\n")

cat(
  "Non-European-only CpGs:",
  nrow(nonEA_coord),
  "\n"
)

cat(
  "CpGs with hg19 coordinates:",
  sum(!is.na(nonEA_coord$CPG_START)),
  "\n"
)

cat(
  "CpGs missing hg19 coordinates:",
  sum(is.na(nonEA_coord$CPG_START)),
  "\n"
)


if (any(is.na(nonEA_coord$CPG_START))) {

  cat("\nMissing CpGs:\n")

  print(
    nonEA_coord[
      is.na(CPG_START),
      .(
        CpG,
        SigPop
      )
    ]
  )
}


cpg <- nonEA_coord[
  !is.na(CPG_START) &
    !is.na(CPG_END)
]

cpg[, CHR := as.character(CHR)]


# ============================================================
# 9. IDENTIFY ST4 SHEET AUTOMATICALLY
# ============================================================

sheet_names <- excel_sheets(
  st4_file
)


cat("\n============================================\n")
cat("EXCEL SHEETS\n")
cat("============================================\n")

print(sheet_names)


st4_sheet <- grep(
  "^ST4$|ST4|Supplementary Table 4|Table 4",
  sheet_names,
  value = TRUE,
  ignore.case = TRUE
)


if (length(st4_sheet) == 0) {

  stop(
    "\nCould not automatically find ST4 sheet.\n",
    "Available sheets:\n",
    paste(
      sheet_names,
      collapse = "\n"
    )
  )
}


st4_sheet <- st4_sheet[1]


cat(
  "\nSelected ST4 sheet:",
  st4_sheet,
  "\n"
)


# ============================================================
# 10. READ RAW ST4
# ============================================================

st4_raw <- read_excel(
  st4_file,
  sheet = st4_sheet,
  col_names = FALSE
)


cat("\n============================================\n")
cat("RAW ST4 PREVIEW\n")
cat("============================================\n")

print(
  head(
    st4_raw,
    20
  )
)


# ============================================================
# 11. FIND HEADER ROW CONTAINING "LOCUS"
# ============================================================

header_row_candidates <- which(
  apply(
    st4_raw,
    1,
    function(x) {

      x <- trimws(
        as.character(x)
      )

      any(
        x == "Locus",
        na.rm = TRUE
      )
    }
  )
)


if (length(header_row_candidates) == 0) {

  stop(
    "\nCould not find a header row containing 'Locus'.\n"
  )
}


header_row <- header_row_candidates[1]


cat(
  "\nDetected ST4 header row:",
  header_row,
  "\n"
)


# ============================================================
# 12. RE-READ ST4 USING CORRECT HEADER
# ============================================================

st4 <- as.data.table(
  read_excel(
    st4_file,
    sheet = st4_sheet,
    skip = header_row - 1
  )
)


# ============================================================
# 13. CLEAN COLUMN NAMES
# ============================================================

clean_names <- trimws(
  gsub(
    "\\s+",
    " ",
    names(st4)
  )
)

setnames(
  st4,
  clean_names
)


cat("\n============================================\n")
cat("ST4 COLUMN NAMES\n")
cat("============================================\n")

print(names(st4))


# ============================================================
# 14. DETECT REQUIRED ST4 COLUMNS
# ============================================================

locus_col <- grep(
  "^Locus$",
  names(st4),
  value = TRUE,
  ignore.case = TRUE
)[1]


chr_col <- grep(
  "^Chromosome$",
  names(st4),
  value = TRUE,
  ignore.case = TRUE
)[1]


interval_col <- grep(
  "Interval.*b37",
  names(st4),
  value = TRUE,
  ignore.case = TRUE
)[1]


index_snv_col <- grep(
  "^Index SNV$",
  names(st4),
  value = TRUE,
  ignore.case = TRUE
)[1]


cat("\nDetected ST4 columns:\n")
cat("Locus:", locus_col, "\n")
cat("Chromosome:", chr_col, "\n")
cat("Interval:", interval_col, "\n")
cat("Index SNV:", index_snv_col, "\n")


if (
  is.na(locus_col) ||
    is.na(chr_col) ||
    is.na(interval_col)
) {

  stop(
    "\nRequired ST4 locus columns could not be identified."
  )
}


# ============================================================
# 15. RETAIN ONE LOCUS-DEFINING ROW PER T2D LOCUS
# ============================================================

loci <- st4[
  !is.na(get(locus_col)) &
    !is.na(get(chr_col)) &
    !is.na(get(interval_col))
]


cat("\n============================================\n")
cat("LOCUS-DEFINING ROWS\n")
cat("============================================\n")

cat(
  "Rows with locus + chromosome + interval:",
  nrow(loci),
  "\n"
)

# Expected:
# 611


# ============================================================
# 16. RENAME IMPORTANT COLUMNS
# ============================================================

setnames(
  loci,
  old = c(
    locus_col,
    chr_col,
    interval_col
  ),
  new = c(
    "Locus",
    "CHR",
    "Interval"
  )
)


if (!is.na(index_snv_col)) {

  index_now <- grep(
    "^Index SNV$",
    names(loci),
    value = TRUE,
    ignore.case = TRUE
  )[1]

  if (!is.na(index_now)) {

    setnames(
      loci,
      index_now,
      "Index_SNV"
    )
  }
}


# ============================================================
# 17. CLEAN T2D LOCUS INTERVALS
# ============================================================

loci[
  ,
  CHR := as.character(CHR)
]

loci[
  ,
  CHR := gsub(
    "^chr",
    "",
    CHR
  )
]


loci[
  ,
  interval_clean :=
    gsub(
      ",",
      "",
      as.character(Interval)
    )
]


loci[
  ,
  interval_clean :=
    gsub(
      "\\s+",
      "",
      interval_clean
    )
]


loci[
  ,
  c(
    "LOCUS_START",
    "LOCUS_END"
  ) :=
    tstrsplit(
      interval_clean,
      "-",
      fixed = TRUE
    )
]


loci[, LOCUS_START := as.integer(LOCUS_START)]
loci[, LOCUS_END   := as.integer(LOCUS_END)]


# ============================================================
# 18. CHECK FOR BAD INTERVALS
# ============================================================

bad_interval <- loci[
  is.na(LOCUS_START) |
    is.na(LOCUS_END)
]


if (nrow(bad_interval) > 0) {

  cat(
    "\nWARNING: Some intervals could not be parsed:\n"
  )

  print(
    bad_interval[
      ,
      .(
        Locus,
        CHR,
        Interval
      )
    ]
  )
}


loci <- loci[
  !is.na(LOCUS_START) &
    !is.na(LOCUS_END)
]


# ============================================================
# 19. CREATE UNIQUE LOCUS ID
# ============================================================

loci[
  ,
  locus_id :=
    paste0(
      CHR,
      ":",
      LOCUS_START,
      "-",
      LOCUS_END
    )
]


loci_unique <- unique(
  loci,
  by = "locus_id"
)


cat("\n============================================\n")
cat("KNOWN T2D LOCI\n")
cat("============================================\n")

cat(
  "Unique Suzuki T2D loci:",
  nrow(loci_unique),
  "\n"
)

# Expected:
# 611


cat("\nFirst 10 loci:\n")

print(
  loci_unique[
    1:min(10, .N),
    .(
      Locus,
      CHR,
      LOCUS_START,
      LOCUS_END,
      locus_id
    )
  ]
)


# ============================================================
# 20. PREPARE CpG AND T2D LOCUS INTERVALS
# ============================================================

cpg[, CPG_START := as.integer(CPG_START)]
cpg[, CPG_END   := as.integer(CPG_END)]
cpg[, CHR       := as.character(CHR)]


loci_unique[, LOCUS_START := as.integer(LOCUS_START)]
loci_unique[, LOCUS_END   := as.integer(LOCUS_END)]
loci_unique[, CHR         := as.character(CHR)]


# ============================================================
# 21. OVERLAP CpGs WITH KNOWN T2D LOCI
# ============================================================

setkey(
  loci_unique,
  CHR,
  LOCUS_START,
  LOCUS_END
)


hit <- foverlaps(
  x = cpg,
  y = loci_unique,
  by.x = c(
    "CHR",
    "CPG_START",
    "CPG_END"
  ),
  by.y = c(
    "CHR",
    "LOCUS_START",
    "LOCUS_END"
  ),
  type = "any",
  nomatch = 0L
)


# ============================================================
# 22. MAIN REVIEWER NUMBERS
# ============================================================

n_total_nonEA <- uniqueN(
  cpg$CpG
)

n_cpg_at_T2D_loci <- uniqueN(
  hit$CpG
)

n_T2D_loci <- uniqueN(
  hit$locus_id
)


pct_cpg_at_T2D_loci <- round(
  100 *
    n_cpg_at_T2D_loci /
    n_total_nonEA,
  2
)


cat("\n\n")
cat("============================================================\n")
cat("MAIN RESULTS\n")
cat("============================================================\n")

cat(
  "Total non-European-only CpGs:",
  n_total_nonEA,
  "\n"
)

cat(
  "Non-European-only CpGs overlapping known T2D loci:",
  n_cpg_at_T2D_loci,
  "\n"
)

cat(
  "Percent of non-European-only CpGs overlapping known T2D loci:",
  pct_cpg_at_T2D_loci,
  "%\n"
)

cat(
  "Unique known T2D loci represented:",
  n_T2D_loci,
  "\n"
)


# ============================================================
# 23. POPULATION-SPECIFIC COUNTS
# ============================================================

AA_hit <- hit[
  AA_sig == TRUE
]

AsA_hit <- hit[
  AsA_sig == TRUE
]

HA_hit <- hit[
  HA_sig == TRUE
]


cat("\n============================================================\n")
cat("POPULATION-SPECIFIC RESULTS\n")
cat("============================================================\n")


cat(
  "AA CpGs overlapping known T2D loci:",
  uniqueN(AA_hit$CpG),
  "\n"
)

cat(
  "AA unique known T2D loci:",
  uniqueN(AA_hit$locus_id),
  "\n\n"
)


cat(
  "AsA CpGs overlapping known T2D loci:",
  uniqueN(AsA_hit$CpG),
  "\n"
)

cat(
  "AsA unique known T2D loci:",
  uniqueN(AsA_hit$locus_id),
  "\n\n"
)


cat(
  "HA CpGs overlapping known T2D loci:",
  uniqueN(HA_hit$CpG),
  "\n"
)

cat(
  "HA unique known T2D loci:",
  uniqueN(HA_hit$locus_id),
  "\n"
)


# ============================================================
# 24. EXACT SIGNIFICANCE-PATTERN SUMMARY
# ============================================================

pattern_overlap <- hit[
  ,
  .(
    n_CpGs = uniqueN(CpG),
    n_T2D_loci = uniqueN(locus_id)
  ),
  by = SigPop
][order(-n_CpGs)]


cat("\n============================================================\n")
cat("SIGNIFICANCE PATTERNS AMONG OVERLAPPING CpGs\n")
cat("============================================================\n")

print(pattern_overlap)


# ============================================================
# 25. LOCUS-LEVEL SUMMARY
# ============================================================

locus_summary <- hit[
  ,
  .(
    n_nonEA_CpGs = uniqueN(CpG),

    CpGs = paste(
      sort(
        unique(CpG)
      ),
      collapse = ","
    ),

    SigPop_patterns = paste(
      sort(
        unique(SigPop)
      ),
      collapse = ";"
    ),

    AA = any(AA_sig),

    AsA = any(AsA_sig),

    HA = any(HA_sig)
  ),
  by = .(
    locus_id,
    Locus,
    CHR,
    LOCUS_START,
    LOCUS_END
  )
]


setorderv(
  locus_summary,
  c(
    "CHR",
    "LOCUS_START"
  )
)


cat(
  "\nUnique known T2D loci in locus summary:",
  nrow(locus_summary),
  "\n"
)


# ============================================================
# 26. CpGs NOT OVERLAPPING KNOWN T2D LOCI
# ============================================================

hit_cpgs <- unique(
  hit$CpG
)


not_at_T2D_locus <- cpg[
  !CpG %in% hit_cpgs
]


cat(
  "\nNon-European-only CpGs NOT overlapping known T2D loci:",
  uniqueN(not_at_T2D_locus$CpG),
  "\n"
)


# ============================================================
# 27. EXPORT 399 NON-EUROPEAN-ONLY CpGs
# ============================================================

fwrite(
  nonEA_coord,
  paste0(
    out_prefix,
    "_399_nonEuropean_only_CpGs.txt"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 28. EXPORT CpG-LOCUS OVERLAP
# ============================================================

fwrite(
  hit,
  paste0(
    out_prefix,
    "_CpG_known_T2D_locus_overlap.txt"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 29. EXPORT UNIQUE T2D LOCUS SUMMARY
# ============================================================

fwrite(
  locus_summary,
  paste0(
    out_prefix,
    "_unique_known_T2D_loci.txt"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 30. EXPORT CpGs OUTSIDE KNOWN T2D LOCI
# ============================================================

fwrite(
  not_at_T2D_locus,
  paste0(
    out_prefix,
    "_CpGs_not_at_known_T2D_loci.txt"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 31. EXPORT POPULATION PATTERN SUMMARY
# ============================================================

fwrite(
  pattern_overlap,
  paste0(
    out_prefix,
    "_population_pattern_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 32. CREATE REVIEWER SUMMARY TABLE
# ============================================================

reviewer_summary <- data.table(
  Metric = c(
    "Non-European-only CpG-T2D associations",
    "CpGs overlapping known T2D loci",
    "Percent overlapping known T2D loci",
    "Unique known T2D loci represented",
    "AA CpGs overlapping known T2D loci",
    "AA unique known T2D loci",
    "AsA CpGs overlapping known T2D loci",
    "AsA unique known T2D loci",
    "HA CpGs overlapping known T2D loci",
    "HA unique known T2D loci"
  ),

  Value = c(
    n_total_nonEA,
    n_cpg_at_T2D_loci,
    pct_cpg_at_T2D_loci,
    n_T2D_loci,
    uniqueN(AA_hit$CpG),
    uniqueN(AA_hit$locus_id),
    uniqueN(AsA_hit$CpG),
    uniqueN(AsA_hit$locus_id),
    uniqueN(HA_hit$CpG),
    uniqueN(HA_hit$locus_id)
  )
)


fwrite(
  reviewer_summary,
  paste0(
    out_prefix,
    "_reviewer_summary.txt"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 33. FINAL SUMMARY FOR REVIEWER
# ============================================================

cat("\n\n")
cat("============================================================\n")
cat("FINAL SUMMARY FOR REVIEWER\n")
cat("============================================================\n")

cat(
  "Non-European-only CpG-T2D associations:",
  n_total_nonEA,
  "\n"
)

cat(
  "CpGs overlapping known T2D loci:",
  n_cpg_at_T2D_loci,
  "\n"
)

cat(
  "Percentage:",
  pct_cpg_at_T2D_loci,
  "%\n"
)

cat(
  "Unique known T2D loci represented:",
  n_T2D_loci,
  "\n"
)


cat("\nPopulation-specific results:\n")

cat(
  "AA:",
  uniqueN(AA_hit$CpG),
  "CpGs /",
  uniqueN(AA_hit$locus_id),
  "known T2D loci\n"
)

cat(
  "AsA:",
  uniqueN(AsA_hit$CpG),
  "CpGs /",
  uniqueN(AsA_hit$locus_id),
  "known T2D loci\n"
)

cat(
  "HA:",
  uniqueN(HA_hit$CpG),
  "CpGs /",
  uniqueN(HA_hit$locus_id),
  "known T2D loci\n"
)


cat("============================================================\n")