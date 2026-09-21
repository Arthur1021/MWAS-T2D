library(future.apply)

# 设置并行核数
plan(multisession, workers =30)

BASE_DIR <- "/mnt/Data1/sliu7/project/MWAS_T2D/coloc_analysis"

XQTL_DIR <- file.path(BASE_DIR, "01_extract_sig_cis_region")
GWAS_DIR <- file.path(BASE_DIR, "02_prepare_GWAS")
OUT_DIR  <- file.path(BASE_DIR, "03_coloc")

dir.create(OUT_DIR, showWarnings = FALSE)

# ancestry → GWAS mapping
gwas_map <- list(
  AA = "AFR",
  CA = "EAS",
  EA = "EUR",
  HA = "AMR"
)

# 收集所有任务
tasks <- list()

for (anc in names(gwas_map)) {
  
  xqtl_path <- file.path(XQTL_DIR, anc)
  gwas_file <- file.path(
    GWAS_DIR,
    paste0("T2D_2023_GWAS_summary_", gwas_map[[anc]], ".txt")
  )
  
  files <- list.files(xqtl_path, pattern = "\\.txt$", full.names = TRUE)
  
  for (f in files) {
    tasks[[length(tasks) + 1]] <- list(
      anc = anc,
      xqtl = f,
      gwas = gwas_file
    )
  }
}

# 并行执行
future_lapply(tasks, function(task) {
  
  anc <- task$anc
  xqtl <- task$xqtl
  gwas <- task$gwas
  
  basename <- tools::file_path_sans_ext(basename(xqtl))
  
  out_dir <- file.path(OUT_DIR, anc)
  dir.create(out_dir, showWarnings = FALSE)
  
  cmd <- paste(
    "Rscript coloc_pipeline.R",
    "--xQTL", xqtl,
    "--GWAS", gwas,
    "--outcome_fil", paste0(basename, "_T2D"),
    "--output_dir", out_dir
  )
  
  message(cmd)
  system(cmd)
})
