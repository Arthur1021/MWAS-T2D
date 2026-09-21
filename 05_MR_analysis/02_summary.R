library(data.table)
library(dplyr)
library(stringr)

races <- c("AA", "CA", "EA", "HA")

for (race in races) {
    
    files <- Sys.glob(paste0("results/", race, "/*_mr.txt"))
    res_list <- list()

    if (length(files) == 0) {
        message("No MR files found for ", race)
        next
    }

    for (file in files) {
        tmp <- fread(file)

        ID <- gsub("_mr.txt", "", basename(file))

        if ("Inverse variance weighted" %in% tmp$method) {
            tmp_res <- tmp[method == "Inverse variance weighted"]
        } else if ("Wald ratio" %in% tmp$method) {
            tmp_res <- tmp[method == "Wald ratio"]
        } else {
            message("No IVW or Wald ratio result for: ", ID)
            next
        }

        tmp_res$ID <- ID
        res_list[[length(res_list) + 1]] <- tmp_res
    }

    if (length(res_list) == 0) {
        message("No valid MR results for ", race)
        next
    }

    res <- rbindlist(res_list, fill = TRUE)

    res <- res[order(pval)]
    res[, FDR := p.adjust(pval, method = "fdr")]

    # original MWAS result
    mwas <- fread(paste0(
        "/mnt/Data1/sliu7/data/sliu/project/MESA_methylation_WGS_data/",
        "10_mwas_on_T2D_rsq_0.01/summary_mwas_on_T2D_",
        race,
        ".txt"
    ))

    mwas <- mwas %>% select(ID, TWAS.Z)

    res <- left_join(res, mwas, by = "ID")
    res <- as.data.table(res)

    res[, consistent := sign(b) == sign(TWAS.Z)]

    fwrite(res, paste0("result_", race, "_fdr.txt"), sep = ",")

    sig_n1 <- sum(res$FDR < 0.05, na.rm = TRUE)
    sig_n2 <- sum(res$FDR < 0.05 & res$consistent, na.rm = TRUE)
    total_n <- nrow(res)
    ratio <- sig_n2 / total_n

    print(paste(race, sig_n1, sig_n2, total_n, round(ratio, 2)))
}