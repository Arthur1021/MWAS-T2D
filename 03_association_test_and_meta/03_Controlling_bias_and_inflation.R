library(data.table)
library(bacon)
library(dplyr)

# Lambda GC function
lambda_gc <- function(p){
    chisq <- qchisq(1 - p, df = 1)
    median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
}

res <- data.frame(
    rsq = numeric(),
    race = character(),
    Total_N = numeric(),
    N_sig_raw = numeric(),
    lambda_raw = numeric(),
    N_sig_bacon = numeric(),
    lambda_bacon = numeric(),
    bacon_bias = numeric(),
    bacon_inflation = numeric(),
    stringsAsFactors = FALSE
)

races <- c("AA", "CA", "EA", "HA")
rsqs <- c("0.01")

dir.create(
    "16_Controlling_bias_and_inflation",
    showWarnings = FALSE,
    recursive = TRUE
)

for (rsq in rsqs) {

    for (race in races) {

        cat("Processing:", race, "rsq =", rsq, "\n")

        infile <- paste0(
            "10_mwas_on_T2D_rsq_",
            rsq,
            "/summary_mwas_on_T2D_",
            race,
            ".txt"
        )

        df <- fread(infile, data.table = FALSE)

        # remove last column if needed
        df <- df[, 1:(ncol(df) - 1)]

        # remove NA values
        df <- df %>%
            filter(
                !is.na(TWAS.P),
                !is.na(TWAS.Z),
                is.finite(TWAS.Z)
            )

        # raw Bonferroni
        df$Bonf <- p.adjust(
            df$TWAS.P,
            method = "bonferroni"
        )

        # BACON correction
        bc <- bacon(df$TWAS.Z)

        df$P.bacon <- pval(bc)

        df$Bonf.bacon <- p.adjust(
            df$P.bacon,
            method = "bonferroni"
        )

        # save significant hits
        df_all <- filter(
            df,
            Bonf < 0.05
        )

        fwrite(
            df_all,
            paste0(
                "16_Controlling_bias_and_inflation/",
                "T2D_MWAS_",
                race,
                "_",
                rsq,
                ".txt"
            ),
            sep = "\t"
        )

        # lambda before bacon
        lambda_raw <- round(
            lambda_gc(df$TWAS.P),
            2
        )

        # lambda after bacon
        lambda_bacon <- round(
            lambda_gc(df$P.bacon),
            2
        )

        # bacon parameters
        bacon_bias <- round(
            bias(bc),
            3
        )

        bacon_inflation <- round(
            inflation(bc),
            3
        )

        # summary row
        res <- rbind(
            res,
            data.frame(
                rsq = as.numeric(rsq),
                race = race,
                Total_N = nrow(df),
                N_sig_raw = sum(df$Bonf < 0.05),
                lambda_raw = lambda_raw,
                N_sig_bacon = sum(df$Bonf.bacon < 0.05),
                lambda_bacon = lambda_bacon,
                bacon_bias = bacon_bias,
                bacon_inflation = bacon_inflation,
                stringsAsFactors = FALSE
            )
        )
    }
}

fwrite(
    res,
    "16_Controlling_bias_and_inflation/bacon_corrected_summary.txt",
    sep = "\t"
)

print(res)
