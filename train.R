# Package
library(BEDMatrix)
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
library(optparse)

###########
# Options#
###########

# Options
option_list <- list(
    make_option("--sumstats", type = "character", default = FALSE, action = "store", help = "List of sumstats"
    ),
    make_option("--reference", type = "character", default = FALSE, action = "store", help = "Path of reference panel"
    ),
    make_option("--output", type = "character", default = FALSE, action = "store", help = "Path of output folder (end it by /)"
    ),
    make_option("--method", type = "character", default = FALSE, action = "store", help = "Method being used"
    )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Pass the options to variables
sumstats  <- opt$sumstats %>% fread(., header = FALSE) %>% as.data.frame()
reference <- opt$reference
output    <- opt$output
METHOD    <- opt$method %>% strsplit(., split = " ") %>% unlist()

##############################
# Source this before train() #
##############################
source("support.R")

# Main function
train <- function(sumstats, reference) {
    ########################
    # Read reference panel #
    ########################

    # Read the bim file of reference panel
    bim.ref <- paste0(reference, ".bim") %>% fread(., showProgress = FALSE) %>% as.data.frame()
    names(bim.ref)[2] <- "SNP"

    # Extract genotype information
    seq.ref <- paste0(reference, ".bed") %>% BEDMatrix(., simple_names = TRUE)

    for (qwe in 1:nrow(sumstats)) {
        ###########
        # Read ss #
        ###########

        ss <- fread(sumstats$V1[qwe], showProgress = FALSE) %>% as.data.frame()

        # Rename some columns
        names(ss)[names(ss) == "MarkerName"] <- "SNP"
        names(ss)[names(ss) == "Allele1"] <- "A1"
        names(ss)[names(ss) == "Allele2"] <- "A2"

        # UPPERCASE the a1/a2
        ss$A1 <- toupper(ss$A1)
        ss$A2 <- toupper(ss$A2)

        # Compute Z-score
        ss["Zscore"] <- ss$Effect / ss$StdErr

        ############
        # Quick QC #
        ############

        # Stop if ss is empty to begin with
        if (nrow(ss) == 0) {
            cat("The summary statistics file is empty. \n")
            next
        }

        # Only keep the SNPs that are IN reference panel
        ss <- ss[ss$SNP %in% bim.ref$SNP, ]

        # Sieve out problematic SNPs
        #list.1  <- bim.ref$V2[duplicated(bim.ref$V2)]
        #list.2  <- bim.ref$V2[nchar(bim.ref$V5) > 1 | nchar(bim.ref$V6) > 1]
        list.3  <- ss$SNP[duplicated(ss$SNP)]
        #problem <- ss$SNP %in% list.1 | ss$SNP %in% list.2 | ss$SNP %in% list.3
        problem <- ss$SNP %in% list.3

        ss <- ss[!problem, ]

        #rm(list.1)
        #rm(list.2)
        rm(list.3)
        rm(problem)

        # Stop if ss has few row left
        if (nrow(ss) <= 1) {
            cat("Too few SNPs left after initial QC. \n")
            next
        }

        #################################################
        # Preprocess genotype matrix of reference panel #
        #################################################

        genotype.ref <- seq.ref[, ss$SNP]

        # Stop if genotype is empty
        if (ncol(genotype.ref) == 0) {
            cat("The summary statistics and the reference panel have no intersection. \n")
            next
        }

        # Compare to see if all (A1 and A2) pairs are well-aligned
        # Always use reference panel's A1/A2 as the standard order
        bim.temp <- subset(bim.ref, SNP %in% ss$SNP)
        ss.temp  <- left_join(ss, bim.temp, by = "SNP")
        problem  <- !(ss.temp$A1 == ss.temp$V5)

        if (sum(problem) != 0) {
            # Flip the problematic pairs
            ss$Zscore[problem] <- -1 * ss$Zscore[problem]

            ss$A1 <- ss.temp$V5
            ss$A2 <- ss.temp$V6
        }

        rm(bim.temp)
        rm(ss.temp)
        rm(problem)

        # Patch up the NAs and centralize genotype matrix
        genotype.ref <- PatchUp(genotype.ref)
        genotype.ref <- Standardize(genotype.ref)

        #####################
        # Compute LD matrix #
        #####################

        matrix.LD <- eigenMapMatMult(t(genotype.ref), genotype.ref) / nrow(genotype.ref)

        ############
        # r vector #
        ############

        # Compute r vector
        r <- ss$Zscore / sqrt(ss$N - 1 + ss$Zscore ^2)

        # Get size
        size <- nrow(ss)

        ##########################
        # Iterate By parameter s #
        ##########################

        result <- list()

        s <- 0

        for (method in METHOD) {
            #########################################
            # Constructing the initial lambda array #
            #########################################

            # Get the maximum for lambda
            z.temp <- numeric(size)
            for (m in 1:size) {
                z.temp[m] <- abs(r[m])
            }

            lambda.max   <- max(z.temp)
            lambda.min   <- lambda.max * 1E-3
            lambda.array <- exp(1) ^ seq(log(lambda.max), log(lambda.min), length = 100)
            rm(z.temp)

            #####################
            # Iterate by lambda #
            #####################

            beta <- t(matrix(0, nrow = 1, ncol = size))

            max.iteration <- 250

            threshold <- 1e-5
            alpha     <- 0.5

            if (method == "SCAD") {
                gamma <- 3.7
            } else {
                gamma <- 3
            }

            if (method == "MCP") {
                res <- MCP(r, matrix.LD, lambda.array, s, gamma, max.iteration, threshold)
            }
            if (method == "LASSO") {
                res <- ElNet(r, matrix.LD, lambda.array, s, 1, max.iteration, threshold)
            }
            if (method == "ElNet") {
                res <- ElNet(r, matrix.LD, lambda.array, s, 0.5, max.iteration, threshold)
            }
            if (method == "MNet") {
                res <- MNet(r, matrix.LD, lambda.array, s, alpha, gamma, max.iteration, threshold)
            }
            if (method == "SCAD") {
                res <- SCAD(r, matrix.LD, lambda.array, s, gamma, max.iteration, threshold)
            }
        }

            ############################################################################################################################
            #         -----.                   osssso            ---------..`                                                          #
            #        `ssssso                   osssso           `sssssssssssss+:`                                                      #
            #        `ssssso       `...`       osssso   ``      `ssssssssssssssss/        `...`                ```          `...`      #
            #        `ssssso   `/osssssss+:`   osssso:osssso/`  `sssss/  `.:osssss+   `:+sssssss+:`   +sssso.+sssso/`    -+sssssss+-   #
            #        `ssssso  /sssssssssssso-  ossssssssssssss- `sssss/     `ssssss` :sssssssssssss-  +ssssssssssssso  .ossso/::ossso` #
            #        `ssssso -sssss:``./sssss` osssso-``./sssss``sssss/      osssss..sssss:``./sssss. +sssss-``osssss` ossss/---:ssss+ #
            #        .ssssso /ssss+    `sssss. ossss/    `sssss.`sssss/    `/ssssso :sssso     sssss- +sssss   /sssss``sssssssssssssso #
            # .:/oso/osssss/ .sssss+--:osssso  osssss/--:osssso `ssssso+++ossssss+` `sssss+:-:+sssso` +sssss   /sssss` +ssss/`  `-`    #
            # .ossssssssss+`  .osssssssssss+`  ossssosssssssso. `ssssssssssssss+-    .osssssssssss+.  +sssss   /sssss` `/sssssssssso/` #
            #   -/+osso+/.      ./+oosoo+:`    /++++.-/+oo+/-   `++++++++//:-.`        .:+oosoo+:.    :+++++   :+++++`   `:/oosoo+/-`  #
            ############################################################################################################################

        write.table(
            res,
            file = paste0(sumstats$V1[qwe], ".model"),
            row.names = TRUE,
            col.names = TRUE,
            quote = FALSE
        )
    }

    return()
}

# Call train()
train(sumstats, reference)
