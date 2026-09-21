library(data.table)
races <- c('AA', 'EA', 'HA')
for (race in races){
    result <- data.frame()
    files <- Sys.glob(paste0('TWAS_for_nearest_genes_from_sig_CpGs/', race, '/results/*.txt'))
    for (file in files){
        tmp <- fread(file)
        result <- rbind(result, tmp)
    }
    result <- result[order(result$TWAS.P),]
    result$FDR <- p.adjust(result$TWAS.P, method = 'fdr')
    print(paste(race, 'Bonf', sum(result$TWAS.P < 0.05/nrow(result)), nrow(result), sum(result$TWAS.P < 0.05/nrow(result))/nrow(result)))
    print(paste(race, 'fdr', sum(result$FDR < 0.05), nrow(result), sum(result$FDR < 0.05)/nrow(result)))
    fwrite(result, paste0('MESA_TWAS_', race, '.txt'), sep = ',')
}