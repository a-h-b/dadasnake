log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))


print("reading DADA2 classifier result")
for(dt in 1:length(snakemake@input)){
   currName <- gsub(".RDS$","",gsub("sequenceTables/tax.dada.","",snakemake@input[[dt]]))
   dadTax <- readRDS(snakemake@input[[dt]])
   dadTax <- dadTax[,-ncol(dadTax)]
   colnames(dadTax) <- paste0(colnames(dadTax),".dada2",currName)
   dadTax$OTU <- rownames(dadTax)
   if(!"lastTab" %in% ls()) lastTab <- dadTax else lastTab <- merge(lastTab,dadTax,by="OTU", all=T)
   lastTab[is.na(lastTab)] <- ""
}
print("Saving table with all taxonomies")
saveRDS(lastTab,snakemake@output[[1]])


