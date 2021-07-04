log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))


print("reading basta result")
currName <- paste0("blast.",snakemake@config[["blast"]][["tax_db"]])
basta <- read.delim(snakemake@input[[1]],stringsAsFactors=F,headers=F)  
if(ncol(basta) == 2) colnames(basta) <- c("OTU","taxonomy") else colnames(basta) <- c("OTU","taxonomy", "bestHit")
rownames(basta) <- basta$OTU
tmpTax <- data.frame(do.call(rbind, strsplit(basta$taxonomy, ";", fixed=TRUE)),stringsAsFactors = F)
rankNum <- 7
if(ncol(tmpTax)>= rankNum+1){
   dimCut <- rankNum+1 
   basta <- cbind(basta,tmpTax[,1:dimCut])
}else{
   dimCut <- ncol(tmpTax) 
   basta <- cbind(basta,tmpTax[,1:dimCut],
                    matrix("",nrow=nrow(basta),
                           ncol=rankNum+1-dimCut,
                           dimnames=list(basta$OTU,
                                     paste0("X",c((dimCut+1):(rankNum+1))))))
}
colnames(basta)[grep("^X",colnames(basta))] <- paste0(c("Domain","Phylum","Class","Order","Family","Genus","Species"),
                                                            ".blast.",currName)
colnames(basta)[colnames(basta)=="taxonomy"] <- paste0("taxonomy.blast.",currName)
print("reading OTU table")
seqTab <- readRDS(snakemake@input[[2]])
seqTab <- merge(seqTab,basta,by="OTU",all.x=T)
seqTab[is.na(seqTab)] <- ""
print("Saving table with all taxonomies")
saveRDS(seqTab,snakemake@output[[1]])
write.table(seqTab,snakemake@output[[2]],,row.names=F,sep="\t",quote=F)

