log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

if(snakemake@params[["what"]]=="ASV"){
  units <- "ASV"
}else{
  units <- "OTU"
}

print("reading basta result")
currName <- paste0("blast.",snakemake@params[["tax_db"]])
basta <- read.delim(snakemake@input[[1]],stringsAsFactors=F,header=F) 
if(units=="ASV" & any(grepl("^OTU",basta[,1]))) basta[,1] <- gsub("OTU","ASV",basta[,1])
if(ncol(basta) == 2) colnames(basta) <- c(units,"taxonomy") else colnames(basta) <- c(units,"taxonomy", "bestHit")
rownames(basta) <- basta[,which(colnames(basta)==units)]
basta$taxonomy <- gsub(";$","",sapply(basta$taxonomy,function(x) paste0(x,paste0(rep(";",7-lengths(regmatches(x, gregexpr(";", x)))),collapse="",sep=""))))
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
                           dimnames=list(basta[,which(colnames(basta)==units)],
                                     paste0("X",c((dimCut+1):(rankNum+1))))))
}
colnames(basta)[grep("^X",colnames(basta))] <- paste0(c("Domain","Phylum","Class","Order","Family","Genus","Species"),
                                                            ".blast.",currName)
colnames(basta)[colnames(basta)=="taxonomy"] <- paste0("taxonomy.blast.",currName)
print("reading OTU table")
seqTab <- readRDS(snakemake@input[[2]])
if(units=="ASV" & any(colnames(seqTab)=="OTU")){
  seqTab$OTU <- gsub("OTU_","ASV_",seqTab$OTU)
  colnames(seqTab)[which(colnames(seqTab)=="OTU")] <- "ASV" 
}
seqTab <- merge(seqTab,basta,by=units,all.x=T)
seqTab[is.na(seqTab)] <- ""
print("Saving table with all taxonomies")
saveRDS(seqTab,snakemake@output[[1]])
write.table(seqTab,snakemake@output[[2]],,row.names=F,sep="\t",quote=F)

