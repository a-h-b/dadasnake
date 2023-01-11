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
print(snakemake@input)

print("reading mothur classifier result")
for(mt in 1:length(snakemake@input)){
  currName <- gsub(".tsv$","",gsub(".+Tables/tax.mothur.","",snakemake@input[[mt]]))
  mothTax <- read.delim(snakemake@input[[mt]],
                        stringsAsFactors = F,quote="",header=F)
  colnames(mothTax) <- c(units,"taxonomy.mothur")
  mothTax$taxonomy.mothur <- gsub("^\"","",mothTax$taxonomy.mothur)
  rownames(mothTax) <- mothTax[,which(colnames(mothTax)==units)]
  mothTax$taxonomy.mothur <- gsub("\\([[:digit:]]+\\)","",mothTax$taxonomy.mothur)
  mothTax$taxonomy.mothur[mothTax$taxonomy.mothur==""] <- ";"
  print("formatted")
  tmpTax <- data.frame(do.call(rbind, strsplit(mothTax$taxonomy.mothur, ";", fixed=TRUE)),stringsAsFactors = F)
  rankNum <- as.numeric(snakemake@params[["rank_num"]])
  if(ncol(tmpTax)>= rankNum+1){
   dimCut <- rankNum+1 
   mothTax <- cbind(mothTax,tmpTax[,1:dimCut])
  }else{
   dimCut <- ncol(tmpTax) 
   mothTax <- cbind(mothTax,tmpTax[,1:dimCut],
                    matrix("",nrow=nrow(mothTax),
                           ncol=rankNum+1-dimCut,
                           dimnames=list(mothTax[,which(colnames(mothTax)==units)],
                                     paste0("X",c((dimCut+1):(rankNum+1))))))
  }
  for(i in grep("^X",colnames(mothTax))){
    mothTax[,i][grep("_unclassified$",mothTax[,i])] <- ""
    mothTax[,i][grep("uncultured",mothTax[,i])] <- ""
  }
  prefixVec <- paste0("^",c("k","p","c","o","f","g","s"),"__")
  if(any(grepl(prefixVec[1],mothTax[,3]))){
    for(i in 1:length(prefixVec)){
      mothTax[,2+i][!grepl(prefixVec[i],mothTax[,2+i])] <- ""
      mothTax[,2+i] <- gsub(prefixVec[i],"",mothTax[,2+i])
    }
  }
  if(rankNum %in% c(6,7)){
    colnames(mothTax)[grep("^X",colnames(mothTax))] <- paste0(c("Domain","Phylum","Class","Order","Family","Genus","Species","SH"),
                                                            ".mothur.",currName)
  }else{
    colnames(mothTax)[grep("^X",colnames(mothTax))] <- paste0("Level_",1:(rankNum+1),
                                                            ".mothur.",currName)
   }
   colnames(mothTax)[colnames(mothTax)=="taxonomy.mothur"] <- paste0("taxonomy.mothur.",currName)
  if(!"lastTab" %in% ls()) lastTab <- mothTax else lastTab <- merge(lastTab,mothTax,by=units, all=T)
  lastTab[is.na(lastTab)] <- ""
}

print("Saving table with all taxonomies")
saveRDS(lastTab,snakemake@output[[1]])


