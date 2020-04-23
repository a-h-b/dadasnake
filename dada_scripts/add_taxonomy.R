log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(Biostrings)


print("reading OTU table")
seqTab <- readRDS(snakemake@input[[1]])
if(snakemake@config[['ITSx']][['do']]){
  print("reading ITSx result")
  iadd <- 1
  seqs <- readDNAStringSet(snakemake@input[[2]])
  seqTab$ITSx <- seqTab$OTU %in% names(seqs)
  rm(seqs)
}else{
  iadd <- 0
}
if(snakemake@config[['taxonomy']][['decipher']][['do']]){
  print("reading DECIPHER result")
  decTax <- readRDS(snakemake@input[[2+iadd]])
  decTax <- decTax[,-ncol(decTax)]
  colnames(decTax) <- paste0(colnames(decTax),".decipher")
  seqTab <- merge(seqTab,decTax,by.x="OTU",by.y=0,all=T)
}
if(snakemake@config[['taxonomy']][['mothur']][['do']]){
  print("reading mothur classifier result")
  mothTax <- read.delim(snakemake@input[[length(snakemake@input)]],
                        stringsAsFactors = F,quote="",header=F)
  colnames(mothTax) <- c("OTU","taxonomy.mothur")
  mothTax$taxonomy.mothur <- gsub("^\"","",mothTax$taxonomy.mothur)
  rownames(mothTax) <- mothTax$OTU
  mothTax$taxonomy.mothur <- gsub("\\([[:digit:]]+\\)","",mothTax$taxonomy.mothur)
  mothTax$taxonomy.mothur[mothTax$taxonomy.mothur==""] <- ";"
  tmpTax <- data.frame(do.call(rbind, strsplit(mothTax$taxonomy.mothur, ";", fixed=TRUE)),stringsAsFactors = F)
  if(ncol(tmpTax)>=7){
   dimCut <- 7 
   mothTax <- cbind(mothTax,tmpTax[,1:dimCut])
  }else{
   dimCut <- ncol(tmpTax) 
   mothTax <- cbind(mothTax,tmpTax[,1:dimCut],
                    matrix("",nrow=nrow(mothTax),ncol=7-dimCut,
                           dimnames=list(mothTax$OTU,paste0("X",c((dimCut+1):7)))))
  } 
  prefixVec <- paste0("^",c("k","p","c","o","f","g","s"),"__")
  for(i in grep("^X",colnames(mothTax))){
    mothTax[,i][grep("_unclassified$",mothTax[,i])] <- ""
    mothTax[,i][grep("uncultured",mothTax[,i])] <- ""
  }
  if(any(grepl(prefixVec[1],mothTax[,3]))){
    for(i in 1:length(prefixVec)){
      mothTax[,2+i][!grepl(prefixVec[i],mothTax[,2+i])] <- ""
      mothTax[,2+i] <- gsub(prefixVec[i],"",mothTax[,2+i])
    }
  }
  colnames(mothTax)[grep("^X",colnames(mothTax))] <- paste0(c("Domain","Phylum","Class","Order","Family","Genus","Species"),
                                                            ".mothur")
  seqTab <- merge(seqTab,mothTax,by="OTU",all=T)
}
print("Saving OTU table with taxonomy")
saveRDS(seqTab,snakemake@output[[2]])
write.table(seqTab,snakemake@output[[1]],row.names=F,sep="\t",quote=F)


