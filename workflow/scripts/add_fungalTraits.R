log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

print("reading OTU table")
seqTab <- readRDS(snakemake@input[[1]])
ori_colnames <- colnames(seqTab)

print("reading FungalTraits table")
fun <- read.delim(snakemake@params[['db']])
colnames(fun) <- paste0(colnames(fun),".FungalTraits")
print(snakemake@params[['level']])
if(is.character(snakemake@params[['level']])){ 
  levelOI <- snakemake@params[['level']]
}else{
  levelOI <- paste0("Level_",snakemake@params[['level']])
} 
print(levelOI)
colOI <- grep(paste0("^",levelOI,".",
                            snakemake@params[['classifier']]),
              colnames(seqTab),value=T)

seqTab <- merge(seqTab,fun,
                by.x=colOI,
                by.y="GENUS.FungalTraits",all.x=T,sort=F)
seqTab <- seqTab[,c(ori_colnames,setdiff(colnames(seqTab),ori_colnames))]


print("Saving OTU table with FungalTraits")
saveRDS(seqTab,snakemake@output[[2]])
write.table(seqTab,snakemake@output[[1]],row.names=F,sep="\t",quote=F)


