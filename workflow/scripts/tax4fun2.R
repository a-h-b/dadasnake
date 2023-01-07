log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

if(!require(Tax4Fun2)){
  system(paste0("wget -O ",snakemake@params[['tmp']],"/Tax4Fun2_1.1.5.tar.gz https://github.com/ZihaoShu/Tax4Fun2/raw/main/Tax4Fun2_1.1.5.tar.gz"))
  install.packages(pkgs = paste0(snakemake@params[['tmp']],"/Tax4Fun2_1.1.5.tar.gz"), repos = NULL, source = TRUE)
  require(Tax4Fun2)
}

source(snakemake@params[["customFunc"]])
if(!file.exists(paste0(snakemake@params[["db"]],"/blast_bin/bin/makeblastdb"))){
  system(paste0("mkdir -p ",snakemake@params[["db"]],"/blast_bin/bin/"))
  system(paste0("ln -s ",condap,"/bin/blastn ",snakemake@params[["db"]],"/blast_bin/bin/blastn"))
  system(paste0("ln -s ",condap,"/bin/makeblastdb ",snakemake@params[["db"]],"/blast_bin/bin/makeblastdb"))
}
if(!file.exists(paste0(snakemake@params[["db"]],"/Ref99NR"){
  system(paste0("wget -O ",snakemake@params[['db']],"/Ref99NR.zip https://cloudstor.aarnet.edu.au/plus/s/DkoZIyZpMNbrzSw/download"))
  system(paste0("unzip ",snakemake@params[['db']],"/Ref99NR.zip"))
}
if(!file.exists(paste0(snakemake@params[["db"]],"/Ref100NR"){
  system(paste0("wget -O ",snakemake@params[['db']],"/Ref100NR.zip https://cloudstor.aarnet.edu.au/plus/s/jIByczak9ZAFUB4/download"))
  system(paste0("unzip ",snakemake@params[['db']],"/Ref100NR.zip"))
}

threads <- snakemake@threads
include_user_data <- snakemake@params[["user_data"]]
if(include_user_data){
 path_to_user_data <- snakemake@params[["user_dir"]]
 name_of_user_data <- snakemake@params[["user_db"]]
}else{
 path_to_user_data <- ""
 name_of_user_data <- ""
}


# run blast for tax4fun2
print("running blast for tax4fun2")
runRefBlast(path_to_otus = snakemake@input[[1]], 
            path_to_reference_data = snakemake@params[["db"]], 
            path_to_temp_folder = snakemake@params[["outputDir"]], 
            database_mode = snakemake@params[["database_mode"]], 
            use_force = T, num_threads = threads,
            include_user_data = include_user_data, 
            path_to_user_data = path_to_user_data,
            name_of_user_data = name_of_user_data)

# functional predictions
print("running functional predictions")
makeFunctionalPredictionCustom(path_to_otu_table = snakemake@input[[2]],
                               path_to_reference_data = snakemake@params[["db"]],
                         path_to_temp_folder = snakemake@params[["outputDir"]],
                         database_mode = snakemake@params[["database_mode"]], 
                         normalize_by_copy_number = snakemake@params[["normalize_by_copy_number"]], 
                         min_identity_to_reference = as.numeric(snakemake@params[["min_identity_to_reference"]]), 
                         normalize_pathways = FALSE,
                         include_user_data = include_user_data,
                         path_to_user_data = path_to_user_data,
                         name_of_user_data = name_of_user_data)
print("done")
