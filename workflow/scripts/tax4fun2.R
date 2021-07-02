log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

if(!require(Tax4Fun2)){
  system(paste0("wget -O ",snakemake@params[['tmp']],"/Tax4Fun2_1.1.5.tar.gz https://github.com/bwemheu/Tax4Fun2/releases/download/1.1.5/Tax4Fun2_1.1.5.tar.gz"))
  install.packages(pkgs = paste0(snakemake@params[['tmp']],"/Tax4Fun2_1.1.5.tar.gz"), repos = NULL, source = TRUE)
  require(Tax4Fun2)
}

source(snakemake@params[["customFunc"]])
if(!file.exists(paste0(snakemake@config[["postprocessing"]][["tax4fun2"]][["db"]],"/blast_bin/bin/makeblastdb"))){
  system(paste0("mkdir -p ",snakemake@config[["postprocessing"]][["tax4fun2"]][["db"]],"/blast_bin/bin/"))
  system(paste0("ln -s ",condap,"/bin/blastn ",snakemake@config[["postprocessing"]][["tax4fun2"]][["db"]],"/blast_bin/bin/blastn"))
  system(paste0("ln -s ",condap,"/bin/makeblastdb ",snakemake@config[["postprocessing"]][["tax4fun2"]][["db"]],"/blast_bin/bin/makeblastdb"))
}

threads <- snakemake@threads
include_user_data <- snakemake@config[["postprocessing"]][["tax4fun2"]][["user_data"]]
if(include_user_data){
 path_to_user_data <- snakemake@config[["postprocessing"]][["tax4fun2"]][["user_dir"]]
 name_of_user_data <- snakemake@config[["postprocessing"]][["tax4fun2"]][["user_db"]]
}else{
 path_to_user_data <- ""
 name_of_user_data <- ""
}


# run blast for tax4fun2
print("running blast for tax4fun2")
runRefBlast(path_to_otus = snakemake@input[[1]], 
            path_to_reference_data = snakemake@config[["postprocessing"]][["tax4fun2"]][["db"]], 
            path_to_temp_folder = snakemake@params[["outputDir"]], 
            database_mode = snakemake@config[["postprocessing"]][["tax4fun2"]][["database_mode"]], 
            use_force = T, num_threads = threads,
            include_user_data = include_user_data, 
            path_to_user_data = path_to_user_data,
            name_of_user_data = name_of_user_data)

# functional predictions
print("running functional predictions")
makeFunctionalPredictionCustom(path_to_otu_table = snakemake@input[[2]],
                               path_to_reference_data = snakemake@config[["postprocessing"]][["tax4fun2"]][["db"]],
                         path_to_temp_folder = snakemake@params[["outputDir"]],
                         database_mode = snakemake@config[["postprocessing"]][["tax4fun2"]][["database_mode"]], 
                         normalize_by_copy_number = snakemake@config[["postprocessing"]][["tax4fun2"]][["normalize_by_copy_number"]], 
                         min_identity_to_reference = as.numeric(snakemake@config[["postprocessing"]][["tax4fun2"]][["min_identity_to_reference"]]), 
                         normalize_pathways = FALSE,
                         include_user_data = include_user_data,
                         path_to_user_data = path_to_user_data,
                         name_of_user_data = name_of_user_data)
print("done")
