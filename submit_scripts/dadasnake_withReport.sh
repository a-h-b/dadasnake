#! /bin/bash -i

CONFIGFILE=$1
VARCONFIG=$2
DIR=$(dirname $VARCONFIG)
JNAME=$3
THREADS=$4

while read var val; do unset $var ; declare $var="$val" ; done < $VARCONFIG
if [ "$SNAKEMAKE_VIA_CONDA" = true ]; then
   CONDA_START="conda activate $DIR/conda/snakemake_env"
   CONDA_END="conda deactivate"
else
   CONDA_START=""
   CONDA_END=""
fi

eval $LOADING_MODULES
eval $CONDA_START

snakemake --cores $THREADS -s $DIR/Snakefile --keep-going --local-cores 1 --cluster-config $DIR/config/$SCHEDULER.config.yaml --cluster "{cluster.call} {cluster.runtime}{resources.runtime} {cluster.mem_per_cpu}{resources.mem} {cluster.threads}{threads} {cluster.partition}" --configfile $CONFIGFILE --config sessionName=$JNAME normalMem=$NORMAL_MEM_EACH bigMem=$BIGMEM_MEM_EACH bigCores=$BIGMEM_CORES settingsLocked=$LOCK_SETTINGS sessionKind=cluster --use-conda --conda-prefix $DIR/conda >> $JNAME.stdout 2>> $JNAME.stderr

snakemake --cores 1 -s $DIR/Snakefile --report report.html --configfile $CONFIGFILE

eval $CONDA_END
