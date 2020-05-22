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

snakemake -j $THREADS -s $DIR/Snakefile --cluster-config $DIR/config/$SCHEDULER.config.yaml --cluster "{cluster.call} {cluster.runtime}{params.runtime} {cluster.mem_per_cpu} {cluster.threads}{threads} {cluster.partition}" --configfile $CONFIGFILE --config sessionName=$JNAME --use-conda --conda-prefix $DIR/conda >> $JNAME.stdout 2>> $JNAME.stderr

eval $CONDA_END
