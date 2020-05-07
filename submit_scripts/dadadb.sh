#! /bin/bash -i

CONFIGFILE=$1
VARCONFIG=$2
DIR=$(dirname $VARCONFIG)
JNAME=$3

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

snakemake -j 5 -s $DIR/DB_snakefile --cluster-config $DIR/config/$SCHEDULER.config.yaml --cluster "{cluster.call} {cluster.runtime}{params.runtime} {cluster.mem_per_cpu} {cluster.threads}{threads} {cluster.partition}" --configfile $CONFIGFILE --config what=$WF email=$EMAIL sessionName=$JNAME region=$REGION evalue=$EVALUE fwd=$FWD rvs=$RVS fwd_rc=$FWD_RC rvs_rc=$RVS_RC mismatch=$MISM tax=$TAXNAME input_tax=$TAX input_DB=$DB output_tax=$TAX_OUT output_DB=$DB_OUT --use-conda --conda-prefix $DIR/conda $DB_OUT $TAX_OUT >> $JNAME.stdout 2>> $JNAME.stder

eval $CONDA_END
