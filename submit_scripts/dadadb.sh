#! /bin/bash -i

VARCONFIG=$1
DIR=$(dirname $VARCONFIG)
JNAME=$2
WF=$3
EMAIL=$4
REGION=$5
EVALUE=$6
FWD=$7
RVS=$8
FWD_RC=$9
RVS_RC=$10
MISM=$11
TAXNAME=$12
TAX=$13
DB=$14
TAX_OUT=$15
DB_OUT=$16

while IFS=$'\t' read var val; do unset $var ; declare $var="$val" ; done < $VARCONFIG
if [ "$SNAKEMAKE_VIA_CONDA" = true ]; then
   CONDA_START="conda activate $DIR/conda/snakemake_env"
   CONDA_END="conda deactivate"
else
   CONDA_START=""
   CONDA_END=""
fi

eval $LOADING_MODULES
eval $CONDA_START

snakemake $SNAKEMAKE_EXTRA_ARGUMENTS --cores 5 --jobs 5 -s $DIR/DB_snakefile --cluster-config $DIR/config/$SCHEDULER.config.yaml --cluster "{cluster.call} {cluster.runtime}{resources.runtime} {cluster.mem_per_cpu}{resources.mem} {cluster.threads}{threads} {cluster.partition} {cluster.stdout}" --config what=$WF email=$EMAIL sessionName=$JNAME region=$REGION evalue=$EVALUE fwd=$FWD rvs=$RVS fwd_rc=$FWD_RC rvs_rc=$RVS_RC mismatch=$MISM tax=$TAXNAME input_tax=$TAX input_DB=$DB output_tax=$TAX_OUT output_DB=$DB_OUT --use-conda --conda-prefix $DIR/conda $DB_OUT $TAX_OUT >> $JNAME.stdout 2>> $JNAME.stder

eval $CONDA_END
