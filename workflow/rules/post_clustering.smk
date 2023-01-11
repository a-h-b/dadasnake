localrules: post_clustering_control

rule post_clustering_control:
    input:
        "clusteredTables/clusteredTab.tsv",
        "clusteredTables/clusteredTab.RDS"
    output:
        touch("post_clustering.done")


toCluster = "sequenceTables/all.seqs.fasta"
if 'taxonomy' not in STEPS:
    toClusterTab = "sequenceTables/all.seqTab.RDS"
else:
    if config['taxonomy']['decipher']['do'] and 'ASV' in config['taxonomy']['decipher']['run_on']:
        toClusterTab = "sequenceTables/all.seqTab.tax.RDS"
    elif config['taxonomy']['mothur']['do'] and 'ASV' in config['taxonomy']['mothur']['run_on']:
        toClusterTab = "sequenceTables/all.seqTab.tax.RDS"
    elif config['taxonomy']['dada']['do'] and 'ASV' in config['taxonomy']['dada']['run_on']:
        toClusterTab = "sequenceTables/all.seqTab.tax.RDS"
    else:
        toClusterTab = "sequenceTables/all.seqTab.RDS"
    if config['blast']['do'] and 'ASV' in config['blast']['run_on'] and config['blast']['run_basta']:
        toClusterTab = "sequenceTables/all.seqTab.tax.blast.RDS"



if config['post_clustering']['method'] == "vsearch":    
    rule vsearch_cluster:
        input:
            toCluster
        output:
            cent="clusteredTables/centroids.fasta",
            cons="clusteredTables/consensus.fasta",
            uc="clusteredTables/cluster_info.tsv"
        params:
            id=config['post_clustering']['cutoff'] if config['post_clustering']['cutoff'] < 1 else config['post_clustering']['cutoff']/100.0,
            strand=config['post_clustering']['strand'] 
        threads: getThreads(12)
        resources:
            runtime="2:00:00",
            mem=config['normalMem']
        conda: os.path.join(ENVDIR, "vsearch_env.yml")
        log: "logs/post_clustering_vsearch.log"
        message: "Running vsearch clustering on {input}."
        shell:
            """
            TMPD=$(mktemp -d -t --tmpdir={TMPDIR} 'XXXXXX') 
            sed "s#^>OTU_#>ASV_#g" {input} > $TMPD/asvs.fasta 
            vsearch --cluster_fast $TMPD/asvs.fasta \
             --centroids {output.cent} --id {params.id} --consout {output.cons} \
             --strand {params.strand} --uc {output.uc} --clusterout_id &> {log}
            awk -i inplace -v OFS="" -F ";|=|_" '{{if(/^>/){{$1=sprintf(">OTU_%0" length($2) "d",$4);print$1}}else{{print $1}}}}' {output.cent}
            awk -i inplace -v OFS="" -F ";|=|_" '{{if(/^>/){{$1=sprintf(">OTU_%0" length($3) "d",$7);print$1}}else{{print $1}}}}' {output.cons}
            """

    rule vsearch_table:
        input:
            "clusteredTables/consensus.fasta",
            toClusterTab,
            "clusteredTables/cluster_info.tsv"
        output:
            "clusteredTables/clusteredTab.tsv",
            "clusteredTables/clusteredTab.RDS"
        threads: 1
        resources:
            runtime="2:00:00",
            mem=config['normalMem']
        conda: os.path.join(ENVDIR, "dada2_env.yml")
        log: "logs/post_clustering_vsearch.log"
        message: "Running vsearch clustering on {input}."
        script:
            os.path.join(SCRIPTSDIR,"post_clustering_gatherVsearch.R")

else:
    rule decipher_cluster:
        input:
            toCluster,
            toClusterTab
        output:
            "clusteredTables/cluster_info.tsv",
            "clusteredTables/consensus.fasta",
            "clusteredTables/clusteredTab.tsv",
            "clusteredTables/clusteredTab.RDS"
        params:
            cutoff=config['post_clustering']['cutoff']
        threads: getThreads(12)
        resources:
            runtime="2:00:00",
            mem=config['normalMem']
        conda: os.path.join(ENVDIR, "dada2_env.yml")
        log: "logs/post_clustering_decipher.log"
        message: "Running decipher clustering on {input}."
        script:
            os.path.join(SCRIPTSDIR,"post_clustering_decipher.R")


