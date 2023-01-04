localrules: post_clustering_control

rule post_clustering_control:
    input:
        "clusteredTables/clusteredTab.tsv",
        "clusteredTables/clusteredTab.RDS"
    output:
        touch("post_clustering.done")


if config['post_clustering']['pre_tax']:
    toClusterTab = "sequenceTables/all.seqTab.RDS"
    if config['post_clustering']['post_ITSx']:
        toCluster = "sequenceTables/ITSx.seqs.fasta"
    else:
        toCluster = "sequenceTables/all.seqs.fasta"
else:
    if config['final_table_filtering']['do'] and 'postprocessing' in STEPS:
        toCluster = "post/filtered.seqs.fasta"
        toClusterTab = "post/filtered.seqTab.RDS"
    else:
        toCluster = "sequenceTables/all.seqs.fasta"
        if 'taxonomy' in STEPS and (config['taxonomy']['decipher']['do'] or config['taxonomy']['mothur']['do'] or config['taxonomy']['dada']['do']):
            if config['blast']['do'] and config['blast']['run_basta']:
                toClusterTab = "sequenceTables/all.seqTab.tax.blast.RDS"
            else: 
                toClusterTab = "sequenceTables/all.seqTab.tax.tsv"
        else:
            toClusterTab = "sequenceTables/all.seqTab.RDS"

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
            vsearch --cluster_fast {input} \
             --centroids {output.cent} --id {params.id} --consout {output.cons} \
             --strand {params.strand} --uc {output.uc} --clusterout_id &> {log}
            """

    rule vsearch_table:
        input:
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


