postConts = []
if config['hand_off']['phyloseq']:
    postConts.append("sequenceTables/all.seqTab.phyloseq.RDS")
    if 'postclustering' in STEPS:
        postConts.append("clusteredTables/clusteredTab.phyloseq.RDS")
if config['postprocessing']['treeing']['do']:
    postConts.append("post/ASV.tree.newick")
    if 'postclustering' in STEPS:
        postConts.append("post/clustered.tree.newick")
if config['postprocessing']['rarefaction_curve']:
    if 'postclustering' in STEPS:
        postConts.append("stats/rarefaction_curves.clustered.pdf")
    postConts.append("stats/rarefaction_curves.ASV.pdf")
if config['postprocessing']['funguild']['do'] and 'taxonomy' in STEPS:
    CLASSIFY=config['postprocessing']['funguild']['classifier_db'].split(".")[0]
    if config['taxonomy'][CLASSIFY]['do'] and 'ASV' in config['taxonomy'][CLASSIFY]['run_on']:
        postConts.append("post/all.seqTab.guilds.tsv")
    if config['taxonomy'][CLASSIFY]['do'] and 'cluster' in config['taxonomy'][CLASSIFY]['run_on']:
        postConts.append("post/clusteredTab.guilds.tsv")
if config['postprocessing']['fungalTraits']['do'] and 'taxonomy' in STEPS:
    CLASSIFY=config['postprocessing']['fungalTraits']['classifier_db'].split(".")[0]
    if config['taxonomy'][CLASSIFY]['do'] and 'ASV' in config['taxonomy'][CLASSIFY]['run_on']:
        postConts.append("post/all.seqTab.traits.RDS")
    if config['taxonomy'][CLASSIFY]['do'] and 'cluster' in config['taxonomy'][CLASSIFY]['run_on']:
        postConts.append("post/clusteredTab.traits.RDS")
if config['postprocessing']['tax4fun2']['do']:
    postConts.append("post/tax4fun2/KOs_per_ASV.txt")

localrules: post_control_noFilter

rule post_control_noFilter:
    input:
        postConts
    output:
        "postprocessing.done"
    shell:
        """
        touch {output}
        """

rule rarefaction_curve_noFilter:
    input:
        lambda wildcards: "sequenceTables/all.seqTab.RDS" if wildcards.kind=="ASV" else "clusteredTables/clusteredTab.RDS",
        "reporting/finalNumbers_perSample.tsv"
    output:
        "stats/rarefaction_curves.{kind}.pdf"
    params:
        thing=lambda wildcards: "sequence variants" if wildcards.kind=="ASV" else "OTU clusters"
    threads: 1
    resources:
        runtime="48:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/rarefaction_curve.{kind}.log"
    script:
        SCRIPTSDIR+"rarefaction_curve.R"

rule guilds_noFilter:
    input:
        lambda wildcards: "sequenceTables/" + wildcards.prefix + ".tax.tsv" if wildcards.prefix=="all.seqTab" else "clusteredTables/" + wildcards.prefix + ".tax.tsv"
    output:
        "post/{prefix}.guilds.tsv"
    threads: 1
    resources:
        runtime="24:00:00",
        mem=config['normalMem']
    params:
        src_path=SCRIPTSDIR,
        out_path="post",
        funguild_db=config['postprocessing']['funguild']['funguild_db'],
        classifier=config['postprocessing']['funguild']['classifier_db']
    log: "logs/funguild.{prefix}.log"
    conda: ENVDIR + "dadasnake_env.yml"
    message: "Running funguild on {input}."
    shell:
        """
        mkdir -p {params.out_path}
        {params.src_path}/Guilds_v1.1.local.2.py -otu {input} -output {output} -path_to_db {params.funguild_db} -taxonomy_name taxonomy.{params.classifier}&> {log} || touch {output}
        """


rule funTraits_noFilter:
    input:
        lambda wildcards: "sequenceTables/" + wildcards.prefix + ".tax.RDS" if wildcards.prefix=="all.seqTab" else "clusteredTables/" + wildcards.prefix + ".tax.RDS"
    output:
        "post/{prefix}.traits.tsv",
        "post/{prefix}.traits.RDS"
    params:
        db=config['postprocessing']['fungalTraits']['db'],
        level=config['postprocessing']['fungalTraits']['level'],
        classifier=config['postprocessing']['fungalTraits']['classifier_db']
    threads: 1
    resources:
        runtime="6:00:00",
        mem=config['normalMem']
    log: "logs/fungalTraits.{prefix}.log"
    conda: ENVDIR + "dada2_env.yml"
    message: "Adding fungalTraits to {input}."
    script:
        SCRIPTSDIR + "add_fungalTraits.R"

rule tax4fun2_noFilter:
    input:
        "sequenceTables/all.seqs.fasta",
        "sequenceTables/all.seqTab.RDS"
    output:
        "post/tax4fun2/KOs_per_ASV.txt",
        "post/tax4fun2/pathway_per_ASV.txt",
        "post/tax4fun2/functional_prediction.txt",
        "post/tax4fun2/pathway_prediction.txt"
    threads: getThreads(4) 
    resources:
        runtime="6:00:00",
        mem=config['normalMem']
    params:
        tmp=TMPDIR,
        outputDir="post/tax4fun2",
        customFunc=SCRIPTSDIR + "/functionalPredictionCustom.R",
        db=config["postprocessing"]["tax4fun2"]["db"],
        user_data=config["postprocessing"]["tax4fun2"]["user_data"],
        user_dir=config["postprocessing"]["tax4fun2"]["user_dir"],
        user_db=config["postprocessing"]["tax4fun2"]["user_db"],
        database_mode=config["postprocessing"]["tax4fun2"]["database_mode"],
        normalize_by_copy_number=config["postprocessing"]["tax4fun2"]["normalize_by_copy_number"],
        min_identity_to_reference=config["postprocessing"]["tax4fun2"]["min_identity_to_reference"]
    log: "logs/tax4fun2.log"
    conda: ENVDIR + "tax4fun2_env.yml"
    message: "Running tax4fun2 on {input}."
    script:
        SCRIPTSDIR + "tax4fun2.R"


if config['hand_off']['phyloseq']:
    if 'taxonomy' not in STEPS:
        physInputs = ["sequenceTables/all.seqTab.RDS","reporting/finalNumbers_perSample.tsv"]
    elif config['taxonomy']['dada']['do'] and 'ASV' in config['taxonomy']['dada']['run_on']:
        physInputs = ["sequenceTables/all.seqTab.tax.RDS","reporting/finalNumbers_perSample.tsv"]
    elif config['taxonomy']['decipher']['do'] and 'ASV' in config['taxonomy']['decipher']['run_on']:
        physInputs = ["sequenceTables/all.seqTab.tax.RDS","reporting/finalNumbers_perSample.tsv"]
    elif config['taxonomy']['mothur']['do'] and 'ASV' in config['taxonomy']['mothur']['run_on']:
        physInputs = ["sequenceTables/all.seqTab.tax.RDS","reporting/finalNumbers_perSample.tsv"]
    else:
        physInputs = ["sequenceTables/all.seqTab.RDS","reporting/finalNumbers_perSample.tsv"]
    if 'taxonomy' in STEPS and config['blast']['do'] and config['blast']['run_basta'] and 'ASV' in config['blast']['run_on']:
        physInputs = ["sequenceTables/all.seqTab.tax.blast.RDS","reporting/finalNumbers_perSample.tsv"]
    if config['postprocessing']['treeing']['do']:
        physInputs.append("post/ASV.tree.newick")
    if 'postclustering' in STEPS:
        physInputs.append("clusteredTables/cluster_info.tsv")
    rule phyloseq_handoff_post:
        input:
            physInputs
        output:
            "sequenceTables/all.seqTab.phyloseq.RDS"
        threads: 1
        params:
            currentStep = "post",
            what = "ASV",
            tree = "add" if config['postprocessing']['treeing']['do'] else "none",
            cluster = "add" if 'postclustering' in STEPS else "none",
            clusterMeth = config['post_clustering']['method']
        resources:
            runtime="4:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "add_R_env.yml"
        log: "logs/phyloseq_hand-off.log"
        script:
            SCRIPTSDIR+"phyloseq_handoff.R"

    if 'taxonomy' not in STEPS:
        physInputsCl = ["clusteredTables/clusteredTab.RDS","reporting/finalNumbers_perSample.tsv"]
    elif config['taxonomy']['dada']['do'] and 'cluster' in config['taxonomy']['dada']['run_on']:
        physInputsCl = ["clusteredTables/clusteredTab.tax.RDS","reporting/finalNumbers_perSample.tsv"]
    elif config['taxonomy']['decipher']['do'] and 'cluster' in config['taxonomy']['decipher']['run_on']:
        physInputsCl = ["clusteredTables/clusteredTab.tax.RDS","reporting/finalNumbers_perSample.tsv"]
    elif config['taxonomy']['mothur']['do'] and 'cluster' in config['taxonomy']['mothur']['run_on']:
        physInputsCl = ["clusteredTables/clusteredTab.tax.RDS","reporting/finalNumbers_perSample.tsv"]
    else:
        physInputsCl = ["clusteredTables/clusteredTab.RDS","reporting/finalNumbers_perSample.tsv"]
    if 'taxonomy' in STEPS and config['blast']['do'] and config['blast']['run_basta'] and 'cluster' in config['blast']['run_on']:
        physInputsCl = ["clusteredTables/clusteredTab.tax.blast.RDS","reporting/finalNumbers_perSample.tsv"]
    if config['postprocessing']['treeing']['do']:
        physInputsCl.append("post/cluster.tree.newick")
    rule phyloseq_handoff_postCl:
        input:
            physInputsCl
        output:
            "clusteredTables/clusteredTab.phyloseq.RDS"
        threads: 1
        params:
            currentStep = "post",
            what = "cluster",
            tree = "add" if config['postprocessing']['treeing']['do'] else "none",
            cluster = "none"
        resources:
            runtime="4:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "add_R_env.yml"
        log: "logs/phyloseq_hand-off.cluster.log"
        script:
            SCRIPTSDIR+"phyloseq_handoff.R"


rule multiAlign_noFilter:
    input:
        lambda wildcards: "sequenceTables/" + wildcards.prefix + ".fasta" if wildcards.prefix=="all.seqs" else "clusteredTables/" + wildcards.prefix + ".fasta"
    output:
        "post/{prefix}.multi.fasta"
    threads: getThreads(10)
    resources:
        runtime="24:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dadasnake_env.yml"
    log: "logs/treeing_multiAlign.{prefix}.log"
    shell:
        """
        clustalo -i {input} -o {output} --outfmt=fasta --threads={threads} --force &> {log} || touch {output}
        """

if config['postprocessing']['treeing']['fasttreeMP'] != "":
    rule treeing_noFilter_fasttreeMP:
        input:
            lambda wildcards: "post/all.seqs.multi.fasta" if wildcards.kind=="ASV" else "post/consensus.multi.fasta"
        output:
            "post/{kind}.tree.newick"
        params:
            ftMP=config['postprocessing']['treeing']['fasttreeMP']
        threads: getThreads(10)
        resources:
            runtime="48:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/treeing.{kind}.log"
        shell:
            """
            {params.ftMP} -nt -gamma -no2nd -fastest -spr 4 \
             -log {log} -quiet {input} > {output} 2> {log}  || touch {output}
            """
else:
    rule treeing_noFilter:
        input:
            lambda wildcards: "post/all.seqs.multi.fasta" if wildcards.kind=="ASV" else "post/consensus.multi.fasta"
        output:
            "post/{kind}.tree.newick"
        threads: 1
        resources:
            runtime="48:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/treeing.{kind}.log"
        shell:
            """
            fasttree -nt -gamma -no2nd -fastest -spr 4 \
             -log {log} -quiet {input} > {output} 2>> {log} || touch {output}
            """




