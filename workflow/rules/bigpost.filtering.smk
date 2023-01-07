postConts = ["post/filtered.seqTab.RDS","reporting/post_finalNumbers_perSample.tsv"]
if 'postclustering' in STEPS:
    postConts.append("post/filtered.clusteredTab.RDS")
if config['hand_off']['phyloseq']:
    postConts.append("post/filtered.seqTab.phyloseq.RDS")
    if 'postclustering' in STEPS:
        postConts.append("post/filtered.clusteredTab.phyloseq.RDS")
if config['postprocessing']['treeing']['do']:
    postConts.append("post/ASV.tree.newick")
    if 'postclustering' in STEPS:
        postConts.append("post/clustered.tree.newick")
if config['postprocessing']['rarefaction_curve']:
    postConts.append("stats/rarefaction_curves.ASV.pdf")
    if 'postclustering' in STEPS:
        postConts.append("stats/rarefaction_curves.clustered.pdf")
if config['postprocessing']['funguild']['do'] and 'taxonomy' in STEPS:
    CLASSIFY=config['postprocessing']['funguild']['classifier_db'].split(".")[0]
    if config['taxonomy'][CLASSIFY]['do'] and 'ASV' in config['taxonomy'][CLASSIFY]['run_on']:
        postConts.append("post/filtered.seqTab.guilds.tsv")
    if config['taxonomy'][CLASSIFY]['do'] and 'cluster' in config['taxonomy'][CLASSIFY]['run_on']:
        postConts.append("post/filtered.clusteredTab.guilds.tsv")
if config['postprocessing']['fungalTraits']['do'] and 'taxonomy' in STEPS:
    CLASSIFY=config['postprocessing']['fungalTraits']['classifier_db'].split(".")[0]
    if config['taxonomy'][CLASSIFY]['do'] and 'ASV' in config['taxonomy'][CLASSIFY]['run_on']:
        postConts.append("post/filtered.seqTab.traits.RDS")
    if config['taxonomy'][CLASSIFY]['do'] and 'cluster' in config['taxonomy'][CLASSIFY]['run_on']:
        postConts.append("post/filtered.clusteredTab.traits.RDS")
if config['postprocessing']['tax4fun2']['do']:
    postConts.append("post/tax4fun2/KOs_per_ASV.txt")

localrules: post_control_Filter

rule post_control_Filter:
    input:
        postConts
    output:
        "postprocessing.done"
    shell:
        """
        touch {output}
        """

filtTabs = ["sequenceTables/all.seqs.fasta"]
if 'taxonomy' not in STEPS:
    filtTabs = ["sequenceTables/all.seqs.fasta","sequenceTables/all.seqTab.RDS"]
elif config['taxonomy']['dada']['do'] and 'ASV' in config['taxonomy']['dada']['run_on']:
    filtTabs = ["sequenceTables/all.seqs.fasta","sequenceTables/all.seqTab.tax.RDS"]
elif config['taxonomy']['decipher']['do'] and 'ASV' in config['taxonomy']['decipher']['run_on']:
    filtTabs = ["sequenceTables/all.seqs.fasta","sequenceTables/all.seqTab.tax.RDS"]
elif config['taxonomy']['mothur']['do'] and 'ASV' in config['taxonomy']['mothur']['run_on']:
    filtTabs = ["sequenceTables/all.seqs.fasta","sequenceTables/all.seqTab.tax.RDS"]
else:
    filtTabs = ["sequenceTables/all.seqs.fasta","sequenceTables/all.seqTab.RDS"]
if 'taxonomy' in STEPS and config['blast']['do'] and config['blast']['run_basta'] and 'ASV' in config['blast']['run_on']:
    filtTabs = ["sequenceTables/all.seqs.fasta","sequenceTables/all.seqTab.tax.blast.RDS"]
rule filtering_table:
    input: 
       filtTabs
    output:
       "post/filtered.seqTab.RDS",
       "post/filtered.seqTab.tsv",   
       "post/filtered.seqs.fasta"
    params:
        what="ASV",
        target_taxa=config['final_table_filtering']['keep_target_taxa'],
        min=config['final_table_filtering']['target_min_length'],
        max=config['final_table_filtering']['target_max_length'] 
    threads: config['bigCores']
    resources:
        runtime="8:00:00",
        mem=config['bigMem']
    log: "logs/post_filtering_table.log"
    conda: ENVDIR + "dada2_env.yml"
    script:
        SCRIPTSDIR+"post_filtering.R"


if 'taxonomy' not in STEPS:
    filtTabsCl = ["clusteredTables/consensus.fasta","clusteredTables/clusteredTab.RDS"]
elif config['taxonomy']['dada']['do'] and 'cluster' in config['taxonomy']['dada']['run_on']:
    filtTabsCl = ["clusteredTables/consensus.fasta","clusteredTables/clusteredTab.tax.RDS"]
elif config['taxonomy']['decipher']['do'] and 'cluster' in config['taxonomy']['decipher']['run_on']:
    filtTabsCl = ["clusteredTables/consensus.fasta","clusteredTables/clusteredTab.tax.RDS"]
elif config['taxonomy']['mothur']['do'] and 'cluster' in config['taxonomy']['mothur']['run_on']:
    filtTabsCl = ["clusteredTables/consensus.fasta","clusteredTables/clusteredTab.tax.RDS"]
else:
    filtTabsCl = ["clusteredTables/consensus.fasta","clusteredTables/clusteredTab.RDS"]
if 'taxonomy' in STEPS and config['blast']['do'] and config['blast']['run_basta'] and 'cluster' in config['blast']['run_on']:
    filtTabsCl = ["clusteredTables/consensus.fasta","clusteredTables/clusteredTab.tax.blast.RDS"]
rule filtering_tableCl:
    input:
       filtTabsCl
    output:
       "post/filtered.clusteredTab.RDS",
       "post/filtered.clusteredTab.tsv",
       "post/filtered.consensus.fasta"
    params:
        what="cluster",
        target_taxa=config['final_table_filtering']['keep_target_taxa'],
        min=config['final_table_filtering']['target_min_length'],
        max=config['final_table_filtering']['target_max_length']
    threads: config['bigCores']
    resources:
        runtime="8:00:00",
        mem=config['bigMem']
    log: "logs/post_filtering_table.log"
    conda: ENVDIR + "dada2_env.yml"
    script:
        SCRIPTSDIR+"post_filtering.R"

rule table_filter_numbers:
    input:
        "reporting/finalNumbers_perSample.tsv",
        "post/filtered.seqTab.RDS"
    output:
        "reporting/post_finalNumbers_perSample.tsv"
    threads: config['bigCores']
    params:
        currentStep = "post"
    resources:
        runtime="12:00:00",
        mem=config['bigMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/countPostfilteredReads.log"
    script:
        SCRIPTSDIR+"report_readNumbers.R"


rule rarefaction_curve_Filter:
    input:
        lambda wildcards: "post/filtered.seqTab.RDS" if wildcards.kind=="ASV" else "post/filtered.clusteredTab.RDS",
        "reporting/finalNumbers_perSample.tsv"
    output:
        "stats/rarefaction_curves.{kind}.pdf"
    params:
        thing=lambda wildcards: "sequence variants" if wildcards.kind=="ASV" else "OTU clusters"
    threads: config['bigCores']
    resources:
        runtime="48:00:00",
        mem=config['bigMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/rarefaction_curve.{kind}.log"
    script:
        SCRIPTSDIR+"rarefaction_curve.R"


rule guilds_Filter:
    input:
        "post/filtered.{tab}.tsv"
    output:
        "post/filtered.{tab}.guilds.tsv"
    threads: 1
    resources:
        runtime="24:00:00",
        mem=config['bigMem']
    params:
        src_path=SCRIPTSDIR,
        out_path="post",
        funguild_db=config['postprocessing']['funguild']['funguild_db'],
        classifier=config['postprocessing']['funguild']['classifier_db']        
    log: "logs/funguild.{tab}.log"
    conda: ENVDIR + "dadasnake_env.yml"
    message: "Running funguild on {input}."
    shell:
        """
        mkdir -p {params.out_path}
        {params.src_path}/Guilds_v1.1.local.2.py -otu {input} -output {output} -path_to_db {params.funguild_db} -taxonomy_name taxonomy.{params.classifier} &> {log} || touch {output}
        """

rule funTraits_Filter:
    input:
        "post/filtered.{tab}.RDS"
    output:
        "post/filtered.{tab}.traits.tsv",
        "post/filtered.{tab}.traits.RDS"
    params:
        db=config['postprocessing']['fungalTraits']['db'],
        level=config['postprocessing']['fungalTraits']['level'],
        classifier=config['postprocessing']['fungalTraits']['classifier_db']
    threads: config['bigCores']
    resources:
        runtime="12:00:00",
        mem=config['bigMem']
    log: "logs/fungalTraits.{tab}.log"
    conda: ENVDIR + "dada2_env.yml"
    message: "Adding fungalTraits to {input}."
    script:
        SCRIPTSDIR + "add_fungalTraits.R"

rule tax4fun2_Filter:
    input:
        "post/filtered.seqs.fasta",
        "post/filtered.seqTab.RDS"
    output:
        "post/tax4fun2/KOs_per_ASV.txt",
        "post/tax4fun2/pathway_per_ASV.txt",
        "post/tax4fun2/functional_prediction.txt",
        "post/tax4fun2/pathway_prediction.txt"
    threads: config['bigCores']
    resources:
        runtime="24:00:00",
        mem=config['bigMem']
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
    physInputs = ["post/filtered.seqTab.RDS","reporting/finalNumbers_perSample.tsv"]
    if config['postprocessing']['treeing']['do']:
        physInputs.append("post/ASV.tree.newick")
    if 'postclustering' in STEPS:
        physInputs.append("clusteredTables/cluster_info.tsv")
    rule phyloseq_handoff_postFilter:
        input:
            physInputs
        output:
            "post/filtered.seqTab.phyloseq.RDS"
        threads: config['bigCores']
        params:
            currentStep = "post",
            what = "ASV",
            tree = "add" if config['postprocessing']['treeing']['do'] else "none",
            cluster = "add" if 'postclustering' in STEPS else "none",
            clusterMeth = config['post_clustering']['method']
        resources:
            runtime="12:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "add_R_env.yml"
        log: "logs/phyloseq_hand-off.log"
        script:
            SCRIPTSDIR+"phyloseq_handoff.R"

    physInputsCl = ["post/filtered.clusteredTab.RDS","reporting/finalNumbers_perSample.tsv"]
    if config['postprocessing']['treeing']['do']:
        physInputsCl.append("post/cluster.tree.newick")
    rule phyloseq_handoff_postFilterCl:
        input:
            physInputsCl
        output:
            "post/filtered.clusteredTab.phyloseq.RDS"
        threads: config['bigCores']
        params:
            currentStep = "post",
            what = "cluster",
            tree = "add" if config['postprocessing']['treeing']['do'] else "none",
            cluster = "none"
        resources:
            runtime="12:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "add_R_env.yml"
        log: "logs/phyloseq_hand-off.cluster.log"
        script:
            SCRIPTSDIR+"phyloseq_handoff.R"

rule multiAlign_Filter:
    input:
        "post/filtered.{seqs}.fasta"
    output:
        "post/filtered.{seqs}.multi.fasta"
    threads: config['bigCores']
    resources:
        runtime="48:00:00",
        mem=config['bigMem']
    conda: ENVDIR + "dadasnake_env.yml"
    log: "logs/treeing_multiAlign.{seqs}.log"
    shell:
        """
        clustalo -i {input} -o {output} --outfmt=fasta --threads={threads} --force &> {log} || touch {output}
        """

if config['postprocessing']['treeing']['fasttreeMP'] != "":
    rule treeing_Filter_fasttreeMP:
        input:
            lambda wildcards: "post/filtered.seqs.multi.fasta" if wildcards.kind=="ASV" else "post/filtered.consensus.multi.fasta"
        output:
            "post/{kind}.tree.newick"
        params:
            ftMP=config['postprocessing']['treeing']['fasttreeMP']
        threads: config['bigCores']
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/treeing.{kind}.log"
        shell:
            """
            {params.ftMP} -nt -gamma -no2nd -fastest -spr 4 \
             -log {log} -quiet {input} > {output} 2>> {log} || touch {output}
            """
else:
    rule treeing_Filter:
        input:
            lambda wildcards: "post/filtered.seqs.multi.fasta" if wildcards.kind=="ASV" else "post/filtered.consensus.multi.fasta"
        output:
            "post/{kind}.tree.newick"
        threads: 1
        resources:
            runtime="120:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/treeing.{kind}.log"
        shell:
            """
            fasttree -nt -gamma -no2nd -fastest -spr 4 \
             -log {log} -quiet {input} > {output} 2>> {log} || touch {output}
            """


