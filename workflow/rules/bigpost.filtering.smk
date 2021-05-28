postConts = ["post/filtered.seqTab.RDS","reporting/post_finalNumbers_perSample.tsv"]
if config['hand_off']['phyloseq']:
    postConts.append("post/filtered.seqTab.phyloseq.RDS")
if config['postprocessing']['treeing']['do']:
    postConts.append("post/tree.newick")
if config['postprocessing']['rarefaction_curve']:
    postConts.append("stats/rarefaction_curves.pdf")
if config['postprocessing']['funguild']['do']:
    CLASSIFY=config['postprocessing']['funguild']['classifier']
    if config['taxonomy'][CLASSIFY]['do']:
        postConts.append("post/filtered.seqTab.guilds.tsv")
if config['postprocessing']['fungalTraits']['do']:
    CLASSIFY=config['postprocessing']['fungalTraits']['classifier']
    if config['taxonomy'][CLASSIFY]['do']:
        postConts.append("post/filtered.seqTab.traits.RDS")
if config['postprocessing']['tax4fun2']['do']:
    postConts.append("post/tax4fun2/KOs_per_OTU.txt")

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
if config['do_taxonomy'] and (config['taxonomy']['decipher']['do'] or config['taxonomy']['mothur']['do'] or config['taxonomy']['dada']['do']):
    filtTabs.append("sequenceTables/all.seqTab.tax.RDS")
else:
    filtTabs.append("sequenceTables/all.seqTab.RDS")
rule bigfiltering_table:
    input: 
       filtTabs
    output:
       "post/filtered.seqTab.RDS",
       "post/filtered.seqTab.tsv",   
       "post/filtered.seqs.fasta"  
    threads: config['bigCores']
    resources:
        runtime="12:00:00",
        mem=config['bigMem']
    log: "logs/bigpost_filtering_table.log"
    conda: ENVDIR + "dada2_env.yml"
    script:
        SCRIPTSDIR+"post_filtering.R"

rule bigtable_filter_numbers:
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
    log: "logs/bigcountPostfilteredReads.log"
    script:
        SCRIPTSDIR+"report_readNumbers.R"


rule bigrarefaction_curve_Filter:
    input:
        "post/filtered.seqTab.RDS",
        "reporting/post_finalNumbers_perSample.tsv"
    output:
        "stats/rarefaction_curves.pdf"
    threads: config['bigCores']
    resources:
        runtime="120:00:00",
        mem=config['bigMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/bigrarefaction_curve.log"
    script:
        SCRIPTSDIR+"rarefaction_curve.R"


rule bigguilds_Filter:
    input:
        "post/filtered.seqTab.tsv"
    output:
        "post/filtered.seqTab.guilds.tsv"
    threads: 1
    resources:
        runtime="120:00:00",
        mem=config['bigMem']
    params:
        src_path=SCRIPTSDIR
    log: "logs/funguild.log"
    conda: ENVDIR + "dadasnake_env.yml"
    message: "Running funguild on {input}."
    shell:
        """
        {params.src_path}/Guilds_v1.1.local.2.py -otu {input} -output {output} -path_to_db {config[postprocessing][funguild][funguild_db]} -taxonomy_name taxonomy.{config[postprocessing][funguild][classifier]} &> {log} || touch {output}
        """

rule bigfunTraits_Filter:
    input:
        "post/filtered.seqTab.RDS"
    output:
        "post/filtered.seqTab.traits.tsv",
        "post/filtered.seqTab.traits.RDS"
    threads: 1
    resources:
        runtime="6:00:00",
        mem=config['bigMem']
    log: "logs/fungalTraits.log"
    conda: ENVDIR + "dada2_env.yml"
    message: "Adding fungalTraits to {input}."
    script:
        SCRIPTSDIR + "add_fungalTraits.R"

rule bigtax4fun2_Filter:
    input:
        "post/filtered.seqs.fasta",
        "post/filtered.seqTab.RDS"
    output:
        "post/tax4fun2/KOs_per_OTU.txt",
        "post/tax4fun2/pathway_per_OTU.txt",
        "post/tax4fun2/functional_prediction.txt",
        "post/tax4fun2/pathway_prediction.txt"
    threads: getThreads(4)
    resources:
        runtime="6:00:00",
        mem=config['bigMem']
    params:
        tmp=TMPDIR,
        outputDir="post/tax4fun2",
        customFunc=SCRIPTSDIR + "/functionalPredictionCustom.R"
    log: "logs/tax4fun2.log"
    conda: ENVDIR + "tax4fun2_env.yml"
    message: "Running tax4fun2 on {input}."
    script:
        SCRIPTSDIR + "tax4fun2.R"


if config['hand_off']['phyloseq']:
    physInputs = ["post/filtered.seqTab.RDS","reporting/post_finalNumbers_perSample.tsv"]
    physOutputs = "post/filtered.seqTab.phyloseq.RDS"
    if config['postprocessing']['treeing']['do']:
        physInputs.append("post/tree.newick")
    rule bigphyloseq_handoff_postFilter:
        input:
            physInputs
        output:
            "post/filtered.seqTab.phyloseq.RDS"
        threads: config['bigCores']
        params:
            currentStep = "post"
        resources:
            runtime="12:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "add_R_env.yml"
        log: "logs/bigphyloseq_hand-off.log"
        script:
            SCRIPTSDIR+"phyloseq_handoff.R"

rule bigmultiAlign_Filter:
    input:
        "post/filtered.seqs.fasta"
    output:
        "post/filtered.seqs.multi.fasta"
    threads: config['bigCores']
    resources:
        runtime="12:00:00",
        mem=config['bigMem']
    conda: ENVDIR + "dadasnake_env.yml"
    log: "logs/treeing_multiAlign.log"
    shell:
        """
        clustalo -i {input} -o {output} --outfmt=fasta --threads={threads} --force &> {log} || touch {output}
        """


if config['postprocessing']['treeing']['fasttreeMP'] != "":
    rule bigtreeing_Filter_fasttreeMP:
        input:
            "post/filtered.seqs.multi.fasta"
        output:
            "post/tree.newick"
        threads: config['bigCores']
        resources:
            runtime="12:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/treeing.log"
        shell:
            """
            {config[postprocessing][treeing][fasttreeMP]} -nt -gamma -no2nd -fastest -spr 4 \
             -log {log} -quiet {input} > {output} 2>> {log} || touch {output}
            """
else:
    rule bigtreeing_Filter:
        input:
            "post/filtered.seqs.multi.fasta"
        output:
            "post/tree.newick"
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/treeing.log"
        shell:
            """
            fasttree -nt -gamma -no2nd -fastest -spr 4 \
             -log {log} -quiet {input} > {output} 2>> {log} || touch {output}
            """


