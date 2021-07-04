postConts = []
if config['hand_off']['phyloseq']:
    postConts.append("sequenceTables/all.seqTab.phyloseq.RDS")
if config['postprocessing']['treeing']['do']:
    postConts.append("post/tree.newick")
if config['postprocessing']['rarefaction_curve']:
    postConts.append("stats/rarefaction_curves.pdf")
if config['postprocessing']['funguild']['do']:
    CLASSIFY=config['postprocessing']['funguild']['classifier']
    if config['taxonomy'][CLASSIFY]['do']:
        postConts.append("post/all.seqTab.guilds.tsv")
if config['postprocessing']['fungalTraits']['do']:
    CLASSIFY=config['postprocessing']['fungalTraits']['classifier']
    if config['taxonomy'][CLASSIFY]['do']:
        postConts.append("post/all.seqTab.traits.RDS")
if config['postprocessing']['tax4fun2']['do']:
    postConts.append("post/tax4fun2/KOs_per_OTU.txt")

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
        "sequenceTables/all.seqTab.RDS",
        "reporting/finalNumbers_perSample.tsv"
    output:
        "stats/rarefaction_curves.pdf"
    threads: 1
    resources:
        runtime="48:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/rarefaction_curve.log"
    script:
        SCRIPTSDIR+"rarefaction_curve.R"



rule guilds_noFilter:
    input:
        "sequenceTables/all.seqTab.tax.tsv"
    output:
        "post/all.seqTab.guilds.tsv"
    threads: 1
    resources:
        runtime="24:00:00",
        mem=config['normalMem']
    params:
        src_path=SCRIPTSDIR
    log: "logs/funguild.log"
    conda: ENVDIR + "dadasnake_env.yml"
    message: "Running funguild on {input}."
    shell:
        """
        {params.src_path}/Guilds_v1.1.local.2.py -otu {input} -output {output} -path_to_db {config[postprocessing][funguild][funguild_db]} -taxonomy_name taxonomy.{config[postprocessing][funguild][classifier]}&> {log} || touch {output}
        """

rule funTraits_noFilter:
    input:
        "sequenceTables/all.seqTab.tax.RDS"
    output:
        "post/all.seqTab.traits.tsv",
        "post/all.seqTab.traits.RDS"
    threads: 1
    resources:
        runtime="6:00:00",
        mem=config['normalMem']
    log: "logs/fungalTraits.log"
    conda: ENVDIR + "dada2_env.yml"
    message: "Adding fungalTraits to {input}."
    script:
        SCRIPTSDIR + "add_fungalTraits.R"

rule tax4fun2_noFilter:
    input:
        "sequenceTables/all.seqs.fasta",
        "sequenceTables/all.seqTab.RDS"
    output:
        "post/tax4fun2/KOs_per_OTU.txt",
        "post/tax4fun2/pathway_per_OTU.txt",
        "post/tax4fun2/functional_prediction.txt",
        "post/tax4fun2/pathway_prediction.txt"
    threads: getThreads(4) 
    resources:
        runtime="6:00:00",
        mem=config['normalMem']
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
    if config['do_taxonomy']:
        physInputs = ["sequenceTables/all.seqTab.tax.RDS","reporting/finalNumbers_perSample.tsv"]
    else:
        physInputs = ["sequenceTables/all.seqTab.RDS","reporting/finalNumbers_perSample.tsv"]
    if config['postprocessing']['treeing']['do']:
        physInputs.append("post/tree.newick")
    rule phyloseq_handoff_post:
        input:
            physInputs
        output:
            "sequenceTables/all.seqTab.phyloseq.RDS"
        threads: 1
        params:
            currentStep = "post"
        resources:
            runtime="4:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "add_R_env.yml"
        log: "logs/phyloseq_hand-off.log"
        script:
            SCRIPTSDIR+"phyloseq_handoff.R"

rule multiAlign_noFilter:
    input:
        "sequenceTables/all.seqs.fasta"
    output:
        "post/all.seqs.multi.fasta"
    threads: getThreads(10)
    resources:
        runtime="24:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dadasnake_env.yml"
    log: "logs/treeing_multiAlign.log"
    shell:
        """
        clustalo -i {input} -o {output} --outfmt=fasta --threads={threads} --force &> {log} || touch {output}
        """

if config['postprocessing']['treeing']['fasttreeMP'] != "":
    rule treeing_noFilter_fasttreeMP:
        input:
            "post/all.seqs.multi.fasta"
        output:
            "post/tree.newick"
        threads: getThreads(10)
        resources:
            runtime="48:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/treeing.log"
        shell:
            """
            {config[postprocessing][treeing][fasttreeMP]} -nt -gamma -no2nd -fastest -spr 4 \
             -log {log} -quiet {input} > {output} 2> {log}  || touch {output}
            """
else:
    rule treeing_noFilter:
        input:
            "post/all.seqs.multi.fasta"
        output:
            "post/tree.newick"
        threads: 1
        resources:
            runtime="48:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/treeing.log"
        shell:
            """
            fasttree -nt -gamma -no2nd -fastest -spr 4 \
             -log {log} -quiet {input} > {output} 2>> {log} || touch {output}
            """




