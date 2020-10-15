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

rule bigrarefaction_curve_noFilter:
    input:
        "sequenceTables/all.seqTab.RDS",
        "reporting/finalNumbers_perSample.tsv"
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



rule bigguilds_noFilter:
    input:
        "sequenceTables/all.seqTab.tax.tsv"
    output:
        "post/all.seqTab.guilds.tsv"
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
        {params.src_path}/Guilds_v1.1.local.2.py -otu {input} -output {output} -path_to_db {config[postprocessing][funguild][funguild_db]} -taxonomy_name taxonomy.{config[postprocessing][funguild][classifier]}&> {log} || touch {output}
        """


if config['hand_off']['phyloseq']:
    physInputs = ["sequenceTables/all.seqTab.RDS","reporting/finalNumbers_perSample.tsv"]
    if config['postprocessing']['treeing']['do']:
        physInputs.append("post/tree.newick")
    rule bigphyloseq_handoff_post:
        input:
            physInputs
        output:
            "sequenceTables/all.seqTab.phyloseq.RDS"
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


rule bigmultiAlign_noFilter:
    input:
        "sequenceTables/all.seqs.fasta"
    output:
        "post/all.seqs.multi.fasta"
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
    rule bigtreeing_noFilter_fasttreeMP:
        input:
            "post/all.seqs.multi.fasta"
        output:
            "post/tree.newick"
        threads: config['bigCores']
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/bigtreeing.log"
        shell:
            """
            {config[postprocessing][treeing][fasttreeMP]} -nt -gamma -no2nd -fastest -spr 4 \
             -log {log} -quiet {input} > {output} 2>> {log} || touch {output}
            """
else:
    rule bigtreeing_noFilter:
        input:
            "post/all.seqs.multi.fasta"
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




