localrules: dada_control

rule dada_control:
    input:
        "sequenceTables/all.seqTab.RDS",
        "sequenceTables/all.seqs.fasta",
        "reporting/finalNumbers_perSample.tsv",
        expand("stats/QC_{step}.{run}.pdf",step=['1','filtered'],run=samples.run.unique())  
    output:
        "dada.done"
    shell:
        """
        touch {output}
        """

def get_sample_perRun(wildcards,prefix,suffix):
    return prefix+samples.loc[samples['run']==wildcards.run, "sample"].unique()+suffix

rule filter_numbers:
    input:
        "reporting/primerNumbers_perLibrary.tsv",
        expand("filtered/{samples.run}/{samples.sample}.fastq.gz", samples=samples.itertuples())
    output:
        report("reporting/filteredNumbers_perLibrary.tsv",category="Reads"),
        report("reporting/filteredNumbers_perSample.tsv",category="Reads")
    threads: 1
    params:
        currentStep = "filtered"
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/countFilteredReads.log"
    script:
        SCRIPTSDIR+"report_readNumbers.single.R"

rule merged_numbers:
    input:
        "reporting/filteredNumbers_perLibrary.tsv",
        "merged/dada_merged.RDS"
    output:
        report("reporting/mergedNumbers_perLibrary.tsv",category="Reads"),
        report("reporting/mergedNumbers_perSample.tsv",category="Reads")
    threads: 1
    params:
        currentStep = "merged",
        pooling = config['dada']['pool']
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/countMergedReads.log"
    script:
        SCRIPTSDIR+"report_readNumbers.single.R"

rule dada_qc1:
    input:
        lambda wildcards: get_sample_perRun(wildcards,"preprocessing/{run}/",".fastq")
    output:
        report("stats/QC_1.{run}.pdf")
    threads: 1
    params:
        path="preprocessing/{run}"
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/DADA2_QC_1.{run}.log"
    message: "Running QC on {params.path}."
    script:
        SCRIPTSDIR+"dada_QC.single.R"

rule dada_qc_filtered:
    input:
        lambda wildcards: get_sample_perRun(wildcards,"filtered/{run}/",".fastq.gz")
    output:
        report("stats/QC_filtered.{run}.pdf")
    threads: 1
    params:
        path="filtered/{run}"
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/DADA2_QC_filtered.{run}.log"
    message: "Running QC on {params.path}."
    script:
        SCRIPTSDIR+"dada_QC.single.R"

rule dada_filter:
    input:
        "preprocessing/{run}/{sample}.fastq"
    output:
        "filtered/{run}/{sample}.fastq.gz"
    threads: 1
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/DADA2_filtering.{run}.{sample}.log"
    message: "Running filtering on {input}."
    script:
        SCRIPTSDIR+"dada_filter.single.R"

rule dada_errors:
    input:
        expand("filtered/{samples.run}/{samples.sample}.fastq.gz", samples=samples.itertuples())
    output:
        "errors/models.RDS",
        "stats/error_models.pdf",
    threads: 1
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/DADA2_errors.log"
    message: "Running error models on {input}."
    script:
        SCRIPTSDIR+"dada_errors.R"

if config['dada']['use_quals']:
    rule dada_dadaSingle_pool:
        input:
            "errors/models.RDS",
            expand("filtered/{samples.run}/{samples.sample}.fastq.gz", samples=samples.itertuples())
        output:
            "merged/dada_merged.RDS"
        threads: getThreads(12)
        resources:
            runtime="24:00:00",
            mem=config['normalMem']
        params:
            pooling=config['dada']['pool']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_read2RDS.log"
        message: "converting fastq to dada-RDS."
        script:
            SCRIPTSDIR+"dada_dadaReads.pool.R"
else:
    rule dada_dadaSingle_pool:
        input:
            expand("filtered/{samples.run}/{samples.sample}.fastq.gz", samples=samples.itertuples())
        output:
            "merged/dada_merged.RDS"
        threads: getThreads(12)
        resources:
            runtime="24:00:00",
            mem=config['normalMem']
        params:
            pooling=config['dada']['pool']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_read2RDS.log"
        message: "converting fastq to dada-RDS."
        script:
            SCRIPTSDIR+"dada_dadaReads.pool.noError.R"



if config["chimeras"]["remove"]:
    rule dada_poolTabs:
        input:
            "merged/dada_merged.RDS"
        output:
            "sequenceTables/all.seqTab.originalFormat.RDS",
            "sequenceTables/all.seqTab.RDS",
            "sequenceTables/all.seqs.fasta",
            "sequenceTables/all.seqTab.tsv",
            "sequenceTables/pre_chimera.seqTab.RDS",
            "sequenceTables/pre_chimera.seqs.fasta",
            "sequenceTables/pre_chimera.seqTab.tsv"
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_poolTabs.log"
        message: "removing chimeras for {input}."
        script:
            SCRIPTSDIR+"dada_poolTabs.R"

    rule nochime_numbers:
        input:
            "reporting/filteredNumbers_perSample.tsv",
            "sequenceTables/pre_chimera.seqTab.RDS",
            "sequenceTables/all.seqTab.originalFormat.RDS"
        output:
            report("reporting/finalNumbers_perSample.tsv",category="Reads")
        threads: 1
        params:
            currentStep = "chimera"
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/countNonchimericReads.log"
        script:
            SCRIPTSDIR+"report_readNumbers.single.R"
else:
    rule dada_poolTabs:
        input:
            "merged/dada_merged.RDS"
        output:
            "sequenceTables/all.seqTab.originalFormat.RDS",
            "sequenceTables/all.seqTab.RDS",
            "sequenceTables/all.seqs.fasta",
            "sequenceTables/all.seqTab.tsv"
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_poolTabs.log"
        message: "writing tables and fasta file for {input}."
        script:
            SCRIPTSDIR+"dada_poolTabs.R"

    rule tabled_numbers:
        input:
            "reporting/filteredNumbers_perSample.tsv",
            "sequenceTables/all.seqTab.originalFormat.RDS"
        output:
            report("reporting/finalNumbers_perSample.tsv",category="Reads")
        threads: 1
        params:
            currentStep = "table"
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/countTabledReads.log"
        script:
            SCRIPTSDIR+"report_readNumbers.single.R"

