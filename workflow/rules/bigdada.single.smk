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
        expand("merged/{samples.run}/{samples.sample}.RDS", samples=samples.itertuples())
    output:
        report("reporting/mergedNumbers_perLibrary.tsv",category="Reads"),
        report("reporting/mergedNumbers_perSample.tsv",category="Reads")
    threads: 1
    params:
        currentStep = "merged"
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
        lambda wildcards: get_sample_perRun(wildcards,"filtered/{run}/",".fastq.gz")
    output:
        "errors/models.{run}.RDS",
        "stats/error_models.{run}.pdf",
    threads: 1
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/DADA2_errors.{run}.log"
    message: "Running error models on {input}."
    script:
        SCRIPTSDIR+"dada_errors.R"

if config['dada']['use_quals']:
    rule dada_dadaSingle:
        input:
            "errors/models.{run}.RDS",
            "filtered/{run}/{sample}.fastq.gz"
        output:
            "merged/{run}/{sample}.RDS"
        threads: 1
        resources:
            runtime="24:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_read2RDS.{run}.{sample}.log"
        message: "converting fastq to dada-RDS for {wildcards.run} {wildcards.sample}."
        script:
            SCRIPTSDIR+"dada_dadaReads.single.R"
else:
    rule dada_dadaSingle:
        input:
            "filtered/{run}/{sample}.fastq.gz"
        output:
            "merged/{run}/{sample}.RDS"
        threads: 1
        resources:
            runtime="24:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_read2RDS.{run}.{sample}.log"
        message: "converting fastq to dada-RDS for {wildcards.run} {wildcards.sample}."
        script:
            SCRIPTSDIR+"dada_dadaReads.single.noError.R"


rule dada_mergeSamples:
    input:
        lambda wildcards: get_sample_perRun(wildcards,"merged/{run}/",".RDS"),
    output:
        "merged/dada_merged.{run}.RDS",
        "sequenceTables/seqTab.{run}.RDS",
        "sequenceTables/seqTab.{run}.tsv"
    threads: 1
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/DADA2_mergeSamples.{run}.log"
    message: "preparing sequence table for {wildcards.run}."
    script:
        SCRIPTSDIR+"dada_gatherMergedReads.R"

if config["chimeras"]["remove"]:
    rule bigdada_mergeruns:
        input:
            expand("sequenceTables/seqTab.{run}.RDS",run=samples.run.unique())
        output:
            "sequenceTables/all.seqTab.originalFormat.RDS",
            "sequenceTables/all.seqTab.RDS",
            "sequenceTables/all.seqs.fasta",
            "sequenceTables/all.seqTab.tsv",
            "sequenceTables/pre_chimera.seqTab.RDS",
            "sequenceTables/pre_chimera.seqs.fasta",
            "sequenceTables/pre_chimera.seqTab.tsv"
        threads: config['bigCores']
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/bigDADA2_mergeRuns.log"
        message: "merging runs and removing chimeras for {input}."
        script:
            SCRIPTSDIR+"dada_mergeRuns.R"

    rule bignochime_numbers:
        input:
            "reporting/mergedNumbers_perSample.tsv",
            "sequenceTables/pre_chimera.seqTab.RDS",
            "sequenceTables/all.seqTab.originalFormat.RDS"
        output:
            report("reporting/finalNumbers_perSample.tsv",category="Reads")
        threads: config['bigCores']
        params:
            currentStep = "chimera"
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/bignochime_numbers.log"
        script:
            SCRIPTSDIR+"report_readNumbers.single.R"

else:
    rule bigdada_mergeruns:
        input:
            expand("sequenceTables/seqTab.{run}.RDS",run=samples.run.unique())
        output:
            "sequenceTables/all.seqTab.originalFormat.RDS",
            "sequenceTables/all.seqTab.RDS",
            "sequenceTables/all.seqs.fasta",
            "sequenceTables/all.seqTab.tsv"
        threads: config['bigCores']
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/bigdada_mergeruns.log"
        message: "merging runs for {input}."
        script:
            SCRIPTSDIR+"dada_mergeRuns.R"

    rule bigtabled_numbers:
        input:
            "reporting/mergedNumbers_perSample.tsv",
            "sequenceTables/all.seqTab.originalFormat.RDS"
        output:
            report("reporting/finalNumbers_perSample.tsv",category="Reads")
        threads: config['bigCores']
        params:
            currentStep = "table"
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/bigtabled_numbers.log"
        script:
            SCRIPTSDIR+"report_readNumbers.single.R"

