localrules: dada_control

rule dada_control:
    input:
        "sequenceTables/all.seqTab.RDS",
        "sequenceTables/all.seqs.fasta",
        "reporting/finalNumbers_perSample.tsv",
        expand("stats/QC_{step}.{run}.{direction}.pdf",step=['1','filtered'],direction=['fwd','rvs'],run=samples.run.unique())  
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
        expand("filtered/{samples.run}/{samples.sample}.{direction}.fastq.gz", samples=samples.itertuples(), direction=["fwd","rvs"])
    output:
        report("reporting/filteredNumbers_perLibrary.tsv",category="Reads"),
        report("reporting/filteredNumbers_perSample.tsv",category="Reads")
    threads: 1
    params:
        currentStep = "filtered",
        mem="8G",
        runtime="12:00:00"
    conda: ENVDIR + "dada_env.yml"
    log: "logs/countFilteredReads.log"
    script:
        SCRIPTSDIR+"report_readNumbers.R"

rule merged_numbers:
    input:
        "reporting/filteredNumbers_perLibrary.tsv",
        expand("merged/{samples.run}/{samples.sample}.RDS", samples=samples.itertuples())
    output:
        report("reporting/mergedNumbers_perLibrary.tsv",category="Reads"),
        report("reporting/mergedNumbers_perSample.tsv",category="Reads")
    threads: 1
    params:
        currentStep = "merged",
        mem="8G",
        runtime="12:00:00"
    conda: ENVDIR + "dada_env.yml"
    log: "logs/countMergedReads.log"
    script:
        SCRIPTSDIR+"report_readNumbers.R"

rule dada_qc1:
    input:
        lambda wildcards: get_sample_perRun(wildcards,"preprocessing/{run}/",".fwd.fastq"),
        lambda wildcards: get_sample_perRun(wildcards,"preprocessing/{run}/",".rvs.fastq")
    output:
        report("stats/QC_1.{run}.fwd.pdf"),
        report("stats/QC_1.{run}.rvs.pdf")
    threads: 1
    params:
        path="preprocessing/{run}",
        mem="8G",
        runtime="12:00:00"
    conda: ENVDIR + "dada_env.yml"
    log: "logs/DADA2_QC_1.{run}.log"
    message: "Running QC on {params.path}."
    script:
        SCRIPTSDIR+"dada_QC.R"

rule dada_qc_filtered:
    input:
        lambda wildcards: get_sample_perRun(wildcards,"filtered/{run}/",".fwd.fastq.gz"),
        lambda wildcards: get_sample_perRun(wildcards,"filtered/{run}/",".rvs.fastq.gz")
    output:
        report("stats/QC_filtered.{run}.fwd.pdf"),
        report("stats/QC_filtered.{run}.rvs.pdf")
    threads: 1
    params:
        path="filtered/{run}",
        mem="8G",
        runtime="12:00:00"
    conda: ENVDIR + "dada_env.yml"
    log: "logs/DADA2_QC_filtered.{run}.log"
    message: "Running QC on {params.path}."
    script:
        SCRIPTSDIR+"dada_QC.R"

rule dada_filter:
    input:
        "preprocessing/{run}/{sample}.fwd.fastq",
        "preprocessing/{run}/{sample}.rvs.fastq"
    output:
        "filtered/{run}/{sample}.fwd.fastq.gz",
        "filtered/{run}/{sample}.rvs.fastq.gz"
    threads: 1
    params:
        mem="8G",
        runtime="12:00:00"
    conda: ENVDIR + "dada_env.yml"
    log: "logs/DADA2_filtering.{run}.{sample}.log"
    message: "Running filtering on {input}."
    script:
        SCRIPTSDIR+"dada_filter.R"


rule dada_errors:
    input:
        lambda wildcards: get_sample_perRun(wildcards,"filtered/{run}/",".{direction}.fastq.gz")
    output:
        "errors/models.{run}.{direction}.RDS",
        "stats/error_models.{run}.{direction}.pdf",
    threads: 1
    params:
        mem="8G",
        runtime="12:00:00"
    conda: ENVDIR + "dada_env.yml"
    log: "logs/DADA2_errors.{run}.{direction}.log"
    message: "Running error models on {input}."
    script:
        SCRIPTSDIR+"dada_errors.R"

if config['dada']['use_quals']:
    rule dada_mergeReadPairs:
        input:
            "errors/models.{run}.fwd.RDS",
            "errors/models.{run}.rvs.RDS",
            "filtered/{run}/{sample}.fwd.fastq.gz",
            "filtered/{run}/{sample}.rvs.fastq.gz"
        output:
            "merged/{run}/{sample}.RDS"
        threads: 1
        params:
            mem="8G",
            runtime="12:00:00"
        conda: ENVDIR + "dada_env.yml"
        log: "logs/DADA2_mergeReadPairs.{run}.{sample}.log"
        message: "merging reads for {wildcards.run} {wildcards.sample}."
        script:
            SCRIPTSDIR+"dada_dadaMergeReads.R"
else:
    rule dada_mergeReadPairs:
        input:
            "filtered/{run}/{sample}.fwd.fastq.gz", 
            "filtered/{run}/{sample}.rvs.fastq.gz"
        output:
            "merged/{run}/{sample}.RDS"
        threads: 1
        params:
            mem="8G",
            runtime="12:00:00"
        conda: ENVDIR + "dada_env.yml"
        log: "logs/DADA2_mergeReadPairs.{run}.{sample}.log"
        message: "merging reads for {wildcards.run} {wildcards.sample}."
        script:
            SCRIPTSDIR+"dada_dadaMergeReads.noError.R"


rule dada_mergeSamples:
    input:
        lambda wildcards: get_sample_perRun(wildcards,"merged/{run}/",".RDS"),
    output:
        "merged/dada_merged.{run}.RDS",
        "sequenceTables/seqTab.{run}.RDS",
        "sequenceTables/seqTab.{run}.tsv"
    threads: 1
    params:
        mem="30G",
        runtime="12:00:00"
    conda: ENVDIR + "dada_env.yml"
    log: "logs/DADA2_mergeSamples.{run}.log"
    message: "preparing sequence table for {wildcards.run}."
    script:
        SCRIPTSDIR+"dada_gatherMergedReads.R"

if config["chimeras"]["remove"]:
    rule dada_mergeruns:
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
        threads: 1
        params:
            mem="8G",
            runtime="12:00:00"
        conda: ENVDIR + "dada_env.yml"
        log: "logs/DADA2_mergeRuns.log"
        message: "merging runs and removing chimeras for {input}."
        script:
            SCRIPTSDIR+"dada_mergeRuns.R"

    rule nochime_numbers:
        input:
            "reporting/mergedNumbers_perSample.tsv",
            "sequenceTables/pre_chimera.seqTab.RDS",
            "sequenceTables/all.seqTab.originalFormat.RDS"
        output:
            report("reporting/finalNumbers_perSample.tsv",category="Reads")
        threads: 1
        params:
            currentStep = "chimera",
            mem="8G",
            runtime="12:00:00"
        conda: ENVDIR + "dada_env.yml"
        log: "logs/countNonchimericReads.log"
        script:
            SCRIPTSDIR+"report_readNumbers.R"
else:
    rule dada_mergeruns:
        input:
            expand("sequenceTables/seqTab.{run}.RDS",run=samples.run.unique())
        output:
            "sequenceTables/all.seqTab.originalFormat.RDS",
            "sequenceTables/all.seqTab.RDS",
            "sequenceTables/all.seqs.fasta",
            "sequenceTables/all.seqTab.tsv"
        threads: 1
        params:
            mem="8G",
            runtime="12:00:00"
        conda: ENVDIR + "dada_env.yml"
        log: "logs/DADA2_mergeRuns.log"
        message: "merging runs for {input}."
        script:
            SCRIPTSDIR+"dada_mergeRuns.R"

    rule tabled_numbers:
        input:
            "reporting/mergedNumbers_perSample.tsv",
            "sequenceTables/all.seqTab.originalFormat.RDS"
        output:
            report("reporting/finalNumbers_perSample.tsv",category="Reads")
        threads: 1
        params:
            currentStep = "table",
            mem="8G",
            runtime="12:00:00"
        conda: ENVDIR + "dada_env.yml"
        log: "logs/countTabledReads.log"
        script:
            SCRIPTSDIR+"report_readNumbers.R"

