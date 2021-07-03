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
        currentStep = "filtered"
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/countFilteredReads.log"
    script:
        SCRIPTSDIR+"report_readNumbers.R"

rule merged_numbers:
    input:
        "reporting/filteredNumbers_perLibrary.tsv",
        expand("merged/dada.RDS", samples=samples.itertuples())
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
        SCRIPTSDIR+"report_readNumbers.R"

rule dada_qc1:
    input:
        lambda wildcards: get_sample_perRun(wildcards,"preprocessing/{run}/",".fwd.fastq.gz"),
        lambda wildcards: get_sample_perRun(wildcards,"preprocessing/{run}/",".rvs.fastq.gz")
    output:
        report("stats/QC_1.{run}.fwd.pdf"),
        report("stats/QC_1.{run}.rvs.pdf")
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
        path="filtered/{run}"
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/DADA2_QC_filtered.{run}.log"
    message: "Running QC on {params.path}."
    script:
        SCRIPTSDIR+"dada_QC.R"


rule fastqc_1:
    input:
        expand("preprocessing/{samples.run}/{samples.sample}.{{dir}}.fastq.gz", samples=samples.itertuples())
    output:
        directory('stats/fastqc_1_{dir}'),
        directory('stats/multiqc_1_{dir}_data'),
        "stats/multiqc_1_{dir}_report.html"
    threads: getThreads(6)
    resources:
        runtime = "8:00:00",
        mem = config['normalMem']
    params:
        outputdir = 'stats/fastqc_1_{dir}'
    conda: ENVDIR + "fastqc.yml"
    log: "logs/fastqc_1_{dir}.log"
    message: "fastqc_1: Running fastQC on raw {wildcards.dir} reads."
    shell:
        """
        mkdir -p {params.outputdir}
        fastqc --noextract -o {params.outputdir} -f fastq {input} -t {threads} -d {TMPDIR} > {log} 2>&1
        multiqc -n {output[2]} -o {output[1]} {params.outputdir} >> {log} 2>&1
        """

rule fastqc_filtered:
    input:
        expand("filtered/{samples.run}/{samples.sample}.{{dir}}.fastq.gz", samples=samples.itertuples())
    output:
        directory('stats/fastqc_filtered_{dir}'),
        directory('stats/multiqc_filtered_{dir}_data'),
        "stats/multiqc_filtered_{dir}_report.html"
    threads: getThreads(6)
    resources:
        runtime = "8:00:00",
        mem = config['normalMem']
    params:
        outputdir = 'stats/fastqc_filtered_{dir}'
    conda: ENVDIR + "fastqc.yml"
    log: "logs/fastqc_filtered_{dir}.log"
    message: "fastqc_filtered: Running fastQC on filtered {wildcards.dir} reads."
    shell:
        """
        mkdir -p {params.outputdir}
        fastqc --noextract -o {params.outputdir} -f fastq {input} -t {threads} -d {TMPDIR} > {log} 2>&1
        multiqc -n {output[2]} -o {output[1]} {params.outputdir} >> {log} 2>&1
        """

rule dada_filter:
    input:
        "preprocessing/{run}/{sample}.fwd.fastq.gz",
        "preprocessing/{run}/{sample}.rvs.fastq.gz"
    output:
        "filtered/{run}/{sample}.fwd.fastq.gz",
        "filtered/{run}/{sample}.rvs.fastq.gz"
    threads: 1
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/DADA2_filtering.{run}.{sample}.log"
    message: "Running filtering on {input}."
    script:
        SCRIPTSDIR+"dada_filter.R"

if config['downsampling']['do']:
    rule downsampling:
        input:
            "reporting/filteredNumbers_perLibrary.tsv",
            "filtered/{run}/{sample}.fwd.fastq.gz",
            "filtered/{run}/{sample}.rvs.fastq.gz"
        output:
            "downsampled/{run}/{sample}.fwd.fastq.gz",
            "downsampled/{run}/{sample}.rvs.fastq.gz"
        threads: 1
        resources:
            runtime="2:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_downsampling.{run}.{sample}.log"
        message: "Downsampling {input}."
        script:
            SCRIPTSDIR+"dada_downsample.R"

    rule dada_errors_pool:
        input:
            expand("downsampled/{samples.run}/{samples.sample}.{{direction}}.fastq.gz", samples=samples.itertuples())
        output:
            "errors/models.{direction}.RDS",
            "stats/error_models.{direction}.pdf",
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_errors.{direction}.log"
        message: "Running error models on {input}."
        script:
            SCRIPTSDIR+"dada_errors.R"
    if config['dada']['use_quals']:
        rule dada_dadaPairs_pool:
            input:
                "errors/models.{direction}.RDS",
                expand("downsampled/{samples.run}/{samples.sample}.{{direction}}.fastq.gz", samples=samples.itertuples())
            output:
                "merged/derep.{direction}.RDS",
                "merged/dada.{direction}.RDS"
            wildcard_constraints:
                direction = "fwd|rvs"
            threads: getThreads(int(config['bigCores']))
            params:
                pooling=config['dada']['pool']
            resources:
                runtime="48:00:00",
                mem=config['bigMem']
            conda: ENVDIR + "dada2_env.yml"
            log: "logs/DADA2_dada.{direction}.log"
            message: "converting to DADA reads for {wildcards.direction}."
            script:
                SCRIPTSDIR+"dada_dadaReads.paired.pool.R"
    else:
        rule dada_dadaPairs_pool:
            input:
                expand("downsampled/{samples.run}/{samples.sample}.{{direction}}.fastq.gz", samples=samples.itertuples())
            output:
                "merged/derep.{direction}.RDS"
                "merged/dada.{direction}.RDS"
            wildcard_constraints:
                direction = "fwd|rvs"
            threads: getThreads(12)
            params:
                pooling=config['dada']['pool']
            resources:
                runtime="24:00:00",
                mem=config['normalMem']
            conda: ENVDIR + "dada2_env.yml"
            log: "logs/DADA2_dada.{direction}.log"
            message: "converting to DADA reads for {wildcards.direction}."
            script:
                SCRIPTSDIR+"dada_dadaReads.paired.pool.noError.R"
else:
    rule dada_errors_pool:
        input:
            expand("filtered/{samples.run}/{samples.sample}.{{direction}}.fastq.gz", samples=samples.itertuples())
        output:
            "errors/models.{direction}.RDS",
            "stats/error_models.{direction}.pdf",
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_errors.{direction}.log"
        message: "Running error models on {input}."
        script:
            SCRIPTSDIR+"dada_errors.R"
    if config['dada']['use_quals']:
        rule dada_dadaPairs_pool:
            input:
                "errors/models.{direction}.RDS",
                expand("filtered/{samples.run}/{samples.sample}.{{direction}}.fastq.gz", samples=samples.itertuples())
            output:
                "merged/derep.{direction}.RDS",
                "merged/dada.{direction}.RDS"
            wildcard_constraints:
                direction = "fwd|rvs"
            threads: getThreads(int(config['bigCores']))
            params:
                pooling=config['dada']['pool']
            resources:
                runtime="48:00:00",
                mem=config['bigMem']
            conda: ENVDIR + "dada2_env.yml"
            log: "logs/DADA2_dada.{direction}.log"
            message: "converting to DADA reads for {wildcards.direction}."
            script:
                SCRIPTSDIR+"dada_dadaReads.paired.pool.R"
    else:
        rule dada_dadaPairs_pool:
            input:
                expand("filtered/{samples.run}/{samples.sample}.{{direction}}.fastq.gz", samples=samples.itertuples())
            output:
                "merged/derep.{direction}.RDS"
                "merged/dada.{direction}.RDS"
            wildcard_constraints:
                direction = "fwd|rvs"
            threads: getThreads(12)
            params:
                pooling=config['dada']['pool']
            resources:
                runtime="24:00:00",
                mem=config['normalMem']
            conda: ENVDIR + "dada2_env.yml"
            log: "logs/DADA2_dada.{direction}.log"
            message: "converting to DADA reads for {wildcards.direction}."
            script:
                SCRIPTSDIR+"dada_dadaReads.paired.pool.noError.R"


rule dada_mergeReadPairs_pool:
    input:
        "merged/dada.fwd.RDS",
        "merged/dada.rvs.RDS",
        "merged/derep.fwd.RDS",
        "merged/derep.rvs.RDS"
    output:
        "merged/dada.RDS"
    threads: 2
    resources:
        runtime="18:00:00",
        mem=config['bigMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/DADA2_mergeReadPairs.log"
    message: "merging reads."
    script:
        SCRIPTSDIR+"dada_onlyMergePoolReads.R"

if config["chimeras"]["remove"]:
    rule dada_mergeruns:
        input:
            "merged/dada.RDS"
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
            runtime="24:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_mergeRuns.log"
        message: "merging runs and removing chimeras for {input}."
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
            SCRIPTSDIR+"report_readNumbers.R"
else:
    rule dada_mergeruns:
        input:
            "merged/dada.RDS"
        output:
            "sequenceTables/all.seqTab.originalFormat.RDS",
            "sequenceTables/all.seqTab.RDS",
            "sequenceTables/all.seqs.fasta",
            "sequenceTables/all.seqTab.tsv"
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_mergeRuns.log"
        message: "merging runs for {input}."
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
            SCRIPTSDIR+"report_readNumbers.R"

