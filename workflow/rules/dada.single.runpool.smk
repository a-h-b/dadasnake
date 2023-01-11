localrules: dada_control

rule dada_control:
    input:
        "sequenceTables/all.seqTab.RDS",
        "sequenceTables/all.seqs.fasta",
        "reporting/finalNumbers_perSample.tsv",
        expand("stats/QC_{step}.{run}.pdf",step=['1','filtered'],run=samples.run.unique()),
        expand("stats/multiqc_{step}_report.html",step=['1','filtered']) 
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
        "reporting/GtailsNumbers_perLibrary.tsv" if  'primers' not in STEPS and config['nextseq_novaseq'] else "reporting/primerNumbers_perLibrary.tsv",
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
        pool = config['dada']['pool']
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/countMergedReads.log"
    script:
        SCRIPTSDIR+"report_readNumbers.single.R"

rule dada_qc1:
    input:
        lambda wildcards: get_sample_perRun(wildcards,"preprocessing/{run}/",".fastq.gz")
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
        "preprocessing/{run}/{sample}.fastq.gz"
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

rule fastqc_1:
    input:
        expand("preprocessing/{samples.run}/{samples.sample}.fastq.gz", samples=samples.itertuples())
    output:
        directory('stats/fastqc_1')
    threads: 1
    resources:
        runtime = "8:00:00",
        mem = config['normalMem']
    conda: ENVDIR + "fastqc.yml"
    log: "logs/fastqc_1.log"
    message: "fastqc_1: Running fastQC on raw reads."
    shell:
        """
        mkdir -p {output[0]}
        fastqc --noextract -o {output[0]} -f fastq -t {threads} -d {TMPDIR} {input} >> {log} 2>&1
        """

rule multiqc:
    input:
        "stats/fastqc_{step}"
    output:
        directory('stats/multiqc_{step}_report_data'),
        "stats/multiqc_{step}_report.html"
    threads: 1
    resources:
        runtime = "8:00:00",
        mem = config['normalMem']
    conda: ENVDIR + "fastqc.yml"
    log: "logs/multiqc_{step}.log"
    message: "multiqc_1: Collecting fastQC on {wildcards.step} reads."
    shell:
        """
        export LC_ALL=en_GB.utf8
        export LANG=en_GB.utf8
        multiqc -n {output[1]} {input} >> {log} 2>&1
        """

rule fastqc_filtered:
    input:
        expand("filtered/{samples.run}/{samples.sample}.fastq.gz", samples=samples.itertuples())
    output:
        directory('stats/fastqc_filtered')
    threads: 1
    resources:
        runtime = "8:00:00",
        mem = config['normalMem']
    params:
        outputdir = 'stats/fastqc_filtered'
    conda: ENVDIR + "fastqc.yml"
    log: "logs/fastqc_filtered.log"
    message: "fastqc_filtered: Running fastQC on filtered reads."
    shell:
        """
        mkdir -p {output[0]}
        fastqc --noextract -o {output[0]} -f fastq -t {threads} -d {TMPDIR} {input} >> {log} 2>&1
        """

if config['downsampling']['do']:
    rule downsampling:
        input:
            "reporting/filteredNumbers_perLibrary.tsv",
            "filtered/{run}/{sample}.fastq.gz"
        output:
            "downsampled/{run}/{sample}.fastq.gz"
        threads: 1
        resources:
            runtime="2:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_downsampling.{run}.{sample}.log"
        message: "Downsampling {input}."
        script:
            SCRIPTSDIR+"dada_downsample.single.R"

    rule dada_errors:
        input:
            lambda wildcards: get_sample_perRun(wildcards,"downsampled/{run}/",".fastq.gz")
        output:
            "errors/models.{run}.RDS",
            "stats/error_models.{run}.pdf",
        threads: 1
        params:
            errorFunctions=SCRIPTSDIR+"errorFunctions.R"
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_errors.{run}.log"
        message: "Running error models on {input}."
        script:
            SCRIPTSDIR+"dada_errors.R"
    if config['dada']['use_quals']:
        rule dada_dadaSingle_pool:
            input:
                "errors/models.{run}.RDS",
                lambda wildcards: get_sample_perRun(wildcards,"downsampled/{run}/",".fastq.gz")
            output:
                "merged/{run}/dada_merged.RDS",
                "merged/{run}/seqTabs.RDS"
            threads: getThreads(12)
            resources:
                runtime="24:00:00",
                mem=config['bigMem']
            params:
                pooling=config['dada']['pool'],
                errorFunctions=SCRIPTSDIR+"errorFunctions.R"
            conda: ENVDIR + "dada2_env.yml"
            log: "logs/DADA2_read2RDS.{run}.log"
            message: "converting fastq to dada-RDS for {wildcards.run}."
            script:
                SCRIPTSDIR+"dada_dadaReads.runpool.R"
    else:
        rule dada_dadaSingle_pool:
            input:
                lambda wildcards: get_sample_perRun(wildcards,"downsampled/{run}/",".fastq.gz")
            output:
                "merged/{run}/dada_merged.RDS",
                "merged/{run}/seqTabs.RDS"
            threads: getThreads(12)
            resources:
                runtime="24:00:00",
                mem=config['normalMem']
            params:
                pooling=config['dada']['pool']
            conda: ENVDIR + "dada2_env.yml"
            log: "logs/DADA2_read2RDS.{run}.log"
            message: "converting fastq to dada-RDS {wildcards.run}."
            script:
                SCRIPTSDIR+"dada_dadaReads.single.runpool.noError.R"
else:
    rule dada_errors:
        input:
            lambda wildcards: get_sample_perRun(wildcards,"filtered/{run}/",".fastq.gz")
        output:
            "errors/models.{run}.RDS",
            "stats/error_models.{run}.pdf",
        threads: 1
        params:
            errorFunctions=SCRIPTSDIR+"errorFunctions.R"
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_errors.{run}.log"
        message: "Running error models on {input}."
        script:
            SCRIPTSDIR+"dada_errors.R"
    if config['dada']['use_quals']:
        rule dada_dadaSingle_pool:
            input:
                "errors/models.{run}.RDS",
                lambda wildcards: get_sample_perRun(wildcards,"filtered/{run}/",".fastq.gz") 
            output:
                "merged/{run}/dada_merged.RDS",
                "merged/{run}/seqTabs.RDS"
            threads: getThreads(12)
            resources:
                runtime="24:00:00",
                mem=config['bigMem']
            params:
                pooling=config['dada']['pool'],
                errorFunctions=SCRIPTSDIR+"errorFunctions.R"
            conda: ENVDIR + "dada2_env.yml"
            log: "logs/DADA2_read2RDS.{run}.log"
            message: "converting fastq to dada-RDS for {wildcards.run}."
            script:
                SCRIPTSDIR+"dada_dadaReads.runpool.R"
    else:
        rule dada_dadaSingle_pool:
            input:
                lambda wildcards: get_sample_perRun(wildcards,"filtered/{run}/",".fastq.gz")
            output:
                "merged/{run}/dada_merged.RDS",
                "merged/{run}/seqTabs.RDS"
            threads: getThreads(12)
            resources:
                runtime="24:00:00",
                mem=config['normalMem']
            params:
                pooling=config['dada']['pool']
            conda: ENVDIR + "dada2_env.yml"
            log: "logs/DADA2_read2RDS.{run}.log"
            message: "converting fastq to dada-RDS {wildcards.run}."
            script:
                SCRIPTSDIR+"dada_dadaReads.single.runpool.noError.R"


if config["chimeras"]["remove"]:
    rule dada_mergeruns:
        input:
            expand("merged/{run}/seqTabs.RDS",run=samples.run.unique())
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
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_poolTabs.log"
        message: "removing chimeras for {input}."
        script:
            SCRIPTSDIR+"dada_mergeRuns.R"

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
    rule dada_mergeruns:
        input:
            expand("merged/{run}/seqTabs.RDS",run=samples.run.unique())
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
        log: "logs/DADA2_poolTabs.log"
        message: "writing tables and fasta file for {input}."
        script:
            SCRIPTSDIR+"dada_mergeRuns.R"

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

