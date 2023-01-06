localrules: dada_control

rule dada_control:
    input:
        "sequenceTables/all.seqTab.RDS",
        "sequenceTables/all.seqs.fasta",
        "reporting/finalNumbers_perSample.tsv",
        expand("stats/QC_{step}.{run}.{direction}.pdf",step=['1','filtered'],direction=['fwd','rvs'],run=samples.run.unique()),
        expand("stats/multiqc_{step}_{direction}_report.html",step=['1','filtered'],direction=['fwd','rvs']) 
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
        directory('stats/fastqc_1_{dir}')
    threads: 1
    resources:
        runtime = "8:00:00",
        mem = config['normalMem']
    conda: ENVDIR + "fastqc.yml"
    log: "logs/fastqc_1_{dir}.log"
    message: "fastqc_1: Running fastQC on raw {wildcards.dir} reads."
    shell:
        """
        mkdir -p {output[0]}
        fastqc --noextract -o {output[0]} -f fastq -t {threads} -d {TMPDIR} {input} >> {log} 2>&1
        """

rule multiqc:
    input:
        "stats/fastqc_{step}_{dir}"
    output:
        directory('stats/multiqc_{step}_{dir}_report_data'),
        "stats/multiqc_{step}_{dir}_report.html"
    wildcard_constraints:
        dir="rvs|fwd"
    threads: 1
    resources:
        runtime = "8:00:00",
        mem = config['normalMem']
    conda: ENVDIR + "fastqc.yml"
    log: "logs/multiqc_{step}_{dir}.log"
    message: "multiqc_1: Collecting fastQC on {wildcards.step} {wildcards.dir} reads."
    shell:
        """
        export LC_ALL=en_GB.utf8
        export LANG=en_GB.utf8
        multiqc -n {output[1]} {input} >> {log} 2>&1
        """

rule fastqc_filtered:
    input:
        expand("filtered/{samples.run}/{samples.sample}.{{dir}}.fastq.gz", samples=samples.itertuples())
    output:
        directory('stats/fastqc_filtered_{dir}')
    threads: 1
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
        mkdir -p {output[0]}
        fastqc --noextract -o {output[0]} -f fastq -t {threads} -d {TMPDIR} {input} >> {log} 2>&1
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

    rule dada_errors:
        input:
            lambda wildcards: get_sample_perRun(wildcards,"downsampled/{run}/",".{direction}.fastq.gz")
        output:
            "errors/models.{run}.{direction}.RDS",
            "stats/error_models.{run}.{direction}.pdf",
        threads: 1
        params:
            errorFunctions=SCRIPTSDIR+"errorFunctions.R"
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_errors.{run}.{direction}.log"
        message: "Running error models on {input}."
        script:
            SCRIPTSDIR+"dada_errors.R"

    if config['dada']['use_quals']:
        rule dada_dadaPairs_runpool:
            input:
                "errors/models.{run}.{direction}.RDS",
                lambda wildcards: get_sample_perRun(wildcards,"downsampled/{run}/",".{direction}.fastq.gz")
            output:
                "merged/{run}/derep.{direction}.RDS",
                "merged/{run}/dada.{direction}.RDS"
            wildcard_constraints:
                direction = "fwd|rvs"
            threads: getThreads(12)
            params:
                pooling=config['dada']['pool'],
                errorFunctions=SCRIPTSDIR+"errorFunctions.R"
            resources:
                runtime="12:00:00",
                mem=config['normalMem']
            conda: ENVDIR + "dada2_env.yml"
            log: "logs/DADA2_dada.{run}.{direction}.log"
            message: "converting to DADA  reads for {wildcards.run} {wildcards.direction}."
            script:
                SCRIPTSDIR+"dada_dadaReads.paired.runpool.R"
    else:
        rule dada_dadaPairs_runpool:
            input:
                lambda wildcards: get_sample_perRun(wildcards,"downsampled/{run}/",".{direction}.fastq.gz")
            output:
                "merged/{run}/derep.{direction}.RDS",
                "merged/{run}/dada.{direction}.RDS"
            wildcard_constraints:
                direction = "fwd|rvs"
            threads: getThreads(12)
            params:
                pooling=config['dada']['pool']
            resources:
                runtime="12:00:00",
                mem=config['normalMem']
            conda: ENVDIR + "dada2_env.yml"
            log: "logs/DADA2_dada.{run}.{direction}.log"
            message: "converting to DADA reads for {wildcards.run} {wildcards.direction}."
            script:
                SCRIPTSDIR+"dada_dadaReads.runpool.noError.R"
else:
    rule dada_errors:
        input:
            lambda wildcards: get_sample_perRun(wildcards,"filtered/{run}/",".{direction}.fastq.gz")
        output:
            "errors/models.{run}.{direction}.RDS",
            "stats/error_models.{run}.{direction}.pdf",
        threads: 1
        params:
            errorFunctions=SCRIPTSDIR+"errorFunctions.R"
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_errors.{run}.{direction}.log"
        message: "Running error models on {input}."
        script:
            SCRIPTSDIR+"dada_errors.R"

    if config['dada']['use_quals']:
        rule dada_dadaPairs_runpool:
            input:
                "errors/models.{run}.{direction}.RDS",
                lambda wildcards: get_sample_perRun(wildcards,"filtered/{run}/",".{direction}.fastq.gz")
            output:
                "merged/{run}/derep.{direction}.RDS",
                "merged/{run}/dada.{direction}.RDS"
            wildcard_constraints:
                direction = "fwd|rvs"
            threads: getThreads(12)
            params:
                pooling=config['dada']['pool'],
                errorFunctions=SCRIPTSDIR+"errorFunctions.R"
            resources:
                runtime="12:00:00",
                mem=config['normalMem']
            conda: ENVDIR + "dada2_env.yml"
            log: "logs/DADA2_dada.{run}.{direction}.log"
            message: "converting to DADA  reads for {wildcards.run} {wildcards.direction}."
            script:
                SCRIPTSDIR+"dada_dadaReads.paired.runpool.R"
    else:
        rule dada_dadaPairs_runpool:
            input:
                lambda wildcards: get_sample_perRun(wildcards,"filtered/{run}/",".{direction}.fastq.gz")
            output:
                "merged/{run}/derep.{direction}.RDS",
                "merged/{run}/dada.{direction}.RDS"
            wildcard_constraints:
                direction = "fwd|rvs"
            threads: getThreads(12)
            params:
                pooling=config['dada']['pool']
            resources:
                runtime="12:00:00",
                mem=config['normalMem']
            conda: ENVDIR + "dada2_env.yml"
            log: "logs/DADA2_dada.{run}.{direction}.log"
            message: "converting to DADA reads for {wildcards.run} {wildcards.direction}."
            script:
                SCRIPTSDIR+"dada_dadaReads.runpool.noError.R"


rule dada_mergeReadPairs_runpool:
    input:
        "merged/{run}/dada.fwd.RDS",
        "merged/{run}/dada.rvs.RDS",
        "merged/{run}/derep.fwd.RDS",
        "merged/{run}/derep.rvs.RDS"
    output:
        "merged/{run}/dada.RDS",
        "merged/{run}/seqTab.RDS"
    threads: 1
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/DADA2_mergeReadPairs.{run}.log"
    message: "merging reads for {wildcards.run}."
    script:
        SCRIPTSDIR+"dada_onlyMergeRunPoolReads.R"

if config["chimeras"]["remove"]:
    rule dada_mergeruns:
        input:
            expand("merged/{run}/seqTab.RDS",run=samples.run.unique())
        output:
            "sequenceTables/all.seqTab.originalFormat.RDS",
            "sequenceTables/all.seqTab.RDS",
            "sequenceTables/all.seqs.fasta",
            "sequenceTables/all.seqTab.tsv",
            "sequenceTables/pre_chimera.seqTab.RDS",
            "sequenceTables/pre_chimera.seqs.fasta",
            "sequenceTables/pre_chimera.seqTab.tsv"
        threads: getThreads(12)
        resources:
            runtime="48:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/DADA2_mergeRuns.log"
        message: "merging runs and removing chimeras for {input}."
        script:
            SCRIPTSDIR+"dada_mergeRuns_runpools.R"

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
            expand("merged/{run}/seqTab.RDS",run=samples.run.unique())
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
            SCRIPTSDIR+"report_readNumbers.R"

