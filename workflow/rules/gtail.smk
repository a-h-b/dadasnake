def get_fastq(wildcards):
    fastqs = config['raw_directory'] + "/" + samples.loc[(wildcards.library,wildcards.run), ["r1_file", "r2_file"]].dropna()
    return fastqs

def get_lib_perRunAndSample(wildcards,prefix,suffix):
    return prefix+samples.loc[(samples['run']==wildcards.run) & (samples['sample']==wildcards.sample), "library"].unique()+suffix


localrules: gtails_control

rule gtails_control:
    input:
        expand("preprocessing/{samples.run}/{samples.sample}.{direction}.fastq.gz", samples=samples.itertuples(), direction=["fwd","rvs"]),
        "reporting/readNumbers.tsv",
        "reporting/GtailsNumbers_perSample.tsv"
    output:
        "gtails.done"
    shell:
        """
        touch {output}
        """

rule combine_or_rename:
    input:
        "reporting/GtailsNumbers_perLibrary.tsv",
        files = lambda wildcards: get_lib_perRunAndSample(wildcards,"preprocessing/{run}/",".{direction}.fastq.gz")    
    output:
        "preprocessing/{run}/{sample}.{direction}.fastq.gz"
    wildcard_constraints:
        direction="(fwd|rvs)",
        sample='|'.join(samples['sample'])
    threads: 1
    log: "logs/combine_or_rename.{run}.{sample}.{direction}.log"
    resources:
        runtime="01:00:00",
        mem=config['normalMem']
    run:
        if len(input) > 2:
            shell("cat {input.files} > {output}")
        else:
            shell("mv {input.files} {output}")

rule input_numbers:
    input:
        "reporting/sample_table.tsv",
        expand("{raw_directory}/{file}", file=samples.r1_file,raw_directory=RAW),
        expand("{raw_directory}/{file}", file=samples.r2_file,raw_directory=RAW)
    output:
        report("reporting/readNumbers.tsv",category="Reads")
    threads: 1
    params:
        currentStep = "raw",
        raw_directory = RAW
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/countInputReads.log"
    script:
        SCRIPTSDIR+"report_readNumbers.R" 


rule gtails_numbers:
    input:
        "reporting/readNumbers.tsv",
        expand("preprocessing/{samples.run}/{samples.library}.{direction}.fastq.gz", samples=samples.itertuples(), direction=["fwd","rvs"])
    output:
        "reporting/GtailsNumbers_perLibrary.tsv",
        "reporting/GtailsNumbers_perSample.tsv"
#        report("reporting/GtailsNumbers_perLibrary.tsv",category="Reads"),
#        report("reporting/GtailsNumbers_perSample.tsv",category="Reads")
    threads: 1
    params:
        currentStep = "primers"
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    log: "logs/countPrimerReads.log"
    conda: ENVDIR + "dada2_env.yml"
    script:
        SCRIPTSDIR+"report_readNumbers.R"



if config['sequencing_direction'] == "fwd_1":
    rule cut_gtails_both:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fwd.fastq.gz",
            "preprocessing/{run}/{library}.rvs.fastq.gz"
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        params: 
            nextseq="--nextseq-trim=2" 
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Assuming forward read is read 1."
        shell:
            """
            cutadapt {params.nextseq} -j {threads}\
             -o {output[0]} -p {output[1]} {input} &> {log}
             """

elif config['sequencing_direction'] == "rvs_1":
    rule cut_gtails_both:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fwd.fastq.gz",
            "preprocessing/{run}/{library}.rvs.fastq.gz"
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        params: 
            nextseq="--nextseq-trim=2" 
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Assuming forward read is read 2."
        shell:
            """
            TMPD=$(mktemp -d -t --tmpdir={TMPDIR} "XXXXXX")
            zcat {input[0]} | fastx_reverse_complement -z -o $TMPD/{wildcards.library}.fwd.rc.fastq.gz &>> {log}\
              || touch $TMPD/{wildcards.library}.fwd.rc.fastq.gz
            zcat {input[1]} | fastx_reverse_complement -z -o $TMPD/{wildcards.library}.rvs.rc.fastq.gz &>> {log}\ 
              || touch $TMPD/{wildcards.library}.rvs.rc.fastq.gz
            cutadapt {params.nextseq} -j {threads} \
             -o {output[1]} -p {output[0]} $TMPD/{wildcards.library}.fwd.rc.fastq.gz $TMPD/{wildcards.library}.rvs.rc.fastq.gz &>> {log}
            """

else:
    rule cut_gtails_both:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fwd.fastq.gz",
            "preprocessing/{run}/{library}.rvs.fastq.gz"
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        params: 
            nextseq="--nextseq-trim=2" 
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Assuming forward read is read 1, because we just don't know. \n Note that this step does not check for the direction, so if your libraries were not sequenced with the same direction, this will not turn them."
        shell:
            """
            cutadapt {params.nextseq} -j {threads}\
             -o {output[0]} -p {output[1]} {input} &> {log}
            """



