def get_fastq(wildcards):
    return RAW+"/"+samples.loc[(wildcards.library,wildcards.run), ["r1_file"]].dropna()

def get_lib_perRunAndSample(wildcards,prefix,suffix):
    return prefix+samples.loc[(samples['run']==wildcards.run) & (samples['sample']==wildcards.sample), "library"].unique()+suffix


localrules: gtails_control

rule gtails_control:
    input:
        expand("preprocessing/{samples.run}/{samples.sample}.fastq.gz", samples=samples.itertuples()),
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
        files = lambda wildcards: get_lib_perRunAndSample(wildcards,"preprocessing/{run}/",".fastq.gz")    
    output:
        "preprocessing/{run}/{sample}.fastq.gz"
    wildcard_constraints:
        sample='|'.join(samples['sample'])
    threads: 1
    log: "logs/combine_or_rename.{run}.{sample}.log"
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
        expand("{raw_directory}/{file}", file=samples.r1_file,raw_directory=RAW)
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
        SCRIPTSDIR+"report_readNumbers.single.R" 


rule gtails_numbers:
    input:
        "reporting/readNumbers.tsv",
        expand("preprocessing/{samples.run}/{samples.library}.fastq.gz", samples=samples.itertuples())
    output:
        report("reporting/GtailsNumbers_perLibrary.tsv",category="Reads"),
        report("reporting/GtailsNumbers_perSample.tsv",category="Reads")
    threads: 1
    params:
        currentStep = "primers"
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    log: "logs/countGtailsReads.log"
    conda: ENVDIR + "dada2_env.yml"
    script:
        SCRIPTSDIR+"report_readNumbers.single.R"


#script to visualize read numbers -> run once for all steps

if config['sequencing_direction'] == "fwd_1":
    rule cut_gtails:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fastq.gz"
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        params:
            nextseq="--nextseq-trim=2"
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Assuming forward read is in read"
        shell:
            """
            cutadapt {params.nextseq} -j {threads} \
             -o {output} {input} &> {log}
            """

elif config['sequencing_direction'] == "rvs_1":
    rule cut_gtails:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fastq.gz",
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        params:
            nextseq="--nextseq-trim=2"
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Assuming read is reverse complement."
        shell:
            """
            TMPD=$(mktemp -d -t --tmpdir={TMPDIR} "XXXXXX")
            zcat {input} | fastx_reverse_complement -z -o $TMPD/{wildcards.library}.rc.fastq.gz &>> {log}\
             || touch $TMPD/{wildcards.library}.rc.fastq.gz
            cutadapt {params.nextseq} -j {threads} \
             -o {output} $TMPD/{wildcards.library}.rc.fastq.gz &> {log}
            """

else:
    rule cut_gtails_both:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fastq.gz"
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        params:
            nextseq="--nextseq-trim=2"
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Assuming the read is forward. \n Note that this step does not check for the direction, so if your libraries were not sequenced with the same direction, this will not turn them."
        shell:
            """
            cutadapt {params.nextseq} -j {threads} \
             -o {output} {input} &> {log}
            """



