def get_fastq(wildcards):
    return RAW+"/"+samples.loc[(wildcards.library,wildcards.run), ["r1_file", "r2_file"]].dropna()

def get_lib_perRunAndSample(wildcards,prefix,suffix):
    return prefix+samples.loc[(samples['run']==wildcards.run) & (samples['sample']==wildcards.sample), "library"].unique()+suffix


localrules: primers_control

rule primers_control:
    input:
        expand("preprocessing/{samples.run}/{samples.sample}.{direction}.fastq", samples=samples.itertuples(), direction=["fwd","rvs"]),
        "reporting/readNumbers.tsv",
        "reporting/primerNumbers_perSample.tsv"
    output:
        "primers.done"
    shell:
        """
        touch {output}
        """

rule combine_or_rename:
    input:
        "reporting/primerNumbers_perLibrary.tsv",
        files = lambda wildcards: get_lib_perRunAndSample(wildcards,"preprocessing/{run}/",".{direction}.fastq")    
    output:
        "preprocessing/{run}/{sample}.{direction}.fastq"
    wildcard_constraints:
        direction="(fwd|rvs)",
        sample='|'.join(samples['sample'])
    threads: 1
    log: "logs/combine_or_rename.{run}.{sample}.{direction}.log"
    params:
        runtime="01:00:00",
        mem="8G"
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
        runtime="12:00:00",
        mem="8G",
        raw_directory = RAW
    conda: "dada_env.yml"
    log: "logs/countInputReads.log"
    script:
        SRC_dir+"report_readNumbers.R" 


rule primer_numbers:
    input:
        "reporting/readNumbers.tsv",
        expand("preprocessing/{samples.run}/{samples.library}.{direction}.fastq", samples=samples.itertuples(), direction=["fwd","rvs"])
    output:
        report("reporting/primerNumbers_perLibrary.tsv",category="Reads"),
        report("reporting/primerNumbers_perSample.tsv",category="Reads")
    threads: 1
    params:
        currentStep = "primers",
        runtime="12:00:00",
        mem="8G"
    log: "logs/countPrimerReads.log"
    conda: "dada_env.yml"
    script:
        SRC_dir+"report_readNumbers.R"



if config['sequencing_direction'] == "fwd_1":
    rule cut_primer_both:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fwd.fastq",
            "preprocessing/{run}/{library}.rvs.fastq"
        threads: 1
        params:
            runtime="12:00:00",
            mem="8G"
        conda: "dada_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Assuming forward primer is in read 1. {config[primers][fwd][sequence]}"
        shell:
            """
            TMPD=$(mktemp -d -t --tmpdir={TMPDIR} 'XXXXXX') 
            FWD_RC=`echo {config[primers][fwd][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`
            RVS_RC=`echo {config[primers][rvs][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`
                
            cutadapt -g {config[primers][fwd][sequence]} -G {config[primers][rvs][sequence]} \
            --no-indels -n {config[primer_cutting][count]} -O {config[primer_cutting][overlap]} \
            -m 1:1 --pair-filter={config[primer_cutting][filter_if_not_match]} \
            -j {threads} -e {config[primer_cutting][perc_mismatch]} --trimmed-only \
            -o $TMPD/{wildcards.library}.fwd.fastq -p $TMPD/{wildcards.library}.rvs.fastq {input} &> {log}

            cutadapt -a $RVS_RC -A $FWD_RC \
             --no-indels -n {config[primer_cutting][count]} \
             -m 1:1 \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             -o {output[0]} -p {output[1]} $TMPD/{wildcards.library}.fwd.fastq $TMPD/{wildcards.library}.rvs.fastq >> {log} 2>&1
             """

elif config['sequencing_direction'] == "rvs_1":
    rule cut_primers_both:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fwd.fastq",
            "preprocessing/{run}/{library}.rvs.fastq"
        threads: 1
        params:
            runtime="12:00:00",
            mem="8G"
        conda: "dada_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Assuming forward primer is in read 2."
        shell:
            """
            TMPD=$(mktemp -d -t --tmpdir={TMPDIR} "XXXXXX")

            FWD_RC=`echo {config[primers][fwd][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`
            RVS_RC=`echo {config[primers][rvs][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`

            cutadapt -g {config[primers][rvs][sequence]} -G {config[primers][fwd][sequence]} \
             --no-indels -n {config[primer_cutting][count]} -O {config[primer_cutting][overlap]} \
              -m 1:1 --pair-filter={config[primer_cutting][filter_if_not_match]} \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} --trimmed-only \
             -o $TMPD/{wildcards.library}.rvs.fastq -p $TMPD/{wildcards.library}.fwd.fastq {input} &> {log}

            cutadapt -a $RVS_RC -A $FWD_RC \
             --no-indels -n {config[primer_cutting][count]} \
              -m 1:1 \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             -o {output[0]} -p {output[1]} $TMPD/{wildcards.library}.fwd.fastq $TMPD/{wildcards.library}.rvs.fastq >> {log}  2>&1
            """

else:
    rule cut_primers_both:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fwd.fastq",
            "preprocessing/{run}/{library}.rvs.fastq"
        threads: 1
        params:
            runtime="12:00:00",
            mem="8G"
        conda: "dada_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Searching for both  primers in both reads."
        shell:
            """
            TMPD=$(mktemp -d -t --tmpdir={TMPDIR} "XXXXXX")
            FWD_RC=`echo {config[primers][fwd][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`
            RVS_RC=`echo {config[primers][rvs][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`

            cutadapt -g {config[primers][fwd][sequence]} -G {config[primers][rvs][sequence]} \
             --no-indels -n {config[primer_cutting][count]} -O {config[primer_cutting][overlap]} \
              -m 1:1 --pair-filter={config[primer_cutting][filter_if_not_match]} \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             --untrimmed-output=$TMPD/{wildcards.library}.fwd_unt.fastq --untrimmed-paired-output=$TMPD/{wildcards.library}.rvs_unt.fastq \
             -o $TMPD/{wildcards.library}.fwd.fastq -p $TMPD/{wildcards.library}.rvs.fastq {input} &> {log}
            
            cutadapt -a $RVS_RC -A $FWD_RC \
             --no-indels -n {config[primer_cutting][count]} \
              -m 1:1 \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             -o {output[0]} -p {output[1]} $TMPD/{wildcards.library}.fwd.fastq $TMPD/{wildcards.library}.rvs.fastq >> {log} 2>&1

            cutadapt -g {config[primers][rvs][sequence]} -G {config[primers][fwd][sequence]} \
             --no-indels -n {config[primer_cutting][count]} -O {config[primer_cutting][overlap]} \
              -m 1:1 --pair-filter={config[primer_cutting][filter_if_not_match]} \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} --trimmed-only \
             -o $TMPD/{wildcards.library}.rvs_tr2.fastq -p $TMPD/{wildcards.library}.fwd_tr2.fastq \
             $TMPD/{wildcards.library}.fwd_unt.fastq $TMPD/{wildcards.library}.rvs_unt.fastq >> {log} 2>&1

            cutadapt -a $RVS_RC -A $FWD_RC \
             --no-indels -n {config[primer_cutting][count]} \
              -m 1:1 \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             -o $TMPD/{wildcards.library}.fwd.final.fastq -p $TMPD/{wildcards.library}.rvs.final.fastq \
             $TMPD/{wildcards.library}.rvs_tr2.fastq $TMPD/{wildcards.library}.fwd_tr2.fastq >> {log} 2>&1

            cat $TMPD/{wildcards.library}.fwd.final.fastq >> {output[0]}
            cat $TMPD/{wildcards.library}.rvs.final.fastq >> {output[1]}
            """



