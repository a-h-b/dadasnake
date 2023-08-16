def get_fastq(wildcards):
    return RAW+"/"+samples.loc[(wildcards.library,wildcards.run), ["r1_file"]].dropna()

def get_lib_perRunAndSample(wildcards,prefix,suffix):
    return prefix+samples.loc[(samples['run']==wildcards.run) & (samples['sample']==wildcards.sample), "library"].unique()+suffix


localrules: primers_control

rule primers_control:
    input:
        expand("preprocessing/{samples.run}/{samples.sample}.fastq.gz", samples=samples.itertuples()),
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


rule primer_numbers:
    input:
        "reporting/readNumbers.tsv",
        expand("preprocessing/{samples.run}/{samples.library}.fastq.gz", samples=samples.itertuples())
    output:
        report("reporting/primerNumbers_perLibrary.tsv",category="Reads"),
        report("reporting/primerNumbers_perSample.tsv",category="Reads")
    threads: 1
    params:
        currentStep = "primers"
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    log: "logs/countPrimerReads.log"
    conda: ENVDIR + "dada2_env.yml"
    script:
        SCRIPTSDIR+"report_readNumbers.single.R"


#script to visualize read numbers -> run once for all steps

if config['primer_cutting']['both_primers_in_read']:
    BOTHMATCH = "--trimmed-only"
else:
    BOTHMATCH = ""

if config['sequencing_direction'] == "fwd_1":
    rule cut_primer_both:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fastq.gz"
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        params:
            both_match=BOTHMATCH,
            nextseq="--nextseq-trim=2" if config['nextseq_novaseq'] else ""
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Assuming forward primer is in read {config[primers][fwd][sequence]}"
        shell:
            """
            TMPD=$(mktemp -d -t --tmpdir={TMPDIR} 'XXXXXX') 
            FWD_RC=`echo {config[primers][fwd][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`
            RVS_RC=`echo {config[primers][rvs][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`
            cutadapt -g {config[primers][fwd][sequence]} \
            {config[primer_cutting][indels]} -n {config[primer_cutting][count]} -O {config[primer_cutting][overlap]} \
            -m 1 {params.nextseq}\
            -j {threads} -e {config[primer_cutting][perc_mismatch]} --trimmed-only \
            -o $TMPD/{wildcards.library}.fastq.gz {input} &> {log}

            cutadapt -a $RVS_RC {config[primer_cutting][indels]} \
             -n {config[primer_cutting][count]} -m 1 \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             {params.both_match} \
             -o {output} $TMPD/{wildcards.library}.fastq.gz >> {log} 2>&1
             """

elif config['sequencing_direction'] == "rvs_1":
    rule cut_primers_both:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fastq.gz",
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        params:
            both_match=BOTHMATCH,
            nextseq="--nextseq-trim=2" if config['nextseq_novaseq'] else ""
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Assuming reverse primer is in read."
        shell:
            """
            TMPD=$(mktemp -d -t --tmpdir={TMPDIR} "XXXXXX")

            FWD_RC=`echo {config[primers][fwd][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`
            RVS_RC=`echo {config[primers][rvs][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`

            cutadapt -g {config[primers][rvs][sequence]} \
            {config[primer_cutting][indels]} -n {config[primer_cutting][count]} -O {config[primer_cutting][overlap]} \
              -m 1 {params.nextseq}\
             -j {threads} -e {config[primer_cutting][perc_mismatch]} --trimmed-only \
             -o $TMPD/{wildcards.library}.fastq.gz {input} &> {log}

            cutadapt -a $FWD_RC \
            {config[primer_cutting][indels]} -n {config[primer_cutting][count]} \
              -m 1 \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             {params.both_match} \
             -o {output} $TMPD/{wildcards.library}.fastq.gz >> {log}  2>&1
            """

else:
    rule cut_primers_both:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fastq.gz"
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        params:
            both_match=BOTHMATCH,
            nextseq="--nextseq-trim=2" if config['nextseq_novaseq'] else ""
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Searching for both primers.\n Note that this step does not check for the direction, so if your libraries were not sequenced with the same direction, this will not turn them."
        shell:
            """
            TMPD=$(mktemp -d -t --tmpdir={TMPDIR} "XXXXXX")
            FWD_RC=`echo {config[primers][fwd][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`
            RVS_RC=`echo {config[primers][rvs][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`

            cutadapt -g {config[primers][fwd][sequence]} \
            {config[primer_cutting][indels]} -n {config[primer_cutting][count]} -O {config[primer_cutting][overlap]} \
              -m 1 {params.nextseq}\
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             --untrimmed-output=$TMPD/{wildcards.library}.unt.fastq.gz \
             -o $TMPD/{wildcards.library}.fastq.gz {input} &> {log}
            
            cutadapt -a $RVS_RC \
            {config[primer_cutting][indels]} -n {config[primer_cutting][count]} \
              -m 1 \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             {params.both_match} \
             -o {output} $TMPD/{wildcards.library}.fastq.gz >> {log} 2>&1

            cutadapt -g {config[primers][rvs][sequence]} \
            {config[primer_cutting][indels]} -n {config[primer_cutting][count]} -O {config[primer_cutting][overlap]} \
              -m 1\
             -j {threads} -e {config[primer_cutting][perc_mismatch]} --trimmed-only \
             -o $TMPD/{wildcards.library}.tr2.fastq.gz \
             $TMPD/{wildcards.library}.unt.fastq.gz >> {log} 2>&1

            zcat $TMPD/{wildcards.library}.tr2.fastq.gz | fastx_reverse_complement -z -o $TMPD/{wildcards.library}.tr2.rc.fastq.gz &>> {log} || touch $TMPD/{wildcards.library}.tr2.rc.fastq.gz

            cutadapt -a {config[primers][fwd][sequence]} \
            {config[primer_cutting][indels]} -n {config[primer_cutting][count]} \
              -m 1 \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             -o $TMPD/{wildcards.library}.final.fastq.gz \
             {params.both_match} \
             $TMPD/{wildcards.library}.tr2.rc.fastq.gz  >> {log} 2>&1

            cat $TMPD/{wildcards.library}.final.fastq.gz >> {output}
            """



