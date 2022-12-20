taxConts = ["sequenceTables/all.seqTab.tsv","sequenceTables/all.seqTab.RDS"]
if config['taxonomy']['decipher']['do'] or config['taxonomy']['mothur']['do'] or config['taxonomy']['dada']['do']:
    taxConts.append("sequenceTables/all.seqTab.tax.tsv")
if config['blast']['do']:
    taxConts.append("sequenceTables/blast_results.tsv")
if config['hand_off']['phyloseq'] and not config['do_postprocessing']:
    taxConts.append("sequenceTables/all.seqTab.phyloseq.RDS")
if config['blast']['do'] and config['blast']['run_basta']:
    taxConts.append("sequenceTables/all.seqTab.tax.blast.RDS")

localrules: taxonomy_control

rule taxonomy_control:
    input:
        taxConts
    output:
        "taxonomy.done"
    shell:
        """
        touch {output}
        """

taxTabs = ["sequenceTables/all.seqTab.RDS"]
if config['ITSx']['do']:
    taxTabs.append("sequenceTables/ITSx.seqs.fasta")
if config['taxonomy']['decipher']['do']:
    taxTabs.append("sequenceTables/tax.decipher.RDS")
if config['taxonomy']['mothur']['do']:
    taxTabs.append("sequenceTables/tax.mothur.RDS")
if config['taxonomy']['dada']['do']:
    taxTabs.append("sequenceTables/tax.dada.RDS")

rule taxonomy_to_OTUtab:
    input:
        taxTabs
    output:
        "sequenceTables/all.seqTab.tax.tsv",
        "sequenceTables/all.seqTab.tax.RDS"
    threads: 1
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/taxonomy.log"
    message: "Combining taxa and OTU tables {input}."
    script:
        SCRIPTSDIR+"add_taxonomy.R"

rule gather_mothur_taxonomy:
    input:
        expand("sequenceTables/tax.mothur.{db}.tsv",db=[db for db in MOTHDB.keys()])
    output:
        "sequenceTables/tax.mothur.RDS"
    threads: 1
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    params:
        rank_num = config["taxonomy"]["mothur"]["rank_number"]
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/taxonomy_mothur.log"
    message: "Combining mothur taxonomy tables {input}."
    script:
        SCRIPTSDIR+"gather_mothur_taxonomy.R"

rule gather_dada_taxonomy:
    input:
        expand("sequenceTables/tax.dada.{db}.RDS",db=[db for db in DADADB.keys()])
    output:
        "sequenceTables/tax.dada.RDS"
    threads: 1
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/taxonomy_dada.log"
    message: "Combining dada taxonomy tables {input}."
    script:
        SCRIPTSDIR+"gather_dada_taxonomy.R"

rule gather_decipher_taxonomy:
    input:
        expand("sequenceTables/tax.decipher.{db}.RDS",db=[db for db in DECIDB.keys()])
    output:
        "sequenceTables/tax.decipher.RDS"
    threads: 1
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/taxonomy_decipher.log"
    message: "Combining decipher taxonomy tables {input}."
    script:
        SCRIPTSDIR+"gather_decipher_taxonomy.R"


if config['taxonomy']['dada']['post_ITSx']:
    rule dada_taxonomy:
        input:
            "sequenceTables/ITSx.seqs.fasta"
        output:
            "sequenceTables/tax.dada.{db}.tsv",
            "sequenceTables/tax.dada.{db}.RDS"
        params: 
            DB = lambda wildcards: DADADB[wildcards.db],
            SPEC = lambda wildcards: DADADB_SPEC[wildcards.db] if config['taxonomy']['dada']['look_for_species'] else ""
        threads: 1
        resources:
            runtime="48:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/dada_taxonomy_postITSx.{db}.log"
        message: "Running DADA2's classifier on {input} against {wildcards.db}."
        script:
            SCRIPTSDIR+"dadatax_ID.R"
else:
    rule dada_taxonomy:
        input:
            "sequenceTables/all.seqs.fasta"
        output:
            "sequenceTables/tax.dada.{db}.tsv",
            "sequenceTables/tax.dada.{db}.RDS"
        params: 
            DB = lambda wildcards: DADADB[wildcards.db],
            SPEC = lambda wildcards: DADADB_SPEC[wildcards.db] if config['taxonomy']['dada']['look_for_species'] else ""
        threads: 1
        resources:
            runtime="48:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/dada_taxonomy.{db}.log"
        message: "Running DADA2's classifier on {input} against {wildcards.db}."
        script:
            SCRIPTSDIR+"dadatax_ID.R"

if config['taxonomy']['decipher']['post_ITSx']:
    rule decipher_taxonomy:
        input:
            "sequenceTables/ITSx.seqs.fasta"
        output:
            "sequenceTables/tax.decipher.{db}.tsv",
            "sequenceTables/tax.decipher.{db}.RDS"
        params:
            DB = lambda wildcards: DECIDB[wildcards.db],
            SPEC = lambda wildcards: DECIDB_SPEC[wildcards.db] if config['taxonomy']['decipher']['look_for_species'] else ""
        threads: 1
        resources:
            runtime="48:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/decipher_taxonomy_postITSx.{db}.log"
        message: "Running decipher on {input} against {wildcards.db}."
        script:
            SCRIPTSDIR+"decipher_ID.R"
else:
    rule decipher_taxonomy:
        input:
            "sequenceTables/all.seqs.fasta"
        output:
            "sequenceTables/tax.decipher.{db}.tsv", 
            "sequenceTables/tax.decipher.{db}.RDS"
        params:
            DB = lambda wildcards: DECIDB[wildcards.db],
            SPEC = lambda wildcards: DECIDB_SPEC[wildcards.db] if config['taxonomy']['decipher']['look_for_species'] else ""
        threads: 1
        resources:
            runtime="48:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/decipher_taxonomy.{db}.log"
        message: "Running decipher on {input} against {wildcards.db}."
        script:
            SCRIPTSDIR+"decipher_ID.R"

if config['taxonomy']['mothur']['post_ITSx']:
    rule mothur_taxonomy_postITSx:
        input:
            "sequenceTables/ITSx.seqs.fasta"
        output:
            "sequenceTables/tax.mothur.{db}.tsv"
        threads: 1 
        resources:
            runtime="48:00:00",
            mem=config['bigMem']
        params:
            outFull= lambda wildcards: "sequenceTables/ITSx.seqs.for_" + wildcards.db,
            inBase= "ITSx.seqs.fasta",
            mothPath = lambda wildcards: MOTHPATH[wildcards.db],
            mothDb = lambda wildcards: MOTHDB[wildcards.db]
        conda: ENVDIR + "mothur_env.yml"
        log: "logs/mothur_taxonomy_postITSx.{db}.log"
        message: "Running mothur classifier on {input} against {wildcards.db}."
        shell:
            """
            ln -sfn {params.inBase} {params.outFull}.fasta
            mothur "#set.dir(tempdefault={params.mothPath});
            classify.seqs(fasta={params.outFull}.fasta, template={params.mothDb}.fasta, taxonomy={params.mothDb}.taxonomy, cutoff={config[taxonomy][mothur][cutoff]}, method=wang, processors={threads})" &> {log}
            mv {params[outFull]}.*.wang.taxonomy {output}
            rm {params.outFull}.fasta
            """
else:
    rule mothur_taxonomy:
        input:
            "sequenceTables/all.seqs.fasta"
        output:
            "sequenceTables/tax.mothur.{db}.tsv"
        threads: 1
        resources:
            runtime="48:00:00",
            mem=config['bigMem']
        params:
            outFull= lambda wildcards: "sequenceTables/all.seqs.for_" + wildcards.db,
            inBase= "all.seqs.fasta",
            mothPath = lambda wildcards: MOTHPATH[wildcards.db],
            mothDb = lambda wildcards: MOTHDB[wildcards.db]
        conda: ENVDIR + "mothur_env.yml"
        log: "logs/mothur_taxonomy.{db}.log"
        message: "Running mothur classifier on {input} against {wildcards.db}."
        shell:
            """
            ln -sfn {params.inBase} {params.outFull}.fasta
            mothur "#set.dir(tempdefault={params.mothPath});
            classify.seqs(fasta={params.outFull}.fasta, template={params.mothDb}.fasta, taxonomy={params.mothDb}.taxonomy, cutoff={config[taxonomy][mothur][cutoff]}, method=wang, processors={threads})" &> {log} 
            mv {params[outFull]}.*.wang.taxonomy {output}
            rm {params.outFull}.fasta
            """

if config['ITSx']['do']:
    rule ITSx:
        input:
            "sequenceTables/all.seqs.fasta"
        output:
            directory("sequenceTables/ITSx_out"),
            "sequenceTables/ITSx.seqs.fasta"
        threads: getThreads(12)
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        log: "logs/ITSx.log"
        conda: ENVDIR + "dadasnake_env.yml"
        message: "Running ITSx on {input}."
        shell:
            """
            PL5=${{PERL5LIB:-}}
            export PERL5LIB=$CONDA_PREFIX/lib/5.26.2/x86_64-linux-thread-multi:$PL5
            mkdir {output[0]}
            ITSx -i {input} --cpu {threads} --detailed_results T --save_regions {config[ITSx][region]} --graphical F \
             -o {output[0]}/ITSx -N {config[ITSx][min_regions]} -t {config[ITSx][query_taxa]} \
             -E {config[ITSx][e_val]} &> {log}
            grep '|{config[ITSx][target_taxon]}|{config[ITSx][region]}' -A 1 \
             --no-group-separator {output[0]}/ITSx.{config[ITSx][region]}.fasta | sed 's/|.*//' > {output[1]}
            """

if config['blast']['do']:
    if not config['blast']['all']:
        rule prepare_blastn:
            input:
                expand("sequenceTables/all.seqTab.{tax}RDS",tax="tax." if (config['taxonomy']['decipher']['do'] or config['taxonomy']['mothur']['do'] or config['taxonomy']['dada']['do']) else "")
            output:
                "sequenceTables/no_anno.seqs.fasta"
            threads: 1
            resources:
                runtime="2:00:00",
                mem=config['normalMem']
            log: "logs/prep_blastn.log"
            conda: ENVDIR + "dada2_env.yml"
            message: "Preparing blast: extracting un-annotated sequences."
            script:
                SCRIPTSDIR+"prepare_blastn.R"
        BLAST_IN="sequenceTables/no_anno.seqs.fasta"
    else:
        BLAST_IN="sequenceTables/all.seqs.fasta"
    if config['blast']['tax2id'] == "" or config['blast']['tax2id'] == "none":
        rule blastn:
            input:
                BLAST_IN
            output:
                "sequenceTables/blast_results.tsv"
            threads: 1
            resources:
                runtime="24:00:00",
                mem=config['bigMem']
            log: "logs/blastn.log"
            conda: ENVDIR + "blast_env.yml"
            message: "Running blastn on {input}."
            shell:
                """
                if [ -s {input} ]; then
                  if [ ! -f "{config[blast][db_path]}/{config[blast][tax_db]}.nin" ]
                    then
                    makeblastdb -dbtype nucl -in {config[blast][db_path]}/{config[blast][tax_db]} \
                     -out {config[blast][db_path]}/{config[blast][tax_db]} &> {log}
                  fi
                  blastn -db {config[blast][db_path]}/{config[blast][tax_db]} \
                   -query {input} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle" -out {output} -max_target_seqs {config[blast][max_targets]} &>> {log}
                else
                  touch {output}
                fi
                """

        rule basta:
            input:
                "sequenceTables/blast_results.tsv"
            output:
                "sequenceTables/basta_results.tsv",
                "sequenceTables/basta_details.txt"
            params: best=  "-b 1" if config['blast']['basta_besthit'] else ""
            threads: 1
            resources:
                runtime="24:00:00",
                mem=config['normalMem']
            log: "logs/basta.log"
            conda: ENVDIR + "basta_env.yml"
            message: "Running basta on {input}."
            shell:
                """
                {config[blast][basta_path]} sequence -e {config[blast][basta_e_val]} -l {config[blast][basta_alen]} -m {config[blast][basta_min]} -i {config[blast][basta_id]} {params.best} -p {config[blast][basta_perchits]} -d {config[blast][basta_db]} -v {output[1]} -b {config[blast][basta_besthit]} {input} {output[0]} gb &>> {log}
                """

        rule format_basta:
            input:
                "sequenceTables/basta_results.tsv",
                "sequenceTables/all.seqTab.tax.RDS"
            output:
                "sequenceTables/all.seqTab.tax.blast.RDS",
                "sequenceTables/all.seqTab.tax.blast.tsv"
            threads: 1
            resources:
                runtime="24:00:00",
                mem=config['bigMem']
            log: "logs/basta_format.log"
            conda: ENVDIR + "dada2_env.yml"
            message: "Formatting basta output {input}."
            script:
                SCRIPTSDIR+"format_basta.R"
    else:
        rule blastn:
            input:
                BLAST_IN
            output:
                "sequenceTables/blast_results.tsv"
            threads: 1
            resources:
                runtime="24:00:00",
                mem=config['bigMem']
            log: "logs/blastn.log"
            conda: ENVDIR + "blast_env.yml"
            message: "Running blastn on {input}."
            shell:
                """
                if [ -s {input} ]; then
                  export BLASTDB={config[blast][db_path]}
                  blastn -query {input} -db {config[blast][db_path]}/{config[blast][tax_db]} \
                   -out sequenceTables/blast_output.{config[blast][tax_db]}.tsv \
                   -evalue {config[blast][e_val]} -max_target_seqs {config[blast][max_targets]} \
                   -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore sgi sacc staxids ssciname scomnames stitle' &> {log}
                 awk -F $"\\t" '{{if (NR==FNR) {{val[$13] = $0; next}} if($1 in val){{print val[$1]"\\t"$0}}}}' \
                  sequenceTables/blast_output.{config[blast][tax_db]}.tsv {config[blast][tax2id]} >> {output}
                else
                  touch {output}
                fi  
                """

if config['hand_off']['biom'] and (config['taxonomy']['decipher']['do'] or config['taxonomy']['mothur']['do'] or config['taxonomy']['dada']['do']):
    rule biom_handoff_tax:
        input:
            "sequenceTables/all.seqTab.tax.RDS",
            "reporting/finalNumbers_perSample.tsv"
        output:
            "sequenceTables/all.seqTab.biom"
        threads: 1
        params:
            currentStep = "taxonomy"
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "add_R_env.yml"
        log: "logs/biom_hand-off.log"
        script:
            SCRIPTSDIR+"biom_handoff.R"

if config['hand_off']['phyloseq'] and (config['taxonomy']['decipher']['do'] or config['taxonomy']['mothur']['do'] or config['taxonomy']['dada']['do']) and not config['do_postprocessing']:
    rule phyloseq_handoff_tax:
        input:
            "sequenceTables/all.seqTab.tax.RDS",
            "reporting/finalNumbers_perSample.tsv"
        output:
            "sequenceTables/all.seqTab.phyloseq.RDS"
        threads: 1
        params:
            currentStep = "taxonomy"
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "add_R_env.yml"
        log: "logs/phyloseq_hand-off.log"
        script:
            SCRIPTSDIR+"phyloseq_handoff.R"



