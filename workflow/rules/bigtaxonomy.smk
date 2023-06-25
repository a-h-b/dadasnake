taxConts = ["sequenceTables/all.seqTab.tsv","sequenceTables/all.seqTab.RDS"]
if config['taxonomy']['decipher']['do'] and 'ASV' in config['taxonomy']['decipher']['run_on']: 
    taxConts.append("sequenceTables/all.seqTab.tax.tsv")
elif config['taxonomy']['mothur']['do'] and 'ASV' in config['taxonomy']['mothur']['run_on']:
    taxConts.append("sequenceTables/all.seqTab.tax.tsv")
elif config['taxonomy']['dada']['do'] and 'ASV' in config['taxonomy']['dada']['run_on']:
    taxConts.append("sequenceTables/all.seqTab.tax.tsv")
if config['blast']['do'] and 'ASV' in config['blast']['run_on']:
    taxConts.append("sequenceTables/blast_results.tsv")
if config['hand_off']['phyloseq'] and "sequenceTables/all.seqTab.tax.tsv" in taxConts and not config['do_postprocessing']:
    taxConts.append("sequenceTables/all.seqTab.phyloseq.RDS")
if config['blast']['do'] and 'ASV' in config['blast']['run_on'] and config['blast']['run_basta']:
    taxConts.append("sequenceTables/all.seqTab.tax.blast.RDS")
if config['taxonomy']['decipher']['do'] and 'cluster' in config['taxonomy']['decipher']['run_on']:
    taxConts.append("clusteredTables/clusteredTab.tax.tsv")
elif config['taxonomy']['mothur']['do'] and 'cluster' in config['taxonomy']['mothur']['run_on']:
    taxConts.append("clusteredTables/clusteredTab.tax.tsv")
elif config['taxonomy']['dada']['do'] and 'cluster' in config['taxonomy']['dada']['run_on']:
    taxConts.append("clusteredTables/clusteredTab.tax.tsv")
if config['blast']['do'] and 'cluster' in config['blast']['run_on']:
    taxConts.append("clusteredTables/blast_results.tsv")
if config['hand_off']['phyloseq'] and "clusteredTables/clusteredTab.tax.tsv" in taxConts and not config['do_postprocessing']:
    taxConts.append("clusteredTables/clusteredTab.phyloseq.RDS")
if config['blast']['do'] and 'cluster' in config['blast']['run_on'] and config['blast']['run_basta']:
    taxConts.append("clusteredTables/clusteredTab.tax.blast.RDS")
if config['ITSx']['do'] and 'ASV' in config['ITSx']['run_on']:
    taxConts.append("sequenceTables/ITSx.seqs.fasta")
if config['ITSx']['do'] and 'cluster' in config['ITSx']['run_on']:
    taxConts.append("clusteredTables/ITSx.seqs.fasta")


localrules: taxonomy_control

rule taxonomy_control:
    input:
        taxConts
    output:
        touch("taxonomy.done")

taxTabs = ["sequenceTables/all.seqTab.RDS"]
if config['ITSx']['do'] and 'ASV' in config['ITSx']['run_on']:
    taxTabs.append("sequenceTables/ITSx.seqs.fasta")
if config['taxonomy']['decipher']['do'] and 'ASV' in config['taxonomy']['decipher']['run_on']:
    taxTabs.append("sequenceTables/tax.decipher.RDS")
if config['taxonomy']['mothur']['do'] and 'ASV' in config['taxonomy']['mothur']['run_on']:
    taxTabs.append("sequenceTables/tax.mothur.RDS")
if config['taxonomy']['dada']['do'] and 'ASV' in config['taxonomy']['dada']['run_on']:
    taxTabs.append("sequenceTables/tax.dada.RDS")

taxTabsCl = ["clusteredTables/clusteredTab.RDS"]
if config['ITSx']['do'] and 'cluster' in config['ITSx']['run_on']:
    taxTabsCl.append("clusteredTables/ITSx.seqs.fasta")
if config['taxonomy']['decipher']['do'] and 'cluster' in config['taxonomy']['decipher']['run_on']:
    taxTabsCl.append("clusteredTables/tax.decipher.RDS")
if config['taxonomy']['mothur']['do'] and 'cluster' in config['taxonomy']['mothur']['run_on']:
    taxTabsCl.append("clusteredTables/tax.mothur.RDS")
if config['taxonomy']['dada']['do'] and 'cluster' in config['taxonomy']['dada']['run_on']:
    taxTabsCl.append("clusteredTables/tax.dada.RDS")


rule taxonomy_to_ASVtab:
    input:
        taxTabs
    output:
        "sequenceTables/all.seqTab.tax.tsv",
        "sequenceTables/all.seqTab.tax.RDS"
    params:
        what="ASV",
        ITSx="add" if config['ITSx']['do'] and 'ASV' in config['ITSx']['run_on'] else "none"
    threads: config['bigCores']
    resources:
        runtime="12:00:00",
        mem=config['bigMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/taxonomy_toASVtable.log"
    message: "Combining taxa and OTU tables {input}."
    script:
        SCRIPTSDIR+"add_taxonomy.R"

rule taxonomy_to_clustertab:
    input:
        taxTabsCl
    output:
        "clusteredTables/clusteredTab.tax.tsv",
        "clusteredTables/clusteredTab.tax.RDS"
    params:
        what="cluster",
        ITSx="add" if config['ITSx']['do'] and 'cluster' in config['ITSx']['run_on'] else "none"
    threads: config['bigCores']
    resources:
        runtime="12:00:00",
        mem=config['bigMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/taxonomy_toClusterTable.log"
    message: "Combining taxa and clustered tables {input}."
    script:
        SCRIPTSDIR+"add_taxonomy.R"


rule gather_mothur_taxonomy:
    input:
        expand("sequenceTables/tax.mothur.{db}.tsv",db=[db for db in MOTHDB.keys()])
    output:
        "sequenceTables/tax.mothur.RDS"
    threads: 1
    resources:
        runtime="120:00:00",
        mem=config['bigMem']
    params:
        rank_num = config["taxonomy"]["mothur"]["rank_number"],
        what="ASV"
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
    params:
        what="ASV"
    threads: config['bigCores']
    resources:
        runtime="12:00:00",
        mem=config['bigMem']
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
    params:
        what="ASV"
    threads: config['bigCores']
    resources:
        runtime="12:00:00",
        mem=config['bigMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/taxonomy_decipher.log"
    message: "Combining decipher taxonomy tables {input}."
    script:
        SCRIPTSDIR+"gather_decipher_taxonomy.R"


rule gather_mothur_taxonomyCl:
    input:
        expand("clusteredTables/tax.mothur.{db}.tsv",db=[db for db in MOTHDB.keys()])
    output:
        "clusteredTables/tax.mothur.RDS"
    threads: config['bigCores']
    resources:
        runtime="12:00:00",
        mem=config['bigMem']
    params:
        rank_num = config["taxonomy"]["mothur"]["rank_number"],
        what="cluster"
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/taxonomy_mothur_clusters.log"
    message: "Combining mothur taxonomy tables {input}."
    script:
        SCRIPTSDIR+"gather_mothur_taxonomy.R"

rule gather_dada_taxonomyCl:
    input:
        expand("clusteredTables/tax.dada.{db}.RDS",db=[db for db in DADADB.keys()])
    output:
        "clusteredTables/tax.dada.RDS"
    params:
        what="cluster"
    threads: config['bigCores']
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/taxonomy_dada_clusters.log"
    message: "Combining dada taxonomy tables {input}."
    script:
        SCRIPTSDIR+"gather_dada_taxonomy.R"

rule gather_decipher_taxonomyCl:
    input:
        expand("clusteredTables/tax.decipher.{db}.RDS",db=[db for db in DECIDB.keys()])
    output:
        "clusteredTables/tax.decipher.RDS"
    params:
        what="cluster"
    threads: config['bigCores']
    resources:
        runtime="12:00:00",
        mem=config['bigMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/taxonomy_decipher_clusters.log"
    message: "Combining decipher taxonomy tables {input}."
    script:
        SCRIPTSDIR+"gather_decipher_taxonomy.R"

if config['taxonomy']['dada']['post_ITSx'] and config['ITSx']['do']:
    rule dada_taxonomy:
        input:
            "sequenceTables/ITSx.seqs.fasta"
        output:
            "sequenceTables/tax.dada.{db}.tsv",
            "sequenceTables/tax.dada.{db}.RDS"
        params: 
            DB = lambda wildcards: DADADB[wildcards.db],
            SPEC = lambda wildcards: DADADB_SPEC[wildcards.db] if config['taxonomy']['dada']['look_for_species'] else "",
            do_species=config['taxonomy']['dada']['look_for_species'],
            what="ASV",
            seed=config['taxonomy']['dada']['seed'],
            db_path=config['taxonomy']['dada']['db_path'],
            ref=config['taxonomy']['dada']['refFasta'],
            tryRC=config['taxonomy']['dada']['tryRC'],
            minBoot=config['taxonomy']['dada']['minBoot']
        threads: config['bigCores']
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
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
            SPEC = lambda wildcards: DADADB_SPEC[wildcards.db] if config['taxonomy']['dada']['look_for_species'] else "",
            do_species=config['taxonomy']['dada']['look_for_species'],
            what="ASV",
            seed=config['taxonomy']['dada']['seed'],
            db_path=config['taxonomy']['dada']['db_path'],
            ref=config['taxonomy']['dada']['refFasta'],
            tryRC=config['taxonomy']['dada']['tryRC'],
            minBoot=config['taxonomy']['dada']['minBoot']
        threads: config['bigCores']
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/dada_taxonomy.{db}.log"
        message: "Running DADA2's classifier on {input} against {wildcards.db}."
        script:
            SCRIPTSDIR+"dadatax_ID.R"

if config['taxonomy']['decipher']['post_ITSx'] and config['ITSx']['do']:
    rule decipher_taxonomy:
        input:
            "sequenceTables/ITSx.seqs.fasta"
        output:
            "sequenceTables/tax.decipher.{db}.tsv",
            "sequenceTables/tax.decipher.{db}.RDS"
        params:
            DB = lambda wildcards: DECIDB[wildcards.db],
            SPEC = lambda wildcards: DECIDB_SPEC[wildcards.db] if config['taxonomy']['decipher']['look_for_species'] else "",
            what="ASV",
            db_path=config['taxonomy']['decipher']['db_path'],
            tax_db=config['taxonomy']['decipher']['tax_db'],
            seed=config['taxonomy']['decipher']['seed'],
            strand=config['taxonomy']['decipher']['strand'],
            threshold=config['taxonomy']['decipher']['threshold'],
            bootstraps=config['taxonomy']['decipher']['bootstraps'],
            do_species=config['taxonomy']['decipher']['look_for_species']
        threads: config['bigCores']
        resources:
            runtime="120:00:00",
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
            SPEC = lambda wildcards: DECIDB_SPEC[wildcards.db] if config['taxonomy']['decipher']['look_for_species'] else "",
            what="ASV",
            db_path=config['taxonomy']['decipher']['db_path'],
            tax_db=config['taxonomy']['decipher']['tax_db'],
            seed=config['taxonomy']['decipher']['seed'],
            strand=config['taxonomy']['decipher']['strand'],
            threshold=config['taxonomy']['decipher']['threshold'],
            bootstraps=config['taxonomy']['decipher']['bootstraps'],
            do_species=config['taxonomy']['decipher']['look_for_species']
        threads: config['bigCores']
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/decipher_taxonomy.{db}.log"
        message: "Running decipher on {input} against {wildcards.db}."
        script:
            SCRIPTSDIR+"decipher_ID.R"

if config['taxonomy']['mothur']['post_ITSx'] and config['ITSx']['do']:
    rule mothur_taxonomy_postITSx:
        input:
            "sequenceTables/ITSx.seqs.fasta"
        output:
            "sequenceTables/tax.mothur.{db}.tsv"
        threads: 1 
        resources:
            runtime="120:00:00",
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
            runtime="120:00:00",
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

if config['taxonomy']['dada']['post_ITSx'] and config['ITSx']['do']:
    rule dada_taxonomyCl:
        input:
            "clusteredTables/ITSx.seqs.fasta"
        output:
            "clusteredTables/tax.dada.{db}.tsv",
            "clusteredTables/tax.dada.{db}.RDS"
        params:
            DB = lambda wildcards: DADADB[wildcards.db],
            SPEC = lambda wildcards: DADADB_SPEC[wildcards.db] if config['taxonomy']['dada']['look_for_species'] else "",
            do_species=config['taxonomy']['dada']['look_for_species'],
            what="cluster",
            seed=config['taxonomy']['dada']['seed'],
            db_path=config['taxonomy']['dada']['db_path'],
            ref=config['taxonomy']['dada']['refFasta'],
            tryRC=config['taxonomy']['dada']['tryRC'],
            minBoot=config['taxonomy']['dada']['minBoot']
        threads: config['bigCores']
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/dada_taxonomy_postITSx_cluster.{db}.log"
        message: "Running DADA2's classifier on {input} against {wildcards.db}."
        script:
            SCRIPTSDIR+"dadatax_ID.R"
else:
    rule dada_taxonomyCl:
        input:
            "clusteredTables/consensus.fasta"
        output:
            "clusteredTables/tax.dada.{db}.tsv",
            "clusteredTables/tax.dada.{db}.RDS"
        params:
            DB = lambda wildcards: DADADB[wildcards.db],
            SPEC = lambda wildcards: DADADB_SPEC[wildcards.db] if config['taxonomy']['dada']['look_for_species'] else "",
            do_species=config['taxonomy']['dada']['look_for_species'],
            what="cluster",
            seed=config['taxonomy']['dada']['seed'],
            db_path=config['taxonomy']['dada']['db_path'],
            ref=config['taxonomy']['dada']['refFasta'],
            tryRC=config['taxonomy']['dada']['tryRC'],
            minBoot=config['taxonomy']['dada']['minBoot']
        threads: config['bigCores']
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/dada_taxonomy_cluster.{db}.log"
        message: "Running DADA2's classifier on {input} against {wildcards.db}."
        script:
            SCRIPTSDIR+"dadatax_ID.R"


if config['taxonomy']['decipher']['post_ITSx'] and config['ITSx']['do']:
    rule decipher_taxonomyCl:
        input:
            "clusteredTables/ITSx.seqs.fasta"
        output:
            "clusteredTables/tax.decipher.{db}.tsv",
            "clursteredTables/tax.decipher.{db}.RDS"
        params:
            DB = lambda wildcards: DECIDB[wildcards.db],
            SPEC = lambda wildcards: DECIDB_SPEC[wildcards.db] if config['taxonomy']['decipher']['look_for_species'] else "",
            what="cluster",
            db_path=config['taxonomy']['decipher']['db_path'],
            tax_db=config['taxonomy']['decipher']['tax_db'],
            seed=config['taxonomy']['decipher']['seed'],
            strand=config['taxonomy']['decipher']['strand'],
            threshold=config['taxonomy']['decipher']['threshold'],
            bootstraps=config['taxonomy']['decipher']['bootstraps'],
            do_species=config['taxonomy']['decipher']['look_for_species']
        threads: config['bigCores']
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/decipher_taxonomy_postITSx_cluster.{db}.log"
        message: "Running decipher on {input} against {wildcards.db}."
        script:
            SCRIPTSDIR+"decipher_ID.R"
else:
    rule decipher_taxonomyCl:
        input:
            "clusteredTables/consensus.fasta"
        output:
            "clusteredTables/tax.decipher.{db}.tsv",
            "clusteredTables/tax.decipher.{db}.RDS"
        params:
            DB = lambda wildcards: DECIDB[wildcards.db],
            SPEC = lambda wildcards: DECIDB_SPEC[wildcards.db] if config['taxonomy']['decipher']['look_for_species'] else "",
            what="cluster",
            db_path=config['taxonomy']['decipher']['db_path'],
            tax_db=config['taxonomy']['decipher']['tax_db'],
            seed=config['taxonomy']['decipher']['seed'],
            strand=config['taxonomy']['decipher']['strand'],
            threshold=config['taxonomy']['decipher']['threshold'],
            bootstraps=config['taxonomy']['decipher']['bootstraps'],
            do_species=config['taxonomy']['decipher']['look_for_species']
        threads: config['bigCores']
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "dada2_env.yml"
        log: "logs/decipher_taxonomy_cluster.{db}.log"
        message: "Running decipher on {input} against {wildcards.db}."
        script:
            SCRIPTSDIR+"decipher_ID.R"


if config['taxonomy']['mothur']['post_ITSx'] and config['ITSx']['do']:
    rule mothur_taxonomy_postITSxCl:
        input:
            "clusteredTables/ITSx.seqs.fasta"
        output:
            "clusteredTables/tax.mothur.{db}.tsv"
        threads: 1
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        params:
            outFull= lambda wildcards: "clusteredTables/ITSx.seqs.for_" + wildcards.db,
            inBase= "ITSx.seqs.fasta",
            mothPath = lambda wildcards: MOTHPATH[wildcards.db],
            mothDb = lambda wildcards: MOTHDB[wildcards.db]
        conda: ENVDIR + "mothur_env.yml"
        log: "logs/mothur_taxonomy_postITSx_cluster.{db}.log"
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
    rule mothur_taxonomyCl:
        input:
            "clusteredTables/consensus.fasta"
        output:
            "clusteredTables/tax.mothur.{db}.tsv"
        threads: 1
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        params:
            outFull= lambda wildcards: "clusteredTables/consensus.for_" + wildcards.db,
            inBase= "consensus.fasta",
            mothPath = lambda wildcards: MOTHPATH[wildcards.db],
            mothDb = lambda wildcards: MOTHDB[wildcards.db]
        conda: ENVDIR + "mothur_env.yml"
        log: "logs/mothur_taxonomy_cluster.{db}.log"
        message: "Running mothur classifier on {input} against {wildcards.db}."
        shell:
            """
            ln -sfn {params.inBase} {params.outFull}.fasta
            mothur "#set.dir(tempdefault={params.mothPath});
            classify.seqs(fasta={params.outFull}.fasta, template={params.mothDb}.fasta, taxonomy={params.mothDb}.taxonomy, cutoff={config[taxonomy][mothur][cutoff]}, method=wang, processors={threads})" &> {log}
            mv {params[outFull]}.*.wang.taxonomy {output}
            rm {params.outFull}.fasta
            """



rule ITSx:
    input:
        "sequenceTables/all.seqs.fasta"
    output:
        directory("sequenceTables/ITSx_out"),
        "sequenceTables/ITSx.seqs.fasta"
    threads: config['bigCores']
    resources:
        runtime="120:00:00",
        mem=config['bigMem']
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
        if [[ -s {output[0]}/ITSx.{config[ITSx][region]}.fasta ]]
        then
           grep '|{config[ITSx][target_taxon]}|{config[ITSx][region]}' -A 1 \
            --no-group-separator {output[0]}/ITSx.{config[ITSx][region]}.fasta \
            | sed 's/|.*//' > {output[1]} 2>> {log}
        else
           echo "no ITS found" >> /tmp/sizelog.log
           touch {output[1]}
        fi
        """

rule ITSxCl:
    input:
        "clusteredTables/consensus.fasta"
    output:
        directory("clusteredTables/ITSx_out"),
        "clusteredTables/ITSx.seqs.fasta"
    threads: config['bigCores']
    resources:
        runtime="120:00:00",
        mem=config['bigMem']
    log: "logs/ITSx_cluster.log"
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

        if [[ -s {output[0]}/ITSx.{config[ITSx][region]}.fasta ]]
        then
           grep '|{config[ITSx][target_taxon]}|{config[ITSx][region]}' -A 1 \
            --no-group-separator {output[0]}/ITSx.{config[ITSx][region]}.fasta \
            | sed 's/|.*//' > {output[1]} 2>> {log}
        else
           echo "no ITS found" >> /tmp/sizelog.log
           touch {output[1]}
        fi
        """

if not config['blast']['all']:
    rule prepare_blastn:
        input:
            expand("sequenceTables/all.seqTab.{tax}RDS",tax="tax." if (config['taxonomy']['decipher']['do'] or config['taxonomy']['mothur']['do'] or config['taxonomy']['dada']['do']) else "")
        output:
            "sequenceTables/no_anno.seqs.fasta"
        params:
            what="ASV"
        threads: config['bigCores']
        resources:
            runtime="12:00:00",
            mem=config['bigMem']
        log: "logs/prep_blastn.log"
        conda: ENVDIR + "dada2_env.yml"
        message: "Preparing blast: extracting un-annotated sequences."
        script:
            SCRIPTSDIR+"prepare_blastn.R"

    BLAST_IN="sequenceTables/no_anno.seqs.fasta"

    rule prepare_blastnCl:
        input:
            expand("clusteredTables/clusteredTab.{tax}RDS",tax="tax." if (config['taxonomy']['decipher']['do'] or config['taxonomy']['mothur']['do'] or config['taxonomy']['dada']['do']) else "")
        output:
            "clusteredTables/no_anno.seqs.fasta"
        params:
            what="cluster"
        threads: config['bigCores']
        resources:
            runtime="12:00:00",
            mem=config['bigMem']
        log: "logs/prep_blastn_cluster.log"
        conda: ENVDIR + "dada2_env.yml"
        message: "Preparing blast: extracting un-annotated sequences."
        script:
            SCRIPTSDIR+"prepare_blastn.R"

    BLAST_IN_CL="clusteredTables/no_anno.seqs.fasta"

else:
    BLAST_IN="sequenceTables/all.seqs.fasta"
    BLAST_IN_CL="clusteredTables/consensus.fasta"
if config['blast']['tax2id'] == "" or config['blast']['tax2id'] == "none":
    rule blastn:
        input:
            BLAST_IN
        output:
            "sequenceTables/blast_results.tsv"
        threads: 1
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        log: "logs/blastn.log"
        conda: ENVDIR + "blast_env.yml"
        message: "Running blastn on {input}."
        shell:
            """
            if [ -s {input} ]; then
              if [ ! -f "{config[blast][db_path]}/{config[blast][tax_db]}.nin" -a ! -f "{config[blast][db_path]}/{config[blast][tax_db]}.nal" ]
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
            "{dir}/blast_results.tsv"
        output:
            "{dir}/basta_results.tsv",
            "{dir}/basta_details.txt"
        params:
            best=  "-b 1" if config['blast']['basta_besthit'] else "",
            e_val=config['blast']['basta_e_val'],
            alen=config['blast']['basta_alen'],
            min=config['blast']['basta_min'],
            id=config['blast']['basta_id'],
            perchits=config['blast']['basta_perchits'],
            besthit=config['blast']['basta_besthit'],
            db_path=config['blast']['basta_db']
        threads: 1
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        log: "logs/basta.{dir}.log"
        conda: ENVDIR + "basta_env.yml"
        message: "Running basta on {input}."
        shell:
            """
            TMPD=$(mktemp -d -t --tmpdir={TMPDIR} "XXXXXX")
            cp -r {params.db_path}/gb_mapping.db $TMPD/
            cp -r {params.db_path}/complete_taxa.db $TMPD/
            basta sequence -e {params.e_val} -l {params.alen} \
              -m {params.min} -i {params.id} {params.best} \
              -p {params.perchits} -d $TMPD -v {output[1]} \
              {input} {output[0]} gb &>> {log}
            rm -rf $TMPD/gb_mapping.db
            rm -rf $TMPD/complete_taxa.db            
            """

    rule format_basta:
        input:
            "sequenceTables/basta_results.tsv",
            "sequenceTables/all.seqTab.tax.RDS"
        output:
            "sequenceTables/all.seqTab.tax.blast.RDS",
            "sequenceTables/all.seqTab.tax.blast.tsv"
        params:
            what="ASV",
            tax_db=config["blast"]["tax_db"]
        threads: 1
        resources:
            runtime="24:00:00",
            mem=config['bigMem']
        log: "logs/basta_format.log"
        conda: ENVDIR + "dada2_env.yml"
        message: "Formatting basta output {input}."
        script:
            SCRIPTSDIR+"format_basta.R"

    rule blastnCl:
        input:
            BLAST_IN_CL
        output:
            "clusteredTables/blast_results.tsv"
        threads: 1
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        log: "logs/blastn_cluster.log"
        conda: ENVDIR + "blast_env.yml"
        message: "Running blastn on {input}."
        shell:
            """
            if [ -s {input} ]; then
              if [ ! -f "{config[blast][db_path]}/{config[blast][tax_db]}.nin" -a ! -f "{config[blast][db_path]}/{config[blast][tax_db]}.nal" ]
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

    rule format_bastaCl:
        input:
            "clusteredTables/basta_results.tsv",
            "clusteredTables/clusteredTab.tax.RDS"
        output:
            "clusteredTables/clusteredTab.tax.blast.RDS",
            "clusteredTables/clusteredTab.tax.blast.tsv"
        params:
            what="cluster",
            tax_db=config["blast"]["tax_db"]
        threads: 1
        resources:
            runtime="24:00:00",
            mem=config['bigMem']
        log: "logs/basta_format_cluster.log"
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
        params:
            outdir = "sequenceTables"
        threads: 1
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        log: "logs/blastn.log"
        conda: ENVDIR + "blast_env.yml"
        message: "Running blastn on {input}."
        shell:
            """
            if [ -s {input} ]; then
              export BLASTDB={config[blast][db_path]}
              blastn -query {input} -db {config[blast][db_path]}/{config[blast][tax_db]} \
               -out {params.outdir}/blast_output.{config[blast][tax_db]}.tsv \
               -evalue {config[blast][e_val]} -max_target_seqs {config[blast][max_targets]} \
               -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore sgi sacc staxids ssciname scomnames stitle' &> {log}
             awk -F $"\\t" '{{if (NR==FNR) {{val[$13] = $0; next}} if($1 in val){{print val[$1]"\\t"$0}}}}' \
              {params.outdir}/blast_output.{config[blast][tax_db]}.tsv {config[blast][tax2id]} >> {output}
            else
              touch {output}
            fi  
            """

    rule blastnCl:
        input:
            BLAST_IN_CL
        output:
            "clusteredTables/blast_results.tsv"
        params:
            outdir = "clusteredTables"
        threads: 1
        resources:
            runtime="120:00:00",
            mem=config['bigMem']
        log: "logs/blastn_cluster.log"
        conda: ENVDIR + "blast_env.yml"
        message: "Running blastn on {input}."
        shell:
            """
            if [ -s {input} ]; then
              export BLASTDB={config[blast][db_path]}
              blastn -query {input} -db {config[blast][db_path]}/{config[blast][tax_db]} \
               -out {params.outdir}/blast_output.{config[blast][tax_db]}.tsv \
               -evalue {config[blast][e_val]} -max_target_seqs {config[blast][max_targets]} \
               -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore sgi sacc staxids ssciname scomnames stitle' &> {log}
             awk -F $"\\t" '{{if (NR==FNR) {{val[$13] = $0; next}} if($1 in val){{print val[$1]"\\t"$0}}}}' \
              {params.outdir}/blast_output.{config[blast][tax_db]}.tsv {config[blast][tax2id]} >> {output}
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
            runtime="24:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "add_R_env.yml"
        log: "logs/biom_hand-off.log"
        script:
            SCRIPTSDIR+"biom_handoff.R"


if config['hand_off']['phyloseq'] and not config['do_postprocessing']:
    if config['taxonomy']['decipher']['do'] and 'ASV' in config['taxonomy']['decipher']['run_on']:
        physInputs = ["sequenceTables/all.seqTab.tax.RDS","reporting/finalNumbers_perSample.tsv"]
    elif config['taxonomy']['mothur']['do'] and 'ASV' in config['taxonomy']['mothur']['run_on']:
        physInputs = ["sequenceTables/all.seqTab.tax.RDS","reporting/finalNumbers_perSample.tsv"]
    elif config['taxonomy']['dada']['do'] and 'ASV' in config['taxonomy']['dada']['run_on']:
        physInputs = ["sequenceTables/all.seqTab.tax.RDS","reporting/finalNumbers_perSample.tsv"]
    else:
        physInputs = ["sequenceTables/all.seqTab.RDS","reporting/finalNumbers_perSample.tsv"]
    if config['blast']['do'] and config['blast']['run_basta'] and 'ASV' in config['blast']['run_on']:
        physInputs = ["sequenceTables/all.seqTab.tax.blast.RDS","reporting/finalNumbers_perSample.tsv"]
    if 'postclustering' in STEPS:
        physInputs.append("clusteredTables/cluster_info.tsv")
    rule phyloseq_handoff_tax:
        input:
            physInputs
        output:
            "sequenceTables/all.seqTab.phyloseq.RDS"
        threads: 1
        params:
            currentStep = "taxonomy",
            what = "ASV",
            tree = "none",
            cluster = "add" if 'postclustering' in STEPS else "none",
            clusterMeth = config['post_clustering']['method']
        resources:
            runtime="12:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "add_R_env.yml"
        log: "logs/phyloseq_hand-off.log"
        script:
            SCRIPTSDIR+"phyloseq_handoff.R"

    if config['taxonomy']['dada']['do'] and 'cluster' in config['taxonomy']['dada']['run_on']:
        physInputsCl = ["clusteredTables/clusteredTab.tax.RDS","reporting/finalNumbers_perSample.tsv"]
    elif config['taxonomy']['decipher']['do'] and 'cluster' in config['taxonomy']['decipher']['run_on']:
        physInputsCl = ["clusteredTables/clusteredTab.tax.RDS","reporting/finalNumbers_perSample.tsv"]
    elif config['taxonomy']['mothur']['do'] and 'cluster' in config['taxonomy']['mothur']['run_on']:
        physInputsCl = ["clusteredTables/clusteredTab.tax.RDS","reporting/finalNumbers_perSample.tsv"]
    else:
        physInputsCl = ["clusteredTables/clusteredTab.RDS","reporting/finalNumbers_perSample.tsv"]
    if config['blast']['do'] and config['blast']['run_basta'] and 'cluster' in config['blast']['run_on']:
        physInputsCl = ["clusteredTables/clusteredTab.tax.blast.RDS","reporting/finalNumbers_perSample.tsv"]
    rule phyloseq_handoff_postCl:
        input:
            physInputsCl
        output:
            "clusteredTables/clusteredTab.phyloseq.RDS"
        threads: 1
        params:
            currentStep = "taxonomy",
            what = "cluster",
            tree = "none",
            cluster = "none"
        resources:
            runtime="12:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "add_R_env.yml"
        log: "logs/phyloseq_hand-off.cluster.log"
        script:
            SCRIPTSDIR+"phyloseq_handoff.R"


