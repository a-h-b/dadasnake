if config['hand_off']['biom'] and (not config['do_taxonomy'] or (not config['taxonomy']['decipher']['do'] and not config['taxonomy']['mothur']['do'] and not config['taxonomy']['dada']['do'])):
    rule big_biom_handoff_preTax:
        input:
            "sequenceTables/all.seqTab.RDS",
            "reporting/finalNumbers_perSample.tsv"
        output:
            "sequenceTables/all.seqTab.biom"
        threads: 1
        params:
            currentStep = "dada"
        resources:
            runtime="12:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "add_R_env.yml"
        log: "logs/biom_hand-off.log"
        script:
            SCRIPTSDIR+"biom_handoff.R"

if config['hand_off']['phyloseq']:
    if not config['do_taxonomy'] and not config['do_postprocessing']:
        physInputs = ["sequenceTables/all.seqTab.tax.RDS","reporting/finalNumbers_perSample.tsv"]
        if 'postclustering' in STEPS:
            physInputs.append("clusteredTables/cluster_info.tsv")
        rule big_phyloseq_handoff_preTax:
            input:
                physInputs
            output:
                "sequenceTables/all.seqTab.phyloseq.RDS"
            threads: 1
            params:
                currentStep = "dada",
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

        physInputsCl = ["clusteredTables/clusteredTab.RDS","reporting/finalNumbers_perSample.tsv"]
        rule big_phyloseq_handoff_preTaxCl:
            input:
                physInputsCl
            output:
                "clusteredTables/clusteredTab.phyloseq.RDS"
            threads: 1
            params:
                currentStep = "dada",
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


