if config['hand_off']['biom'] and (not config['do_taxonomy'] or (not config['taxonomy']['decipher']['do'] and not config['taxonomy']['mothur']['do'] and not config['taxonomy']['dada']['do'])):
    rule bigbiom_handoff_preTax:
        input:
            "sequenceTables/all.seqTab.RDS",
            "reporting/finalNumbers_perSample.tsv"
        output:
            "sequenceTables/all.seqTab.biom"
        threads: config['bigCores']
        params:
            currentStep = "dada"
        resources:
            runtime="12:00:00",
            mem=config['bigMem']
        conda: ENVDIR + "add_R_env.yml"
        log: "logs/bigbiom_hand-off.log"
        script:
            SCRIPTSDIR+"biom_handoff.R"

if config['hand_off']['phyloseq']:
    if not config['do_taxonomy'] or (not config['taxonomy']['decipher']['do'] and not config['taxonomy']['mothur']['do'] and not config['taxonomy']['dada']['do']):
        if not config['do_postprocessing']:
            rule bigphyloseq_handoff_preTax:
                input:
                    "sequenceTables/all.seqTab.RDS",
                    "reporting/finalNumbers_perSample.tsv"
                output:
                    "sequenceTables/all.seqTab.phyloseq.RDS"
                threads: config['bigCores']
                params:
                    currentStep = "dada"
                resources:
                    runtime="12:00:00",
                    mem=config['bigMem']
                conda: ENVDIR + "add_R_env.yml"
                log: "logs/bigphyloseq_hand-off.log"
                script:
                    SCRIPTSDIR+"phyloseq_handoff.R"

