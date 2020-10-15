if config['hand_off']['biom'] and (not config['do_taxonomy'] or (not config['taxonomy']['decipher']['do'] and not config['taxonomy']['mothur']['do'] and not config['taxonomy']['dada']['do'])):
    rule biom_handoff_preTax:
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
            mem=config['normalMem']
        conda: ENVDIR + "add_R_env.yml"
        log: "logs/biom_hand-off.log"
        script:
            SCRIPTSDIR+"biom_handoff.R"

if config['hand_off']['phyloseq']:
    if not config['do_taxonomy'] or (not config['taxonomy']['decipher']['do'] and not config['taxonomy']['mothur']['do'] and not config['taxonomy']['dada']['do']):
        if not config['do_postprocessing']:
            rule phyloseq_handoff_preTax:
                input:
                    "sequenceTables/all.seqTab.RDS",
                    "reporting/finalNumbers_perSample.tsv"
                output:
                    "sequenceTables/all.seqTab.phyloseq.RDS"
                threads: 1
                params:
                    currentStep = "dada"
                resources:
                    runtime="12:00:00",
                    mem=config['normalMem']
                conda: ENVDIR + "add_R_env.yml"
                log: "logs/phyloseq_hand-off.log"
                script:
                    SCRIPTSDIR+"phyloseq_handoff.R"

