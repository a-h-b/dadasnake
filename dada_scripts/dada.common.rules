if config['hand_off']['biom'] and (not config['do_taxonomy'] or (not config['taxonomy']['decipher']['do'] and not config['taxonomy']['mothur']['do'])):
    rule biom_handoff_preTax:
        input:
            "sequenceTables/all.seqTab.RDS",
            "reporting/finalNumbers_perSample.tsv"
        output:
            "sequenceTables/all.seqTab.biom"
        threads: 1
        params:
            currentStep = "dada",
            mem="8G",
            runtime="12:00:00"
        conda: "dada_env.yml"
        log: "logs/biom_hand-off.log"
        script:
            SRC_dir+"biom_handoff.R"

if config['hand_off']['phyloseq']:
    if not config['do_taxonomy'] or (not config['taxonomy']['decipher']['do'] and not config['taxonomy']['mothur']['do']):
        if not config['do_postprocessing']:
            rule phyloseq_handoff_preTax:
                input:
                    "sequenceTables/all.seqTab.RDS",
                    "reporting/finalNumbers_perSample.tsv"
                output:
                    "sequenceTables/all.seqTab.phyloseq.RDS"
                threads: 1
                params:
                    currentStep = "dada",
                    mem="8G",
                    runtime="12:00:00"
                conda: "dada_env.yml"
                log: "logs/phyloseq_hand-off.log"
                script:
                    SRC_dir+"phyloseq_handoff.R"

