# include configuration file
configfile: srcdir("config.default.yaml")

SRC_dir = srcdir("dada_scripts/")

include:
    "dada_scripts/get_config.smk"

workdir:
    OUTPUTDIR

f = open('full.config.yaml', 'w+')
yaml.dump(config, f, allow_unicode=True,default_flow_style=False)

if config['paired']:
    if 'primers' in STEPS:
        include:
            "dada_scripts/cutadapt.smk"
    else:
        include:
            "dada_scripts/copying.smk"
    if 'dada' in STEPS:
        if not config['dada']['pool']:
            include:
                "dada_scripts/dada.paired.smk"
        else:
            include:
                "dada_scripts/dada.paired.pool.smk"
else:
    if 'primers' in STEPS:
        include:
            "dada_scripts/cutadapt.single.smk"
    if 'dada' in STEPS:
        if not config['dada']['pool']:
            include:
                "dada_scripts/dada.single.smk"
        else:
            include:
                "dada_scripts/dada.single.pool.smk"
if 'dada' in STEPS:
    include:
        "dada_scripts/dada.common.smk"
if 'taxonomy' in STEPS:
    include:
        "dada_scripts/taxonomy.smk"
if 'postprocessing' in STEPS:
    if config['final_table_filtering']['do']:
        include:
            "dada_scripts/post.filtering.smk"
    else:
        include:
            "dada_scripts/post.no_filtering.smk"



inputs = []
if 'primers' in STEPS:
    inputs.append('primers.done')
if 'dada' in STEPS:
    inputs.append('dada.done')
if 'taxonomy' in STEPS:
    inputs.append('taxonomy.done')
if 'postprocessing' in STEPS:
    inputs.append('postprocessing.done')
if config['hand_off']['biom']:
    inputs.append('sequenceTables/all.seqTab.biom')

if EMAIL == "":
    onsuccess:
        shell("mkdir -p job.errs.outs; mv slurm* job.errs.outs || touch job.errs.outs; mv snakejob.* job.errs.outs || touch job.errs.outs; mv *log job.errs.outs || touch job.errs.outs; mv *logfile job.errs.outs || touch job.errs.outs")
else:
    onsuccess:
        shell('mkdir -p job.errs.outs; mv slurm* job.errs.outs || touch job.errs.outs; mv snakejob.* job.errs.outs || touch job.errs.outs; mv *log job.errs.outs || touch job.errs.outs; mv *logfile job.errs.outs || touch job.errs.outs; echo "$(date) {config[sessionName]}" | mail -s "dadasnake finished" {EMAIL} ')
    onerror:
        shell('echo "$(date) {config[sessionName]}" | mail -s "dadasnake exited with error" {EMAIL} ')
    onstart:
        shell('echo "$(date) {config[sessionName]}" | mail -s "dadasnake started" {EMAIL} ')

localrules: ALL, SamplesPrint

# master command
rule ALL:
    input:
        inputs
    output:
        touch('workflow.done')

rule SamplesPrint:
    input:
        sam_path
    output:
        "reporting/sample_table.tsv"
    run:
        samples.to_csv(path_or_buf=output[0],sep="\t",index=False,index_label=False)
