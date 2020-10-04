# include configuration file
configfile: srcdir("config/config.default.yaml")

SCRIPTSDIR = srcdir("workflow/scripts/")
ENVDIR = srcdir("workflow/envs/")

include:
    "workflow/rules/get_config.smk"

workdir:
    OUTPUTDIR

f = open('full.config.yaml', 'w+')
yaml.dump(config, f, allow_unicode=True,default_flow_style=False)

if config['paired']:
    if 'primers' in STEPS:
        include:
            "workflow/rules/cutadapt.smk"
    else:
        include:
            "workflow/rules/copying.smk"
    if 'dada' in STEPS:
        if not config['dada']['pool']:
            include:
                "workflow/rules/dada.paired.smk"
        else:
            include:
                "workflow/rules/dada.paired.pool.smk"
else:
    if 'primers' in STEPS:
        include:
            "workflow/rules/cutadapt.single.smk"
    else:
        include:
            "workflow/rules/copying.single.smk"
    if 'dada' in STEPS:
        if not config['dada']['pool']:
            include:
                "workflow/rules/dada.single.smk"
        else:
            include:
                "workflow/rules/dada.single.pool.smk"
if 'dada' in STEPS:
    include:
        "workflow/rules/dada.common.smk"
if 'taxonomy' in STEPS:
    include:
        "workflow/rules/taxonomy.smk"
if 'postprocessing' in STEPS:
    if config['final_table_filtering']['do']:
        include:
            "workflow/rules/post.filtering.smk"
    else:
        include:
            "workflow/rules/post.no_filtering.smk"



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
        shell("mkdir -p job.errs.outs &>> logs/cleanup.log; ( mv slurm* job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv snakejob.* job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *log job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *logfile job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log")
else:
    onsuccess:
        shell('mkdir -p job.errs.outs &>> logs/cleanup.log; ( mv slurm* job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; (  mv snakejob.* job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *log job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *logfile job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; echo "$(date) {config[sessionName]}" | mail -s "dadasnake finished" {EMAIL} ')
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
