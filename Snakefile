# include configuration file
configfile: srcdir("config.default.yaml")

SRC_dir = srcdir("dada_scripts/")

include:
    "dada_scripts/get_config.rules"

workdir:
    OUTPUTDIR

f = open(OUTPUTDIR+'/full.config.yaml', 'w+')
yaml.dump(config, f, allow_unicode=True,default_flow_style=False)

if config['paired']:
    if 'primers' in STEPS:
        include:
            "dada_scripts/cutadapt.rules"
    else:
        include:
            "dada_scripts/copying.rules"
    if 'dada' in STEPS:
        include:
            "dada_scripts/dada.paired.rules"
else:
    if 'primers' in STEPS:
        include:
            "dada_scripts/cutadapt.single.rules"
    if 'dada' in STEPS:
        include:
            "dada_scripts/dada.single.rules"
if 'dada' in STEPS:
    include:
        "dada_scripts/dada.common.rules"
if 'taxonomy' in STEPS:
    include:
        "dada_scripts/taxonomy.rules"
if 'postprocessing' in STEPS:
    if config['final_table_filtering']['do']:
        include:
            "dada_scripts/post.filtering.rules"
    else:
        include:
            "dada_scripts/post.no_filtering.rules"



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
        shell("mv snakejob.* *log job.errs.outs || (  mkdir job.errs.outs && mv snakejob.* *log job.errs.outs )")
else:
    onsuccess:
        shell("mv snakejob.* *log job.errs.outs || (  mkdir job.errs.outs && mv snakejob.* *log job.errs.outs ); date | mail -s 'dadasnake finished' {EMAIL} ")
    onerror:
        shell("date | mail -s 'dadasnake exited with error' {EMAIL} ")
    onstart:
        shell("date | mail -s 'dadasnake started' {EMAIL} ")


# master command
rule ALL:
    input:
        inputs
    output:
        touch('workflow.done')
    threads: 1
    params:
        runtime="00:10:00",
        mem="8G"

#rule ConfigPrint:


rule SamplesPrint:
    input:
        config['sample_table']
    output:
        "reporting/sample_table.tsv"
    threads: 1
    params:
        runtime="00:10:00",
        mem="8G"
    run:
        samples.to_csv(path_or_buf=output[0],sep="\t",index=False,index_label=False)
