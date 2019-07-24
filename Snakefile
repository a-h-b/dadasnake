# include configuration file
#default configuration file (hard-coded, not pretty for portability)
configfile: "/data/project/metaamp/TOOLS/dada_pipe/config.default.yaml"

include:
    "dada_scripts/get_config.rules"

workdir:
    OUTPUTDIR

if not config['skip_db']:
    include:
        "dada_scripts/cut_db.rules"
if 'primers' in STEPS:
    include:
        "dada_scripts/cutadapt.rules"
if 'dada' in STEPS:
    include:
        "dada_scripts/dada.rules"
if 'taxonomy' in STEPS:
    include:
        "dada_scripts/taxonomy.rules"
if 'postprocessing' in STEPS:
    include:
        "dada_scripts/post.rules"


inputs = []
if not config['skip_db']:
    inputs.append('cut_db.done')
if 'primers' in STEPS:
    inputs.append('primers.done')
if 'dada' in STEPS:
    inputs.append('dada.done')
if 'taxonomy' in STEPS:
    inputs.append('taxonomy.done')
if 'postprocessing' in STEPS:
    inputs.append('postprocessing.done')


# master command
rule ALL:
    input:
        inputs
    output:
        touch('workflow.done')
