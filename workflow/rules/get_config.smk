import os
import shutil
import gzip
import yaml
from re import match
from copy import deepcopy
import subprocess
import pandas as pd
from snakemake.utils import validate

# default executable for snakemake
shell.executable("bash")


#validate(config,schema="schemas/config.schema.yaml")

print(config['sessionKind'])
if config['sessionKind'] == "cluster":
    if str(config['settingsLocked']).lower() in ['true','t']:
        config['settingsLocked'] = True
        print("The person who set up dadasnake disabled changing resource settings. The original settings will be used:")
        print("normalMem: " + config['normalMem'])
        print("bigMem: " + config['bigMem'])
        print("bigCores: " + str(config['bigCores']))
    else:
        config['settingsLocked'] = False
    if config['normalMem'] == "":
        if config['settingsLocked']:
            raise Exception("You're attempting to submit dadasnake to a cluster and haven't specified the memory requirements. Ask the person who set up dadasnake to specify NORMAL_MEM_EACH in " + ROOTDIR + "VARIABLE_CONFIG")
        else:
            raise Exception("You're attempting to submit dadasnake to a cluster and haven't specified the memory requirements. Set normalMem in your config file.")
    if config['bigCores'] == "" or int(config['bigCores']) == 0:
        if config['big_data']:
            if config['settingsLocked']:
                raise Exception("You're attempting to submit a big dataset in dadasnake to a cluster, but haven't specified more than 0 cores to be used. Ask the person who set up dadasnake to change BIGMEM_CORES in " + ROOTDIR + "VARIABLE_CONFIG to run a big data set, or set big_data to false in your config file.") 
            else:
                raise Exception("You're attempting to submit a big dataset in dadasnake to a cluster, but haven't specified more than 0 cores to be used. Change bigCores in the config file to run a big data set, or set big_data to false in your config file.")
        else:
            config['bigMem'] = config['normalMem']
            config['bigCores'] = workflow.nodes
            print("You haven't specified more than 0 bigmem cores, all rules will be performed on normal cores with " + config['normalMem'] + ".")
    elif config['bigMem'] == "":
        if config['big_data']:
            if config['settingsLocked']:
                raise Exception("You're attempting to submit a big dataset in dadasnake to a cluster, but haven't specified more than 0 cores to be used. Ask the person who installed dadasnake to change BIGMEM_MEM_EACH in " + ROOTDIR + "VARIABLE_CONFIG to run a big data set, or set big_data to false in your config file, if your dataset isn't that big.")
            else:
                raise Exception("You're attempting to submit a big dataset in dadasnake to a cluster, but haven't specified more than 0 cores to be used. Change bigMem to an adequate number in the config file to run a big data set, or set big_data to false in your config file.")
        else:
            config['bigMem'] = config['normalMem']
            config['bigCores'] = workflow.cores
            print("You haven't specified the size of your bigmem cores, all rules will be performed on normal cores with " + config['normalMem'] + ".")
# get dryruns to have realistic output:
elif config['sessionKind'] == "dryrun":
    if str(config['settingsLocked']).lower() in ['true','t']:
        config['settingsLocked'] = True
        print("The person who set up dadasnake disabled changing resource settings. The original settings will be used:")
        print("normalMem: " + config['normalMem'])
        print("bigMem: " + config['bigMem'])
        print("bigCores: " + str(config['bigCores']))
    else:
        config['settingsLocked'] = False
    if config['normalMem'] == "":
        if config['settingsLocked']:
            print("You will not be able to submit dadasnake to a cluster unless you ask the person who set up dadasnake to specify NORMAL_MEM_EACH in " + ROOTDIR + "VARIABLE_CONFIG")
        else:
            print("You will not be able to submit dadasnake to a cluster unless you set normalMem in your config file.")
    if config['bigCores'] == "" or int(config['bigCores']) == 0:
        if config['big_data']:
            if config['settingsLocked']:
                print("You will not be able to submit a big dataset in dadasnake to a cluster, unless you ask the person who set up dadasnake to change BIGMEM_CORES in " + ROOTDIR + "VARIABLE_CONFIG to run a big data set, or set big_data to false in your config file.")
            else:
                print("You will not be able to submit a big dataset in dadasnake to a cluster, unless you change bigCores in the config file to run a big data set, or set big_data to false in your config file.")
        else:
            config['bigMem'] = config['normalMem']
            config['bigCores'] = workflow.cores
            print("You haven't specified more than 0 bigmem cores, in cluster mode, all rules would be performed on normal cores with " + config['normalMem'] + ".")
    elif config['bigMem'] == "":
        if config['big_data']:
            if config['settingsLocked']:
                print("You will not be able to submit a big dataset in dadasnake to a cluster, unless you ask the person who installed dadasnake to change BIGMEM_MEM_EACH in " + ROOTDIR + "VARIABLE_CONFIG to run a big data set, or set big_data to false in your config file, if your dataset isn't that big.")
            else:
                print("You will not be able to submit a big dataset in dadasnake to a cluster, unless you change bigMem to an adequate number in the config file to run a big data set, or set big_data to false in your config file.")
        else:
            config['bigMem'] = config['normalMem']
            config['bigCores'] = workflow.cores
            print("You haven't specified the size of your bigmem cores, in cluster mode, all rules would be performed on normal cores with " + config['normalMem'] + ".")
    config['bigCores'] = workflow.cores
#for non-cluster runs:
else:
    config['bigCores'] = workflow.cores
print("Final resource settings:")
print("maxCores: " + str(workflow.cores))
#print("normalMem: " + config['normalMem'])
#print("bigMem: " + config['bigMem'])
#print("bigCores: " + str(config['bigCores']))



# get parameters from the command line
OUTPUTDIR = os.path.expandvars(config['outputdir'])

EMAIL = config['email']
if EMAIL != "":
    if not re.fullmatch(r"^\w+([\.-]?\w+)*@\w+([\.-]?\w+)*(\.\w{2,3})+$", EMAIL):
        EMAIL = ""
        print("Your email address is not valid, you will not receive notifications.")

if os.path.isabs(os.path.expandvars(config['sample_table'])):
    sam_path = os.path.expandvars(config['sample_table'])
else:
    sam_path = os.getcwd() + "/" + os.path.expandvars(config['sample_table'])
try:
    samples = pd.read_table(sam_path)
except:
    print("Sample table was not found. Please enter the absolute path and file name in the config file.")
    raise
if 'run' in samples.columns:
    samples['run'] = samples['run'].astype(str)
    for r in samples['run']:
        if not re.match(r"[0-9a-zA-Z]",r):
            raise Exception('please start run names with a letter or number')
else:
    samples['run'] = ["run1"] * samples.shape[0]
    print("adding column with run info")
if samples[['library','run']].duplicated().any():
    raise Exception('names in library should be unique within runs.')
samples = samples.set_index(["library","run"],drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels]) 
samples['library'] = samples['library'].astype(str)
for lib in samples['library']:
    if not re.match(r"[0-9a-zA-Z]",lib):
        raise Exception('please start library names with a letter or number')
if 'sample' not in samples.columns:
    samples['sample'] = samples['library']
    print("adding column with sample names based on library names")
else:
    samples['sample'] = samples['sample'].astype(str)
    for r in samples['run'].unique():
        tmpsr=samples[samples['run']==r]
        for i in range(tmpsr.shape[0]):
            csam = tmpsr['sample'][i]
            if not re.match(r"[0-9a-zA-Z]",csam):
                raise Exception('please start sample names with a letter or number')
            for j in range(tmpsr.shape[0]):
                clib = tmpsr['library'][j]
                if csam==clib and i==j:
                    csl = tmpsr.loc[tmpsr["sample"]==csam,"library"]
                    if len(csl) > 1:
                        raise Exception('names in library should differ from names of samples that have multiple libraries in the same run.')
                elif csam==clib:
                    raise Exception('names in library should differ from unrelated sample names.')    
if 'r1_file' not in samples.columns:
    raise Exception("You haven't provided file names for read 1 - column should be named r1_file.")
if config['paired'] and 'r2_file' not in samples.columns:
    raise Exception("You haven't provided file names for read 2 - column should be named r2_file.")

if os.path.isabs(os.path.expandvars(config['raw_directory'])):
    RAW = os.path.expandvars(config['raw_directory'])
else:
    RAW = os.getcwd() + "/" + os.path.expandvars(config['raw_directory'])
config['raw_directory'] = RAW


if config['big_data'] and config['dada']['pool']:
    raise Exception("No workflow is defined for pooled analysis of large datasets. Use per-sample analysis for large datasets by setting dada:pool: to true.")

#validate(samples, schema="schemas/sample.schema.yaml")

yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
yaml.add_representer(tuple, lambda dumper, data: dumper.represent_sequence('tag:yaml.org,2002:seq', data))

PRELIM_STEPS = ''
if config['do_primers']:
    PRELIM_STEPS += "primers "
if config['do_dada']:
    PRELIM_STEPS += "dada "
    if config['do_taxonomy']:
        PRELIM_STEPS += "taxonomy "
    if config['do_postprocessing']:
        PRELIM_STEPS += "postprocessing "

STEPS = PRELIM_STEPS.split()

# temporary directory will be stored inside the OUTPUTDIR directory
# unless an absolute path is set
TMPDIR = os.path.expandvars(config['tmp_dir'])
if not os.path.isabs(TMPDIR):
    TMPDIR = os.path.abspath(os.path.join(OUTPUTDIR, TMPDIR))
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)

if config['do_taxonomy']:
    if config['taxonomy']['dada']['do']:
        DADADB_OLD = os.path.join(config['taxonomy']['dada']['db_path'], config['taxonomy']['dada']['refFasta'])
        DADADB_MULT = config['taxonomy']['dada']['ref_dbs_full'].split()
        DADADB_NAMES = config['taxonomy']['dada']['db_short_names'].split()
        if config['taxonomy']['dada']['look_for_species']:
               DADA_SPEC = config['taxonomy']['dada']['spec_db'].split()
        else:
               DADA_SPEC = []
        if DADADB_MULT:
            if DADADB_OLD:
                print("Warning: you have given configuration for both ref_dbs_full and db_path / refFasta in taxonomy: dada . Only the ref_dbs_full database will be used.")
                DADADB_OLD = ""
            if DADADB_NAMES:
                DADADB = dict(zip(DADADB_NAMES,DADADB_MULT))
                if DADA_SPEC:
                    if len(DADA_SPEC) == 1:
                        DADADB_SPEC = dict(zip(DADADB_NAMES,[DADA_SPEC] * len(DADADB_MULT)))
                    else:
                        DADADB_SPEC = dict(zip(DADADB_NAMES,DADA_SPEC))    
                else:
                    DADADB_SPEC = {}
            else:
                DADADB_FILES = [os.path.basename(s) for s in DADADB_MULT]
                DADADB = dict(zip(DADADB_FILES,DADADB_MULT))
                if DADA_SPEC:
                    if len(DADA_SPEC) == 1:
                        DADADB_SPEC = dict(zip(DADADB_FILES,[DADA_SPEC] * len(DADADB_MULT)))
                    else:
                        DADADB_SPEC = dict(zip(DADADB_FILES,DADA_SPEC)) 
                else:
                    DADADB_SPEC = {}
        elif DADADB_OLD:
            if DADADB_NAMES:
                DADADB = dict(zip(DADADB_NAMES,[DADADB_OLD]))
                if DADA_SPEC:
                    DADADB_SPEC = dict(zip(DADADB_NAMES,DADA_SPEC))
                else:
                    DADADB_SPEC = {}        
            else:
                DADADB = dict(zip([os.path.basename(DADADB_OLD)],[DADADB_OLD]))
                if DADA_SPEC:
                    DADADB_SPEC = dict(zip([os.path.basename(DADADB_OLD)],DADA_SPEC))
                else:
                    DADADB_SPEC = {}        
        else:
            DADADB = {}
            DADADB_SPEC = {}
    else:
        DADADB = {}
        DADADB_SPEC = {}
    if config['taxonomy']['decipher']['do']:
        DECIDB_OLD = os.path.join(config['taxonomy']['decipher']['db_path'], config['taxonomy']['decipher']['tax_db'])
        DECIDB_MULT = config['taxonomy']['decipher']['ref_dbs_full'].split()
        DECIDB_NAMES = config['taxonomy']['decipher']['db_short_names'].split()
        if config['taxonomy']['decipher']['look_for_species']:
               DECI_SPEC = config['taxonomy']['decipher']['spec_db'].split()
        else:
               DECI_SPEC = []
        if DECIDB_MULT:
            if DECIDB_OLD:
                print("Warning: you have given configuration for both ref_dbs_full and db_path / tax_db in taxonomy: decipher . Only the ref_dbs_full database will be used.")
                DECIDB_OLD = ""
            if DECIDB_NAMES:
                DECIDB = dict(zip(DECIDB_NAMES,DECIDB_MULT))
                if DECI_SPEC:
                    if len(DECI_SPEC) == 1:
                        DECIDB_SPEC = dict(zip(DECIDB_NAMES,[DECI_SPEC] * len(DECIDB_MULT)))
                    else:
                        DECIDB_SPEC = dict(zip(DECIDB_NAMES,DECI_SPEC))
                else:
                    DECIDB_SPEC = {}
            else:
                DECIDB_FILES = [os.path.basename(s) for s in DECIDB_MULT]
                DECIDB = dict(zip(DECIDB_FILES,DECIDB_MULT))
                if DECI_SPEC:
                    if len(DECI_SPEC) == 1:
                        DECIDB_SPEC = dict(zip(DECIDB_FILES,[DECI_SPEC] * len(DECIDB_MULT)))
                    else:
                        DECIDB_SPEC = dict(zip(DECIDB_FILES,DECI_SPEC))
                else:
                    DECIDB_SPEC = {}
        elif DECIDB_OLD:
            if DECIDB_NAMES:
                DECIDB = dict(zip(DECIDB_NAMES,[DECIDB_OLD]))
                if DECI_SPEC:
                    DECIDB_SPEC = dict(zip(DECIDB_NAMES,DECI_SPEC))
                else:
                    DECIDB_SPEC = {}
            else:
                DECIDB = dict(zip([os.path.basename(DECIDB_OLD)],[DECIDB_OLD]))
                if DECI_SPEC:
                    DECIDB_SPEC = dict(zip([os.path.basename(DECIDB_OLD)],DECI_SPEC))
                else:
                    DECIDB_SPEC = {}
        else:
            DECIDB = {}
            DECIDB_SPEC = {}
    else:
        DECIDB = {}
        DECIDB_SPEC = {}
    if config['taxonomy']['mothur']['do']:
        MOTHDB_OLD = os.path.join(config['taxonomy']['mothur']['db_path'], config['taxonomy']['mothur']['tax_db'])
        MOTHDB_MULT = config['taxonomy']['mothur']['ref_dbs_full'].split()
        MOTHDB_NAMES = config['taxonomy']['mothur']['db_short_names'].split()
        if MOTHDB_MULT:
            if MOTHDB_OLD:
                print("Warning: you have given configuration for both ref_dbs_full and db_path / tax_db in taxonomy: mothur . Only the ref_dbs_full database will be used.")
                MOTHDB_OLD = ""
            if MOTHDB_NAMES:
                MOTHPATH = dict(zip(MOTHDB_NAMES,[os.path.dirname(s) for s in MOTHDB_MULT]))
                MOTHDB = dict(zip(MOTHDB_NAMES,[os.path.basename(s) for s in MOTHDB_MULT]))
            else:
                MOTHDB_FILES = [os.path.basename(s) for s in MOTHDB_MULT]
                MOTHPATH = dict(zip(MOTHDB_FILES,[os.path.dirname(s) for s in MOTHDB_MULT]))
                MOTHDB = dict(zip(MOTHDB_FILES,MOTHDB_FILES))
        elif MOTHDB_OLD:
            if MOTHDB_NAMES:
                MOTHPATH = dict(zip(MOTHDB_NAMES,[os.path.dirname(MOTHDB_OLD)]))
                MOTHDB = dict(zip(MOTHDB_NAMES,[os.path.basename(MOTHDB_OLD)]))
            else:
                MOTHPATH = dict(zip([os.path.basename(MOTHDB_OLD)],[os.path.dirname(MOTHDB_OLD)]))
                MOTHDB = dict(zip([os.path.basename(MOTHDB_OLD)],[os.path.basename(MOTHDB_OLD)]))
        else:
            MOTHDB = {}
            MOTHPATH = {}
    else:
        MOTHDB = {}
        MOTHPATH = {}

#print(MOTHDB)
