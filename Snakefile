
##### Configuration #####
configfile: str(workflow.current_basedir) + "/config.yaml"

if config.get("input_path"):
    config["input_path"] = config["input_path"].rstrip('/')

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip('/')


#Wildcards
SAMPLE_DIRS = [x.split('_')[0] for x in glob_wildcards(config["input_path"] + "/{sample}_R1_001.fastq.gz").sample]
SAMPLES = glob_wildcards(config["input_path"] + "/{sample}_R1_001.fastq.gz").sample


##### Create Outputs #####
rule all:
    input:
        trimmed = expand(config["output_path"] + "/{sample}/{sample}_val_1.fq",sample=SAMPLES)

##### Modules #####
include: "rules/readAlign-Illumina-NimaGen.smk"

