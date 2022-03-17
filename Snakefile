
##### Configuration #####
configfile: str(workflow.current_basedir) + "/config.yaml"

if config.get("input_path"):
    config["input_path"] = config["input_path"].rstrip('/')

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip('/')


##### Get Sample IDs from Sample Sheet #####
sample_sheet = open(config["sample_sheet"], 'r')
is_sample_id = False
SAMPLES = ['Undetermined_S0']   #Inserting Undetermined manually since it won't be in sample sheet
sample_count = 0
for row in sample_sheet.readlines():
    if is_sample_id == False:
        if "Sample_ID" in row:
            is_sample_id = True
            continue
        else:
            continue
    elif is_sample_id == True:
        sample_count += 1
        SAMPLES.append(row.split(',')[0] + '_S' + str(sample_count))
sample_sheet.close()


#Testing set
#Uncomment below line to only run pipeline with 10 samples
#SAMPLES = SAMPLES[30:40]


#Shortened sample name
SAMPLE_DIRS = [id.split('_')[0] for id in SAMPLES]
#SAMPLE_DIRS = ['_'.join(id.split('_')[0:2]) for id in SAMPLES]   #Used for waste water samples to include sample condition

#Dictionary of shortened sample name to full sample name
SAMPLES_DICT = {sample.split('_')[0]:sample for sample in SAMPLES}
#SAMPLES_DICT = {'_'.join(sample.split('_')[0:2]):sample for sample in SAMPLES}

##### Request Outputs #####
rule all:
    input:
#Summary outputs
        sum_stats = config["summary_dir"] + "/Stats.csv",
        sum_amp_cov = config["summary_dir"] + "/ampliconDepth.csv",
        collected_amp_dep = expand(config["summary_dir"] + "/CoveragePlots/{sample_dir}-Amplicon-Depth.pdf", sample_dir=SAMPLE_DIRS),
        ambig_pos_and_dep = config["summary_dir"] + "/ambig_nuc_pos_and_dep.csv",
        masked_ambig_nucs = config["summary_dir"] + "/masked_ambig_nuc_pos.csv",
        pangolin_report =  config["summary_dir"] + "/lineage_report.csv",
#Climb Dir
        climb_bam = expand(config["batch_dir"] + "/ClimbSeq/{sample}/{sample}.bam", sample=SAMPLE_DIRS),
#Testing freyja rules
#        freyja_summarised_plot =  config["summary_dir"] + "/All-freyja-summarised-plot.pdf",
#        freyja_lineage_plot =  config["summary_dir"] + "/All-freyja-lineage-plot.pdf",
#        temp_aggregated = config["summary_dir"] + "/All-freyja-demixed_temp.tsv"


##### Modules #####
include: "rules/readAlign-Illumina-NimaGen.smk"
include: "rules/getSampleStats.smk"
include: "rules/getSummaryStats.smk"
include: "rules/postReadAlign.smk"
