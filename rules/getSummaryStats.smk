#Global Wildcards
SAMPLE_DIRS = [x.split('_')[0] for x in glob_wildcards(config["input_path"] + "/{sample}_R1_001.fastq.gz").sample]
SAMPLES = glob_wildcards(config["input_path"] + "/{sample}_R1_001.fastq.gz").sample


rule get_summary_stats:
    input:
        all_stats = expand(config["output_path"] + "/{sample}/{sample}.stat", sample=SAMPLES)
    output:
        sum_stats = config["output_path"] + "/Summary/Stats.csv"        
    log:
        "logs/sum_stats.log"
    shell:
        r"""
        echo "Sample,Total reads,Mapped reads,Mapped reads(>30nt), Mapped read %,Ref name,Ref Length,Ref Coverage,Ref coverage %,Average Depth,Min Depth,Max Depth" > {output.sum_stats};

        sample_array=({input.all_stats})

        for sample_stat in ${{sample_array[@]}}
        do
        sample=$(echo $sample_stat | cut -d / -f2 | cut -d _ -f1)   #MAY NEED TO CHANGE THIS LINE WHEN I FIGURE OUT HOW TO SHORTEN THE SAMPLE DIRECTORY NAMES
        echo -n "$sample,";

        cat $sample_stat |awk '{{for(i=1;i<=NF;i++) printf "%s, ",$i; print ""; }}'

        done >> {output.sum_stats} \
        2>{log}
        """
