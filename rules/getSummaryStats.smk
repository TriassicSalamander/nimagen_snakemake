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


rule get_summary_amplicon_coverage:
    input:
        all_amp_cov = expand(config["output_path"] + "/{sample}/{sample}.amplicon.cov", sample=SAMPLES)
    output:
        sum_amp_cov = config["output_path"] + "/Summary/ampliconDepth.csv"
    log:
        "logs/sum_amp_cov.log"
    shell:
        r"""
        sample_array=({input.all_amp_cov})

        samples=()
        for sample_stat in ${{sample_array[@]}}
        do
        sample=$(echo $sample_stat | cut -d / -f2 | cut -d _ -f1)   #MAY NEED TO CHANGE THIS SECTION WHEN I FIGURE OUT HOW TO SHORTEN THE SAMPLE DIRECTORY NAMES
        samples+=( $sample )
        done 2>{log}

        echo ${{samples[*]}} | sed 's/ /,/g' > {output.sum_amp_cov} 2>{log}

        paste {input.all_amp_cov} | awk '{{for(i=4;i<=NF;i+=4) printf "%s,",$i; print "";}}' >> {output.sum_amp_cov} 2>{log}
        """


rule collect_plots:
    input:
        all_plots = expand(config["output_path"] + "/{sample}/{sample}-{plot_type}.pdf", sample=SAMPLES, plot_type=["Amplicon-Depth","coverage"])
    output:
        collected_plots = expand(config["output_path"] + "/Summary/CoveragePlots/{sample}-{plot_type}.pdf", sample=SAMPLES, plot_type=["Amplicon-Depth","coverage"])
    params:
        plot_directory = config["output_path"] + "/Summary/CoveragePlots/"
    log:
        "logs/collect_plots.log"
    shell:
        r"""
        sample_array=({input.all_plots})

        for sample_plots in ${{sample_array[@]}}
        do
        cp $sample_plots {params.plot_directory}
        done 2>{log}
        """
