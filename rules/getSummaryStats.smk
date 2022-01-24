
rule get_summary_stats:
    input:
        all_stats = expand(config["samples_dir"] + "/{sample_dir}/{sample}.stat", zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES)
    output:
        sum_stats = config["summary_dir"] + "/Stats.csv"        
    log:
        "logs/sumStats.log"
    shell:
        r"""
        echo "Sample,Total reads,Mapped reads,Mapped reads(>30nt), Mapped read %,Ref name,Ref Length,Ref Coverage,Ref coverage %,Average Depth,Min Depth,Max Depth" > {output.sum_stats};

        sample_array=({input.all_stats})

        for sample_stat in ${{sample_array[@]}}
        do
        sample=$(echo $sample_stat | rev | cut -d / -f 1 | rev | cut -d _ -f 1)
        echo -n "$sample,";

        cat $sample_stat |awk '{{for(i=1;i<=NF;i++) printf "%s, ",$i; print ""; }}'

        done >> {output.sum_stats} \
        2>{log}
        """


rule get_summary_amplicon_coverage:
    input:
        all_amp_cov = expand(config["samples_dir"] + "/{sample_dir}/{sample}.amplicon.cov", zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES)
    output:
        sum_amp_cov = config["summary_dir"] + "/ampliconDepth.csv"
    log:
        "logs/sumAmpCov.log"
    shell:
        r"""
        sample_array=({input.all_amp_cov})

        samples=()
        for sample_stat in ${{sample_array[@]}}
        do
        sample=$(echo $sample_stat | rev | cut -d / -f 1 | rev | cut -d _ -f 1)
        samples+=( $sample )
        done 2>{log}

        echo ${{samples[*]}} | sed 's/ /,/g' > {output.sum_amp_cov} 2>{log}

        paste {input.all_amp_cov} | awk '{{for(i=4;i<=NF;i+=4) printf "%s,",$i; print "";}}' >> {output.sum_amp_cov} 2>{log}
        """


rule collect_plots:
    input:
        amp_dep_plots = expand(config["samples_dir"] + "/{sample_dir}/{sample}-Amplicon-Depth.pdf", zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES),
        cov_plots = expand(config["samples_dir"] + "/{sample_dir}/{sample}-coverage.pdf", zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES)
    output:
        collected_amp_dep = expand(config["summary_dir"] + "/CoveragePlots/{sample}-Amplicon-Depth.pdf", zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES),
        collected_cov = expand(config["summary_dir"] + "/CoveragePlots/{sample}-coverage.pdf", zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES)
    params:
        plot_directory = config["summary_dir"] + "/CoveragePlots/"
    log:
        "logs/collectPlots.log"
    shell:
        r"""
        amp_dep_array=({input.amp_dep_plots})
        cov_array=({input.cov_plots})
        all_plots=(${{amp_dep_array[@]}} ${{cov_array[@]}})

        for sample_plots in ${{all_plots[@]}}
        do
        cp $sample_plots {params.plot_directory}
        done 2>{log}
        """


rule get_ambiguous_nucleotide_positions_and_N_counts:
    input:
        aligned_consensus = config["summary_dir"] + "/All-consensus_aligned.fa"
    output:
        ambig_nucs = config["summary_dir"] + "/ambig_nuc_pos.csv",
        N_counts = config["summary_dir"] + "/consensus_N_counts.csv"
    params:
        script = config["scripts"] + "/getAmbiguousPositions.py",
        ref_name = config["ref_genome"].split('/')[-1][:-3]
    log:
        "logs/getAmbPosAndNCounts.log"
    shell:
        r"""
        python {params.script} \
        {input.aligned_consensus} \
        {params.ref_name} \
        {output.ambig_nucs} \
        {output.N_counts} \
        2>{log}
        """


rule get_ambiguous_position_depth:
    input:
        ambig_nucs = rules.get_ambiguous_nucleotide_positions_and_N_counts.output.ambig_nucs,
        all_depth = expand(config["samples_dir"] + "/{sample_dir}/{sample}_Depth_trimx2.tsv", zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES)
    output:
        ambig_pos = temp(config["summary_dir"] + "/ambig_pos"),
        ambig_pos_dep = config["summary_dir"] + "/ambig_pos_depth.csv"
    params:
        sample_path_prefix = config["samples_dir"]
    log:
        "logs/getAmbPosDep.log"
    shell:
        r"""
        echo -e 'Sequence ID,Position,Depth' > {output.ambig_pos_dep}

        awk 'BEGIN {{FS=","}} {{OFS=","}} {{print $1,$3}}' {input.ambig_nucs} > {output.ambig_pos}

        while read sample
        do

            while read position
            do
                awk -v sample="$sample" -v pos="$position" '{{OFS=","}} $2 == pos {{print sample,pos,$3}}' {params.sample_path_prefix}/$sample/$sample*_Depth_trimx2.tsv >> {output.ambig_pos_dep}
            done < <(awk -v sample="$sample" '{{FS=","}} $1 == sample {{print $2}}' {output.ambig_pos})   #use list of ambiguous positions for given sample as input for loop

        done < <(awk '{{FS=","}}(NR>1) {{print $1}}' {input.ambig_nucs} | uniq) 2>{log}   #use unique list of samples as input for loop
        """


rule get_ambiguous_position_counts:
    input:
        ambig_nucs = rules.get_ambiguous_nucleotide_positions_and_N_counts.output.ambig_nucs
    output:
        ambig_pos_count = config["summary_dir"] + "/ambig_pos_count.csv"
    params:
        script = config["scripts"] + "/getAmbPosCounts.py"
    log:
        "logs/getAmbPosCount.log"
    shell:
        r"""
        python {params.script} \
        {input.ambig_nucs} \
        {output.ambig_pos_count} \
        2>{log}
        """


#Not sure this is necessary as info is aslready in Stats.csv
rule get_reference_coverage:
    input:
        all_bams = expand(config["samples_dir"] + "/{sample_dir}/{sample}_trimx2_sorted.bam", zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES),
        all_depth = expand(config["samples_dir"] + "/{sample_dir}/{sample}_Depth_trimx2.tsv", zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES)
    output:
        ref_cov = config["summary_dir"] + "/ref_cov.csv"
    log:
        "logs/getRefCov.log"
    shell:
        r"""
        echo -e 'Sequence ID,Ref Length,Ref Coverage,Ref Coverage %' > {output.ref_cov}

        bam_array=({input.all_bams})

        for sample_bam in ${{bam_array[@]}}
        do

            sample_dir=$(dirname $sample_bam)
            sample=$(echo $sample_dir | rev | cut -d / -f 1 | rev | cut -d _ -f 1)
            refLength=$(samtools idxstats $sample_bam | awk 'NR==1{{print $(NF-2);}}')
            refCoverage=$(awk '$NF< 10{{i++}};END{{print NR-i;}}' $sample_dir/Depth_trimx2.tsv)
            covPerc=$(echo $refCoverage $refLength | awk '{{printf "%0.2f", $1 / $2 * 100}}')

            echo $sample $refLength $refCoverage $covPerc | awk '{{OFS=","}} {{print $1,$2,$3+0,$4}}' >> {output.ref_cov}

        done 2>{log}
        """


rule get_masked_ambiguous_nucleotide_positions_and_N_counts:
    input:
        aligned_masked_consensus = config["summary_dir"] + "/All-masked-consensus_aligned.fa"
    output:
        masked_ambig_nucs = config["summary_dir"] + "/masked_ambig_nuc_pos.csv",
        masked_N_counts = config["summary_dir"] + "/masked_consensus_N_counts.csv"
    params:
        script = config["scripts"] + "/getAmbiguousPositions.py",
        ref_name = config["ref_genome"].split('/')[-1][:-3]
    log:
        "logs/getMaskedAmbPosAndNCounts.log"
    shell:
        r"""
        python {params.script} \
        {input.aligned_masked_consensus} \
        {params.ref_name} \
        {output.masked_ambig_nucs} \
        {output.masked_N_counts} \
        2>{log}
        """


rule assign_lineages:
    input:
        masked_consensus =  config["summary_dir"] + "/All-masked-consensus.fa"
    output:
        pangolin_report =  config["summary_dir"] + "/lineage_report.csv"
    params:
        threads = config["threads"],
        out_dir = config["summary_dir"]
    log:
        "logs/pangolin.log"
    shell:
        r"""
        pangolin \
        {input.masked_consensus} \
        -t {params.threads} \
        -o {params.out_dir} \
        2>{log}
        """
