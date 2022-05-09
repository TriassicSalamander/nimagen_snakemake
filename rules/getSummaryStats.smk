
#Rule which collects the stat files for each sample.
rule get_summary_stats:
    input:
        all_stats = expand(config["samples_dir"] + "/{sample_dir}/{sample_dir}.stat", sample_dir=SAMPLE_DIRS)
    output:
        sum_stats = temp(config["summary_dir"] + "/Stats_temp.csv")
    log:
        "logs/sumStats.log"
    shell:
        r"""
        #Create column headers.
        echo "Sample,Total reads,Mapped reads,Mapped reads(>30nt), Mapped read %,Ref name,Ref Length,Ref Coverage,Ref coverage %,Average Depth,Min Depth,Max Depth" > {output.sum_stats};

        sample_array=({input.all_stats})   #Create array of all the stats files

        IFS=$'\n'
        sorted_samples=($(sort <<<"${{sample_array[*]}}"))   #These three lines are for sorting the array.
        unset IFS

        for sample_stat in ${{sorted_samples[@]}}   #Loop through the sorted array
        do
        sample=$(echo $sample_stat | rev | cut -d / -f 1 | rev | cut -d . -f 1)   #Get sample name from the file path. This will be used for the 'Sample' column.
        echo -n "$sample,";

        cat $sample_stat |awk '{{for(i=1;i<=NF;i++) printf "%s, ",$i; print ""; }}'   #Print the stat file to output as csv.

        done >> {output.sum_stats} \
        2>{log}
        """


#Rule which collects the amplicon depth files for each sample.
rule get_summary_amplicon_coverage:
    input:
        all_amp_cov = expand(config["samples_dir"] + "/{sample_dir}/{sample_dir}.amplicon.cov", sample_dir=SAMPLE_DIRS)
    output:
        sum_amp_cov = config["summary_dir"] + "/ampliconDepth.csv"
    log:
        "logs/sumAmpCov.log"
    shell:
        r"""
        sample_array=({input.all_amp_cov})

        IFS=$'\n'
        sorted_samples=($(sort <<<"${{sample_array[*]}}"))
        unset IFS

        samples=()   #Create empty array which will be populated with sample names.
        for sample_cov in ${{sorted_samples[@]}}
        do
        sample=$(echo $sample_cov | rev | cut -d / -f 1 | rev | cut -d . -f 1)
        samples+=( $sample )
        done 2>{log}

        echo ${{samples[*]}} | sed 's/ /,/g' > {output.sum_amp_cov} 2>{log}   #Print sample names as csv to output. These will be the headers.

        paste ${{sorted_samples[@]}} | awk '{{for(i=4;i<=NF;i+=4) printf "%s,",$i; print "";}}' >> {output.sum_amp_cov} 2>{log}   #From all the amplicon depth files, take the depth column and print to output.
        """


#Rule for copying all of the sample pdfs to the CoveragePlots directory.
rule collect_plots:
    input:
        amp_dep_plots = expand(config["samples_dir"] + "/{sample_dir}/{sample_dir}-Amplicon-Depth.pdf", sample_dir=SAMPLE_DIRS),
        cov_plots = expand(config["samples_dir"] + "/{sample_dir}/{sample_dir}-coverage.pdf", sample_dir=SAMPLE_DIRS)
    output:
        collected_amp_dep = expand(config["summary_dir"] + "/CoveragePlots/{sample_dir}-Amplicon-Depth.pdf", sample_dir=SAMPLE_DIRS),
        collected_cov = expand(config["summary_dir"] + "/CoveragePlots/{sample_dir}-coverage.pdf", sample_dir=SAMPLE_DIRS)
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


#Rule which calls a script to determine ambiguous base and position for samples.
#The script also outputs the number of 'N's in each sample consensus.
rule get_ambiguous_nucleotide_positions_and_N_counts:
    input:
        aligned_consensus = config["summary_dir"] + "/All-consensus_aligned.fa"
    output:
        ambig_nucs = temp(config["summary_dir"] + ("/ambig_nuc_pos.csv")),
        N_counts = temp(config["summary_dir"] + "/consensus_N_counts.csv")
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


#Rule which determines the depth and position for ambiguous nucleotides.
rule get_ambiguous_position_depth:
    input:
        ambig_nucs = rules.get_ambiguous_nucleotide_positions_and_N_counts.output.ambig_nucs,
        all_depth = expand(config["samples_dir"] + "/{sample_dir}/{sample_dir}_Depth_trimx2.tsv", sample_dir=SAMPLE_DIRS)
    output:
        ambig_pos = temp(config["summary_dir"] + "/ambig_pos"),
        ambig_pos_dep = temp(config["summary_dir"] + "/ambig_pos_depth.csv")
    params:
        sample_path_prefix = config["samples_dir"]
    log:
        "logs/getAmbPosDep.log"
    shell:
        r"""
        #Print columnn headers to output.
        echo -e 'Sequence ID,Position,Depth' > {output.ambig_pos_dep}

        awk 'BEGIN {{FS=","}} {{OFS=","}} {{print $1,$3}}' {input.ambig_nucs} > {output.ambig_pos}   #Print sample name and position columns to temp output.

        while read sample   #Loop through sample names.
        do

            while read position   #Loop through ambiguous positions.
            do
                awk -v sample="$sample" -v pos="$position" '{{OFS=","}} $2 == pos {{print sample,pos,$3}}' {params.sample_path_prefix}/$sample/$sample*_Depth_trimx2.tsv >> {output.ambig_pos_dep}   #Get depth at position for given sample and print to output.
            done < <(awk -v sample="$sample" '{{FS=","}} $1 == sample {{print $2}}' {output.ambig_pos})   #use list of ambiguous positions for given sample as input for loop

        done < <(awk '{{FS=","}}(NR>1) {{print $1}}' {input.ambig_nucs} | uniq) 2>{log}   #use unique list of samples as input for loop
        """


#Rule which combines outputs from get_ambiguous_nucleotide_positions_and_N_counts and get_ambiguous_position_depth
#to give a csv file with columns: Sequence ID, Position, Ambiguous Base and Depth
rule combine_ambig_nuc_pos_and_depth:
    input:
        ambig_pos = rules.get_ambiguous_nucleotide_positions_and_N_counts.output.ambig_nucs,
        ambig_dep = rules.get_ambiguous_position_depth.output.ambig_pos_dep
    output:
        ambig_pos_and_dep = config["summary_dir"] + "/ambig_nuc_pos_and_dep.csv"
    log:
        "logs/combineAmbigPosDep.log"
    run:
        import pandas as pd

        pos_df = pd.read_csv(input.ambig_pos)
        dep_df = pd.read_csv(input.ambig_dep)
        combined_df = pd.merge(pos_df,dep_df)
        combined_df.to_csv(output.ambig_pos_and_dep, index=False)


#Rule which calls script which counts the number of ambiguos nucleotides per sample.
rule get_ambiguous_position_counts:
    input:
        ambig_nucs = rules.get_ambiguous_nucleotide_positions_and_N_counts.output.ambig_nucs
    output:
        ambig_pos_count = temp(config["summary_dir"] + "/ambig_pos_count.csv")
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


#Rule which adds ambiguous nucleotide count and N counts to the summary stats file.
rule add_N_and_ambig_counts_to_stats:
    input:
        stats = rules.get_summary_stats.output.sum_stats,
        N_counts = rules.get_ambiguous_nucleotide_positions_and_N_counts.output.N_counts,
        ambig_counts = rules.get_ambiguous_position_counts.output.ambig_pos_count
    output:
        stats = config["summary_dir"] + "/Stats.csv"
    run:
        import pandas as pd

        stats_df = pd.read_csv(input.stats, index_col=False)

        N_count_df = pd.read_csv(input.N_counts, index_col=False)
        N_count_df.rename(columns={'Sequence ID':'Sample', 'N Count':'Consensus N Count'}, inplace=True)   #Rename columns to match summary stats file.

        ambig_count_df = pd.read_csv(input.ambig_counts, index_col=False)
        ambig_count_df.rename(columns={'Sequence ID':'Sample', 'Amb Count':'Ambiguous Nucleotide Count'}, inplace=True)

        stats_N_count = pd.merge(stats_df, N_count_df, how='left', on='Sample')
        stats_N_ambig_count = pd.merge(stats_N_count, ambig_count_df, how='left', on='Sample')

        stats_N_ambig_count.fillna(0, downcast='infer', inplace=True)   #Samples without 'N's or ambiguous nucleotides will have n/a values, so these need to be replaced with '0's.

        stats_N_ambig_count.to_csv(output.stats, index=False)


#Rule which creates a pangolin lineage report for the masked consensus files.
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
        conda activate pangolin   #Move into the pangolin environment, which is separate from the nimagen_snakemake environment.

        pangolin --update   #Make sure pangolin is up to date

        pangolin \
        {input.masked_consensus} \
        -t {params.threads} \
        -o {params.out_dir} \
        2>{log}

        conda deactivate   #Exit the environment, returning to nimagen_snakemake.
        """


rule collect_freyja_demixed:
    input:
        all_demixed = expand(config["samples_dir"] + "/{sample_dir}/{sample_dir}_freyja_demixed.tsv", sample_dir=SAMPLE_DIRS)
    output:
        collected_demixed = expand(config["summary_dir"] + "/FreyjaDemixed/{sample_dir}_freyja_demixed.tsv", sample_dir=SAMPLE_DIRS)
    params:
        demixed_dir = config["summary_dir"] + "/FreyjaDemixed/"
    log:
        "logs/collectFreyjaDemixed.log"
    shell:
        r"""
        demixed_array=({input.all_demixed})

        for demixed in ${{demixed_array[@]}}
        do
        cp $demixed {params.demixed_dir}
        done 2>{log}
        """


rule aggregate_freyja_demixed:
    input:
        collected_demixed = rules.collect_freyja_demixed.output.collected_demixed
    output:
        aggregated_demmixed = config["summary_dir"] + "/All-freyja-demixed.tsv"
    params:
        demixed_dir = config["summary_dir"] + "/FreyjaDemixed/"
    log:
        "logs/aggregateFreyja.log"
    shell:
        r"""
        freyja aggregate \
        {params.demixed_dir} \
        --output {output.aggregated_demmixed} \
        2>{log}
        """


rule plot_aggregate_freyja:
    input:
        aggregated_demixed = rules.aggregate_freyja_demixed.output.aggregated_demmixed
    output:
        temp_aggregated = config["summary_dir"] + "/All-freyja-demixed_temp.tsv",
        freyja_summarised_plot = config["summary_dir"] + "/All-freyja-summarised-plot.pdf",
        freyja_lineage_plot = config["summary_dir"] + "/All-freyja-lineage-plot.pdf"
    log:
        "logs/freyjaPlot.log"
    run:
        import pandas as pd

        #First filter input to only retain non-blank demixed samples.
        demixed_df = pd.read_csv(input.aggregated_demixed, sep="\t", index_col=0)
        demixed_df.dropna(inplace=True)
        demixed_df.to_csv(output.temp_aggregated, sep='\t')

        shell(r"""
            freyja plot \
            {output.temp_aggregated} \
            --output {output.freyja_summarised_plot} \
            2>{log}

            freyja plot \
            {output.temp_aggregated} \
            --lineages \
            --output {output.freyja_lineage_plot} \
            2>{log}
        """)
