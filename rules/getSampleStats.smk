
rule get_depth:
    input:
        bam = rules.sort_ivar_trimmed.output.ivar_trim_sorted
    output:
        depth = config["samples_dir"] + "/{sample_dir}/{sample_dir}_Depth_trimx2.tsv"
    log:
        "logs/samDepth/{sample_dir}.log"
    shell:
        r"""
        samtools depth \
        -a \
        -d 0 \
        -o {output.depth} \
        {input.bam} \
        2>{log}
        """

#Creates variant and depth files required for demixing samples with freyja.
#Could maybe save time by using depth file already created by pipeline.
rule get_freyja_variants_and_depth:
    input:
        bam = rules.sort_ivar_trimmed.output.ivar_trim_sorted
    output:
        variants = config["samples_dir"] + "/{sample_dir}/{sample_dir}_freyja_variants.tsv",
        depth = config["samples_dir"] + "/{sample_dir}/{sample_dir}_freyja_depth.tsv"
    params:
        ref = config["ref_genome"]
    log:
        "logs/freyjaVariants/{sample_dir}.log"
    shell:
        r"""
        freyja variants \
        {input.bam} \
        --variants {output.variants} \
        --depths {output.depth} \
        --ref {params.ref} \
        2>{log}
        """


#Demixes samples to get lineage abundances.
#Should be added as an option to config as only some runs will have mixed samples.
rule get_freyja_demixed:
    input:
        variants = rules.get_freyja_variants_and_depth.output.variants,
        depth = rules.get_freyja_variants_and_depth.output.depth
    output:
        demixed = config["samples_dir"] + "/{sample_dir}/{sample_dir}_freyja_demixed.tsv"
    log:
        "logs/freyjaDemix/{sample_dir}.log"
    shell:
        r"""
        #Only demix if depth file isn't empty. Freyja throws an error otherwise.
        #If it is empty, create a blank demixed output.
        if [ -s {input.depth} ]
        then
            freyja demix \
            {input.variants} \
            {input.depth} \
            --output {output.demixed} \
            2>{log}
        else
            printf "\t{input.variants}\nsummarized\nlineages\nabundances\nresid" > {output.demixed}
        fi
        """


#Calls script which determines the depth of each amplicon.
#Output is a tsv file with columns representing: amplicon number, start position, end position and depth.
rule get_amplicon_depth:
    input:
        bam = rules.sort_ivar_trimmed.output.ivar_trim_sorted
    output:
        amp_dep = config["samples_dir"] + "/{sample_dir}/{sample_dir}.amplicon.cov"
    params:
        amp_dep_script = config["scripts"] + "/ampliconDepth-NimaGen.sh",
        PE_primer = config["PE_primer_bed"]
    log:
        "logs/SampleAmpDepth/{sample_dir}.log"
    shell:
        r"""
        {params.amp_dep_script} \
        {params.PE_primer} \
        {input.bam} \
        > {output.amp_dep} \
        2>{log}
        """


#Calls script which determines a number of statistics.
#Output is a space-delimited file with columns representing: Total reads,Mapped reads,Mapped reads(>30nt), Mapped read %,Ref name,Ref Length,Ref Coverage,Ref coverage %,Average Depth,Min Depth,Max Depth
rule get_stats:
    input:
        trimmed1 = lambda wildcards: config["samples_dir"] + "/{sample_dir}/" + SAMPLES_DICT[wildcards.sample_dir] + "_val_1.fq",
        trimmed2 = lambda wildcards: config["samples_dir"] + "/{sample_dir}/" + SAMPLES_DICT[wildcards.sample_dir] + "_val_2.fq",
        bam = rules.sort_ivar_trimmed.output.ivar_trim_sorted
    output:
        stats = config["samples_dir"] + "/{sample_dir}/{sample_dir}.stat"
    params:
        stats_script = config["scripts"] + "/getAlignmentStats"
    log:
        "logs/SampleStats/{sample_dir}.log"
    shell:
        r"""
        {params.stats_script} \
        -1 {input.trimmed1} \
        -2 {input.trimmed2} \
        -b {input.bam} \
        > {output.stats} \
        2>{log}
        """


rule get_coverage_plot:
    input:
        amp_dep = rules.get_amplicon_depth.output.amp_dep
    output:
        cov_plot = config["samples_dir"] + "/{sample_dir}/{sample_dir}-coverage.pdf"
    params:
        cov_plot_script = config["scripts"] + "/getCoverage.R",
        out_prefix = config["samples_dir"] + "/{sample_dir}/{sample_dir}"
    log:
        "logs/covplot/{sample_dir}.log"
    shell:
        r"""
        Rscript {params.cov_plot_script} \
        {input.amp_dep} \
        {params.out_prefix} \
        2>{log}
        """


rule get_amplicon_depth_plot:
    input:
        amp_dep = rules.get_amplicon_depth.output.amp_dep
    output:
        dep_plot = config["samples_dir"] + "/{sample_dir}/{sample_dir}-Amplicon-Depth.pdf"
    params:
        dep_plot_script = config["scripts"] + "/getAmpCoverageNG.R",
        out_prefix = config["samples_dir"] + "/{sample_dir}/{sample_dir}"
    log:
        "logs/ampDepPlot/{sample_dir}.log"
    shell:
        r"""
        Rscript {params.dep_plot_script} \
        {input.amp_dep} \
        {params.out_prefix} \
        2>{log}
        """
