
rule get_depth:
    input:
        bam = rules.sort_ivar_trimmed.output.ivar_trim_sorted
    output:
        depth = config["samples_dir"] + "/{sample_dir}/{sample_dir}_Depth_trimx2.tsv"
    log:
        "logs/samdepth/{sample_dir}.log"
    shell:
        r"""
        samtools depth \
        -a \
        -d 0 \
        -o {output.depth} \
        {input.bam} \
        2>{log}
        """


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
