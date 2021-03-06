
#May need to use the --force flag.
rule demultiplex_samples:
    input:
        run_dir = config["run_directory"]
    output:
        fastqs = expand(config["fastq_directory"] + "/{sample}_R{pair}_001.fastq.gz", sample=SAMPLES, pair=[1,2])
    params:
        fastq_dir = config["fastq_directory"],
        sample_sheet = config["sample_sheet"],
        bcl_threads = config["bcl-convert_threads"]
    log:
        "logs/bclConvert.log"
    shell:
        r"""
        bcl-convert \
        --bcl-input-directory {input.run_dir} \
        --output-directory {params.fastq_dir} \
        --sample-sheet {params.sample_sheet} \
        --no-lane-splitting true \
        --bcl-num-conversion-threads {params.bcl_threads} \
        2>{log}

        chmod -R 774 {params.fastq_dir} 2>>{log}   #To prevent overwrite by users without group permissions.
        """


#Input determined by lambda function, where sample_dir wildcard from output is used as a key for the 
#SAMPLES_DICT dictionary defined in the Snakefile.
#Trimmed output are marked as temporary and are removed once the pipeline is finished executing.
#This saves on space.
#--nextseq parameter used to account for 2-colour chemistry
rule trim_reads:
    input:
        fq1 = lambda wildcards: config["fastq_directory"] + "/" + SAMPLES_DICT[wildcards.sample_dir] + "_R1_001.fastq.gz",
        fq2 = lambda wildcards: config["fastq_directory"] + "/" + SAMPLES_DICT[wildcards.sample_dir] + "_R2_001.fastq.gz"
    output:
        trimmed1 = temp(config["samples_dir"] + "/{sample_dir}/{sample}_val_1.fq"),
        trimmed2 = temp(config["samples_dir"] + "/{sample_dir}/{sample}_val_2.fq")
    params:
        out_dir = config["samples_dir"] + "/{sample_dir}"
    log:
        "logs/trimGalore/{sample_dir}/{sample}.log"
    shell:
        r"""
        trim_galore \
        --nextseq 30 \
        --dont_gzip \
        --length 50 \
        -o {params.out_dir} \
        --basename {wildcards.sample} \
        --paired {input.fq1} {input.fq2} \
        > /dev/null  2>{log}
        """


#prinseq is used to remove low quality bases from 3' end.
rule prinseq_trim:
    input:
        trimmed1 = lambda wildcards: config["samples_dir"] + "/{sample_dir}/" + SAMPLES_DICT[wildcards.sample_dir] + "_val_1.fq",
        trimmed2 = lambda wildcards: config["samples_dir"] + "/{sample_dir}/" + SAMPLES_DICT[wildcards.sample_dir] + "_val_2.fq"
    output:
        prin_trim1 = config["samples_dir"] + "/{sample_dir}/{sample_dir}_prinseq_1.fastq",
        prin_trim2 = config["samples_dir"] + "/{sample_dir}/{sample_dir}_prinseq_2.fastq"
    params:
        threads = config["threads"],
        out_path = config["samples_dir"] + "/{sample_dir}/"
    log:
        "logs/prinseq/{sample_dir}.log"
    shell:
        r"""
        prinseq++ \
        -trim_qual_right 30 \
        -fastq {input.trimmed1} \
        -fastq2 {input.trimmed2} \
        -threads {params.threads} \
        -out_good {output.prin_trim1} \
        -out_good2 {output.prin_trim2} \
        -out_bad {params.out_path}{wildcards.sample_dir}_bad_1.fastq \
        -out_bad2 {params.out_path}{wildcards.sample_dir}_bad_2.fastq \
        -out_single {params.out_path}{wildcards.sample_dir}_single_1.fastq \
        -out_single2 {params.out_path}{wildcards.sample_dir}_single_2.fastq \
        -min_len 50 \
        2>{log}
        """


#Reads are algined to reference and sorted.
rule align_reads:
    input:
        prin_trim1 = rules.prinseq_trim.output.prin_trim1,
        prin_trim2 = rules.prinseq_trim.output.prin_trim2
    output:
        bam = config["samples_dir"] + "/{sample_dir}/{sample_dir}.bam"
    params:
        ref = config["ref_genome"],
        threads = config["threads"]
    log:
        "logs/bwamem/{sample_dir}.log"
    shell:
        r"""
        bwa mem \
        -t {params.threads} \
        {params.ref} \
        {input.prin_trim1} \
        {input.prin_trim2} \
        2>{log} \
	| samtools view -F4 -bS - \
        2>>{log} \
        | samtools sort \
        -@{params.threads} \
        -o {output.bam} -
        2>>{log}
        """


rule bam_index:
    input:
        bam = rules.align_reads.output.bam
    output:
        index = config["samples_dir"] + "/{sample_dir}/{sample_dir}.bam.bai"
    log:
        "logs/samindex/{sample_dir}.log"
    shell:
        r"""
        samtools index \
        {input.bam} \
        {output.index} \
        2>{log}
        """


#ivar used to trim primers
#-x parameter is used to define the offset of a read relative to the primer.
rule ivar_trim:
    input:
        bam = rules.align_reads.output.bam,
        index = rules.bam_index.output.index
    output:
        ivar_trimmed = config["samples_dir"] + "/{sample_dir}/{sample_dir}_trimx2.bam"
    params:
        out_path = config["samples_dir"] + "/{sample_dir}/{sample_dir}_trimx2",
        primer = config["primer_bed"]
    log:
        "logs/ivartrim/{sample_dir}.log"
    shell:
        r"""
        ivar trim \
        -i {input.bam} \
        -p {params.out_path} \
        -b {params.primer} \
        -m 30 \
        -x 2 \
        2>{log}
        """


rule sort_ivar_trimmed:
    input:
        ivar_trimmed = rules.ivar_trim.output.ivar_trimmed
    output:
        ivar_trim_sorted = config["samples_dir"] + "/{sample_dir}/{sample_dir}_trimx2_sorted.bam"
    params:
        threads = config["threads"]
    log:
        "logs/samsortIvar/{sample_dir}.log"
    shell:
        r"""
        samtools sort \
        -@{params.threads} \
        {input.ivar_trimmed} \
        -o {output.ivar_trim_sorted} \
        2>{log}
        """


rule index_ivar_trimmed:
    input:
        bam = rules.sort_ivar_trimmed.output.ivar_trim_sorted
    output:
        index = config["samples_dir"] + "/{sample_dir}/{sample_dir}_trimx2.bam.bai"
    log:
        "logs/samindexIvar/{sample_dir}.log"
    shell:
        r"""
        samtools index \
        {input.bam} \
        {output.index} \
        2>{log}
        """


#Consensus fasta is created from ivar trimmed reads.
rule generate_consensus:
    input:
        ivar_bam = rules.sort_ivar_trimmed.output.ivar_trim_sorted,
        ivar_bam_index = rules.index_ivar_trimmed.output.index
    output:
        consensus = config["samples_dir"] + "/{sample_dir}/{sample_dir}_trimx2_ivar_consensus.fa"
    params:
        out_path = config["samples_dir"] + "/{sample_dir}/{sample_dir}_trimx2_ivar_consensus"
    log:
        "logs/ivarConsensus/{sample_dir}.log"
    shell:
        r"""
        samtools mpileup \
        -B \
        -A \
        -d 10000000 \
        -Q 0 \
        {input.ivar_bam} \
        |ivar consensus \
        -p {params.out_path} \
        -n N \
        -t 0.6 \
        -m 10 \
        2>{log}
        """
