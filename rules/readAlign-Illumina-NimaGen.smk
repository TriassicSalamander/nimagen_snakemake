
##### Get Sample List #####
sample_sheet = open(config["sample_sheet"], 'r')
is_sample_id = False
SAMPLE_IDS = []
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
        SAMPLE_IDS.append(row.split(',')[0] + '_S' + str(sample_count))
sample_sheet.close()


#GLOBAL WILDCARDS, not currently in use...
SAMPLE_DIRS = [x.split('_')[0] for x in glob_wildcards(config["fastq_directory"] + "/{sample}_R1_001.fastq.gz").sample]
SAMPLES = glob_wildcards(config["fastq_directory"] + "/{sample}_R1_001.fastq.gz").sample


#Flag --first-tile-only being used for testing, REMEMBER TO REMOVE.
#bcl-convert is complaining about the output directory existing, though I can't see it.
#Using --force to get around this for now.
rule demultiplex_samples:
    input:
        run_dir = config["run_directory"]
    output:
        fastqs = expand(config["fastq_directory"] + "/{sample}_R{pair}_001.fastq.gz", sample=SAMPLE_IDS, pair=[1,2])
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
        --force \
        --first-tile-only true \
        --bcl-num-conversion-threads {params.bcl_threads} \
        2>{log}

        chmod -R 774 {params.fastq_dir} 2>{log}
        """


rule trim_reads:
    input:
        fq1 = config["fastq_directory"] + "/{sample}_R1_001.fastq.gz",
        fq2 = config["fastq_directory"] + "/{sample}_R2_001.fastq.gz"
#        fq1 = expand(config["input_path"] + "/{sample}_R1_001.fastq.gz", sample=SAMPLES),
#        fq2 = expand(config["input_path"] + "/{sample}_R2_001.fastq.gz", sample=SAMPLES)
    output:
        trimmed1 = config["output_path"] + "/{sample}/{sample}_val_1.fq",
        trimmed2 = config["output_path"] + "/{sample}/{sample}_val_2.fq"
#        trimmed1 = expand(config["output_path"] + "/{sample_dir}/{sample}_val_1.fq", zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES),
#        trimmed2 = expand(config["output_path"] + "/{sample_dir}/{sample}_val_2.fq", zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES)
    params:
        out_path = config["output_path"] + "/{sample}/",
#        out_path = expand(config["output_path"] + "/{sample_dir}/", sample_dir=SAMPLE_DIRS),
#        out_base = expand("{sample}", sample=SAMPLES)
    log:
        "logs/trim_galore/{sample}.log"
    shell:
        r"""
        trim_galore \
        --nextseq 30 \
        --dont_gzip \
        --length 50 \
        -o {params.out_path} \
        --basename {wildcards.sample} \
        --paired {input.fq1} {input.fq2} \
        > /dev/null  2>{log}
        """


rule prinseq_trim:
    input:
        trimmed1 = rules.trim_reads.output.trimmed1,
        trimmed2 = rules.trim_reads.output.trimmed2
    output:
        prin_trim1 = config["output_path"] + "/{sample}/{sample}_prinseq_1.fastq",
        prin_trim2 = config["output_path"] + "/{sample}/{sample}_prinseq_2.fastq"
    params:
        out_path = config["output_path"] + "/{sample}/"
    log:
        "logs/prinseq/{sample}.log"
    shell:
        r"""
        prinseq-lite.pl \
        -trim_qual_right 30 \
        -fastq {input.trimmed1} \
        -fastq2 {input.trimmed2} \
        -out_good {params.out_path}{wildcards.sample}_prinseq \
        -out_bad {params.out_path}{wildcards.sample}_bad \
        -min_len 50 \
        -log {log}
        """


rule align_reads:
    input:
        prin_trim1 = rules.prinseq_trim.output.prin_trim1,
        prin_trim2 = rules.prinseq_trim.output.prin_trim2
    output:
        bam = config["output_path"] + "/{sample}/{sample}.bam"
    params:
        ref = config["ref_genome"],
        threads = config["threads"],
        out_path = config["output_path"] + "/{sample}/"
    log:
        "logs/bwamem/{sample}.log"
    shell:
        r"""
        bwa mem \
        -t {params.threads} \
        {params.ref} \
        {input.prin_trim1} \
        {input.prin_trim2} \
        2>{log} \
	| samtools view -F4 -bS - \
        2>{log}
        | samtools sort \
        -@{params.threads} \
        -o {output.bam} -
        2>{log}
        """


rule bam_index:
    input:
        bam = rules.align_reads.output.bam
    output:
        index = config["output_path"] + "/{sample}/{sample}.bam.bai"
    log:
        "logs/samindex/{sample}.log"
    shell:
        r"""
        samtools index \
        {input.bam} \
        {output.index} \
        2>{log}
        """


rule ivar_trim:
    input:
        bam = rules.align_reads.output.bam,
        index = rules.bam_index.output.index
    output:
        ivar_trimmed = config["output_path"] + "/{sample}/{sample}_trimx2.bam"
    params:
        out_path = config["output_path"] + "/{sample}/{sample}_trimx2",
        primer = config["primer_bed"]
    log:
        "logs/ivartrim/{sample}.log"
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
        ivar_trim_sorted = config["output_path"] + "/{sample}/{sample}_trimx2_sorted.bam"
    params:
        threads = config["threads"]
    log:
        "logs/samsort_ivar/{sample}.log"
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
        index = config["output_path"] + "/{sample}/{sample}_trimx2.bam.bai"
    log:
        "logs/samindex_ivar/{sample}.log"
    shell:
        r"""
        samtools index \
        {input.bam} \
        {output.index} \
        2>{log}
        """


rule generate_consensus:
    input:
        ivar_bam = rules.sort_ivar_trimmed.output.ivar_trim_sorted,
        ivar_bam_index = rules.index_ivar_trimmed.output.index
    output:
        consensus = config["output_path"] + "/{sample}/{sample}_trimx2_ivar_consensus.fa"
    params:
        out_path = config["output_path"] + "/{sample}/{sample}_trimx2_ivar_consensus"
    log:
        "logs/ivar_consensus/{sample}.log"
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
