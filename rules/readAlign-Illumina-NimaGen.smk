#bioinfo_pipe_version: v2
#The following changes were made to address the issue where a high number of ambiguous nucleotides were>#trim_galore -q20 was changed to --nextseq30
#prinseq added with -trim_qual_right 30
# -x 2 parameter added to ivar trim
#Further downstream, ambiguous nucleotides in amplicon regions with particularly high number of ambigui>#This is done in /home4/sd261/bin/postReadAlign, which is called from autoNimagen.sh
#These steps were taken for batches 187, 189 and 190, where the ambiguity issue was observed.

#fq1=$(ls *_R1*.fastq)
#fq2=$(ls *_R2*.fastq)
#
#sample=$(ls *_R1*.fastq|cut -d_ -f1);
#
#trim_galore --nextseq 30 --dont_gzip --length 50 --paired $fq1 $fq2 > /dev/null  2>&1
#
#trimmed1=$(ls *R1*val_1.fq)
#trimmed2=$(ls *R2*val_2.fq)
#
#prinseq-lite.pl -trim_qual_right 30 -fastq $trimmed1 -fastq2 $trimmed2 -out_good ${sample}_prinseq -out_bad ${sample}_bad -min_len 50
#
#bwa mem -t 10 /home3/bvv2t/bin/SARC-CoV-2-Scripts/Ref/MN908947.3.fa ${sample}_prinseq_1.fastq ${sample}_prinseq_2.fastq | samtools view -F4 -bS - | samtools sort -@10 -o $sample.bam -
#samtools index ${sample}.bam
#
#ivar trim -i ${sample}.bam -p ${sample}_trimx2 -b /home3/bvv2t/bin/SARC-CoV-2-Scripts/Primers/nCoV-2019-NimaGen-$1.bed -m 30 -x 2
#samtools sort -@10 ${sample}_trimx2.bam -o ${sample}_trim_sort.bam
#mv ${sample}_trim_sort.bam ${sample}_trimx2.bam
#samtools index ${sample}_trimx2.bam
#
#samtools mpileup -B -A -d 10000000 -Q 0 ${sample}_trimx2.bam |ivar consensus -p ${sample}_trimx2_ivar_consensus -n N -t 0.6 -m 10
##perl ~/bin/Sequence-manipulation/CountAmbig.pl -in ${sample}_trimx2_ivar_consensus.fa -detail > Ambiguous_trimx2.tsv
##samtools depth -a -d 0 -o Depth_trimx2.tsv ${sample}_trimx2.bam


#GLOBAL WILDCARDS, not currently in use...
SAMPLE_DIRS = [x.split('_')[0] for x in glob_wildcards(config["input_path"] + "/{sample}_R1_001.fastq.gz").sample]
SAMPLES = glob_wildcards(config["input_path"] + "/{sample}_R1_001.fastq.gz").sample


rule trim_reads:
    input:
        fq1 = config["input_path"] + "/{sample}_R1_001.fastq.gz",
        fq2 = config["input_path"] + "/{sample}_R2_001.fastq.gz"
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
    log:
        "logs/prinseq/{sample}.log"
    params:
        out_path = config["output_path"] + "/{sample}/"
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
