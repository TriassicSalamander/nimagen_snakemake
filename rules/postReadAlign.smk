
rule aggregate_consensus:
    input:
        consensus_files = expand(config["samples_dir"] + "/{sample_dir}/{sample}_trimx2_ivar_consensus.fa",zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES)
    output:
        aggregated_consensus = config["summary_dir"] + "/All-consensus.fa"
    log:
        "logs/aggregateConsensus.log"
    shell:
        r"""
        sample_array=({input.consensus_files})

        for sample_consensus in ${{sample_array[@]}}
        do
        fasta_header=$(echo $sample_consensus | rev | cut -d / -f 1 | rev | cut -d _ -f1)
        awk -v fasta_header="$fasta_header" '{{if($1~/>/){{print ">"fasta_header;}} else print}}' $sample_consensus
        done >{output.aggregated_consensus} \
        2>{log}
        """


rule align_consensus:
    input:
        aggregated_consensus = rules.aggregate_consensus.output.aggregated_consensus
    output:
        consensus_with_ref = temp(config["summary_dir"] + "/All-consensus_with_ref.fa"),
        aligned_consensus = config["summary_dir"] + "/All-consensus_aligned.fa"
    params:
        ref = config["ref_genome"]
    log:
        "logs/alignConsensus.log"
    shell:
        r"""
        cp {input.aggregated_consensus} {output.consensus_with_ref} 2>{log}

        cat {params.ref} >> {output.consensus_with_ref} 2>{log}

        mafft {output.consensus_with_ref} > {output.aligned_consensus} 2>{log}
        """


rule align_individual_consensus:
    input:
        sample_consensus = config["samples_dir"] + "/{sample_dir}/{sample}_trimx2_ivar_consensus.fa"
    output:
        consensus_with_ref = temp(config["samples_dir"] + "/{sample_dir}/{sample}_consensus_with_ref.fa"),
        aligned_consensus = config["samples_dir"] + "/{sample_dir}/{sample}_consensus_with_ref_aligned.fa"
    params:
        ref = config["ref_genome"]
    log:
        "logs/alignIndConsensus/{sample_dir}/{sample}.log"
    shell:
        r"""
        cp {input.sample_consensus} {output.consensus_with_ref} 2>{log}

        cat {params.ref} >> {output.consensus_with_ref} 2>{log}

        mafft {output.consensus_with_ref} > {output.aligned_consensus} 2>{log}
        """


rule mask_ambiguous_nucleotides:
    input:
        aligned_ind_consensus = rules.align_individual_consensus.output.aligned_consensus
    output:
        masked_consensus = config["samples_dir"] + "/{sample_dir}/{sample}_aligned_masked_consensus.fa"
    params:
        script = config["scripts"] + "/maskAmbNucs.py",
        mask_regions = config["ambig_regions"]
    log:
        "logs/maskAmbigNucs/{sample_dir}/{sample}.log"
    shell:
        r"""
        python {params.script} \
        {params.mask_regions} \
        {input.aligned_ind_consensus} \
        {output.masked_consensus} \
        2>{log}
        """


rule unalign_masked_consensus:
    input:
        masked_consensus = rules.mask_ambiguous_nucleotides.output.masked_consensus
    output:
        unaligned_consensus = config["samples_dir"] + "/{sample_dir}/{sample}_masked_consensus.fa"
    params:
        script = config["scripts"] + "/removeRefAndGaps.py",
        ref_name = config["ref_genome"].split('/')[-1][:-3]
    log:
       "logs/removeRefandGaps/{sample_dir}/{sample}.log"
    shell:
        r"""
        python {params.script} \
        {input.masked_consensus} \
        {params.ref_name} \
        {output.unaligned_consensus} \
        2>{log}
        """


rule aggregate_masked_consensus:
    input:
        masked_consensus_files = expand(config["samples_dir"] + "/{sample_dir}/{sample}_masked_consensus.fa", zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES)
    output:
        aggregated_masked_consensus = config["summary_dir"] + "/All-masked-consensus.fa"
    log:
        "logs/aggregateMaskedConsensus.log"
    shell:
        r"""
        sample_array=({input.masked_consensus_files})

        for sample_consensus in ${{sample_array[@]}}
        do
        fasta_header=$(echo $sample_consensus | rev | cut -d / -f 1 | rev | cut -d _ -f 1)
        awk -v fasta_header="$fasta_header" '{{if($1~/>/){{print ">"fasta_header;}} else print}}' $sample_consensus
        done >{output.aggregated_masked_consensus} \
        2>{log}
        """


rule align_masked_consensus:
    input:
        aggregated_masked_consensus = rules.aggregate_masked_consensus.output.aggregated_masked_consensus
    output:
        masked_consensus_with_ref = temp(config["summary_dir"] + "/All-masked-consensus_with_ref.fa"),
        aligned_masked_consensus = config["summary_dir"] + "/All-masked-consensus_aligned.fa"
    params:
        ref = config["ref_genome"]
    log:
        "logs/alignMaskedConsensus.log"
    shell:
        r"""
        cp {input.aggregated_masked_consensus} {output.masked_consensus_with_ref} 2>{log}

        cat {params.ref} >> {output.masked_consensus_with_ref} 2>{log}

        mafft {output.masked_consensus_with_ref} > {output.aligned_masked_consensus} 2>{log}
        """


rule make_climb_dir:
    input:
        sample_bams = expand(config["samples_dir"] + "/{sample_dir}/{sample}_trimx2_sorted.bam", zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES),
        masked_consensus = expand(config["samples_dir"] + "/{sample_dir}/{sample}_masked_consensus.fa", zip, sample_dir=SAMPLE_DIRS, sample=SAMPLES)
    output:
        climb_bam = expand(config["batch_dir"] + "/ClimbSeq/{sample}/{sample}.bam", sample=SAMPLE_DIRS),
        climb_fasta = expand(config["batch_dir"] + "/ClimbSeq/{sample}/{sample}.fa", sample=SAMPLE_DIRS)
    params:
        sample_path = config["samples_dir"],
        climb_dir = config["batch_dir"] + "/ClimbSeq"
    log:
        "logs/makeClimbDir.log"
    shell:
        r"""
        bam_array=({input.sample_bams})

        for sample_bam in ${{bam_array[@]}}
        do

        newName=$(echo $sample_bam | rev | cut -d / -f 1 | rev | cut -d _ -f 1)
        cp $sample_bam {params.climb_dir}/$newName/$newName.bam
        awk '{{if($1~/>/){{split(FILENAME,a,"_"); print ">"a[1];}} else print}}' {params.sample_path}/$newName/$newName*_masked_consensus.fa > {params.climb_dir}/$newName/$newName.fa

        done 2>{log}
        """
