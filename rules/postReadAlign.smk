
#Rule which collects the consensus sequences from all of the sample folders into a single fasta file.
rule aggregate_consensus:
    input:
        consensus_files = expand(config["samples_dir"] + "/{sample_dir}/{sample_dir}_trimx2_ivar_consensus.fa", sample_dir=SAMPLE_DIRS)
    output:
        aggregated_consensus = config["summary_dir"] + "/All-consensus.fa"
    log:
        "logs/aggregateConsensus.log"
    shell:
        r"""
        sample_array=({input.consensus_files})   #Create array of all the sample consensus fastas

        IFS=$'\n'
        sorted_samples=($(sort <<<"${{sample_array[*]}}"))   #These three lines are for sorting the array.
        unset IFS

        for sample_consensus in ${{sorted_samples[@]}}   #Loop through the sorted array
        do
        fasta_header=$(echo $sample_consensus | rev | cut -d / -f 1 | cut -d _ -f 4- | rev)   #Get sample name from the file path. This will be used as the fasta header.
        awk -v fasta_header="$fasta_header" '{{if($1~/>/){{print ">"fasta_header;}} else print}}' $sample_consensus   #If a line begins with '>', print the new header. Otherwise print the sequence.
        done >{output.aggregated_consensus} \
        2>{log}
        """


#Rule which adds the reference to the collected fastas and then aligns them using mafft.
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


#Rule which adds the reference to the consensus fasta in each sample folder and then aligns each one using mafft.
rule align_individual_consensus:
    input:
        sample_consensus = config["samples_dir"] + "/{sample_dir}/{sample_dir}_trimx2_ivar_consensus.fa"
    output:
        consensus_with_ref = temp(config["samples_dir"] + "/{sample_dir}/{sample_dir}_consensus_with_ref.fa"),
        aligned_consensus = config["samples_dir"] + "/{sample_dir}/{sample_dir}_consensus_aligned.fa"
    params:
        ref = config["ref_genome"]
    log:
        "logs/alignIndConsensus/{sample_dir}.log"
    shell:
        r"""
        cp {input.sample_consensus} {output.consensus_with_ref} 2>{log}

        cat {params.ref} >> {output.consensus_with_ref} 2>{log}

        mafft {output.consensus_with_ref} > {output.aligned_consensus} 2>{log}
        """


#Rule which calls a script to mask regions with 'N's.
#This may be done if a given region, such as an amplicon, is known to contain a high number of ambiguous nucleotides.
#Regions are defined in a tsv file where start and end are columns 2 and 3, respectively.
rule mask_ambiguous_nucleotides:
    input:
        aligned_ind_consensus = rules.align_individual_consensus.output.aligned_consensus
    output:
        masked_consensus = config["samples_dir"] + "/{sample_dir}/{sample_dir}_masked_consensus_aligned.fa"
    params:
        script = config["scripts"] + "/maskAmbNucs.py",
        mask_regions = config["ambig_regions"]
    log:
        "logs/maskAmbigNucs/{sample_dir}.log"
    shell:
        r"""
        python {params.script} \
        {params.mask_regions} \
        {input.aligned_ind_consensus} \
        {output.masked_consensus} \
        2>{log}
        """


#Rule which calls a script to unalign the individual consensus sequences after the ambiguous regions have been masked.
#Unalign in this case means removing the reference sequence and gaps from the fasta.
rule unalign_masked_consensus:
    input:
        masked_consensus = rules.mask_ambiguous_nucleotides.output.masked_consensus
    output:
        unaligned_consensus = config["samples_dir"] + "/{sample_dir}/{sample_dir}_masked_consensus.fa"
    params:
        script = config["scripts"] + "/removeRefAndGaps.py",
        ref_name = config["ref_genome"].split('/')[-1][:-3]
    log:
       "logs/removeRefandGaps/{sample_dir}.log"
    shell:
        r"""
        python {params.script} \
        {input.masked_consensus} \
        {params.ref_name} \
        {output.unaligned_consensus} \
        2>{log}
        """


#Rule which collects the masked consensus sequences from all of the sample folders into a single fasta file.
#Identical to the aggregate_consensus rule, except for the input.
rule aggregate_masked_consensus:
    input:
        masked_consensus_files = expand(config["samples_dir"] + "/{sample_dir}/{sample_dir}_masked_consensus.fa", sample_dir=SAMPLE_DIRS)
    output:
        aggregated_masked_consensus = config["summary_dir"] + "/All-masked-consensus.fa"
    log:
        "logs/aggregateMaskedConsensus.log"
    shell:
        r"""
        sample_array=({input.masked_consensus_files})

        IFS=$'\n'
        sorted_samples=($(sort <<<"${{sample_array[*]}}"))
        unset IFS

        for sample_consensus in ${{sorted_samples[@]}}
        do
        fasta_header=$(echo $sample_consensus | rev | cut -d / -f 1 | cut -d _ -f 3- | rev)
        awk -v fasta_header="$fasta_header" '{{if($1~/>/){{print ">"fasta_header;}} else print}}' $sample_consensus
        done >{output.aggregated_masked_consensus} \
        2>{log}
        """


#Rule for making the directory that will be uploaded to the climb server.
rule make_climb_dir:
    input:
        sample_bams = expand(config["samples_dir"] + "/{sample_dir}/{sample_dir}_trimx2_sorted.bam", sample_dir=SAMPLE_DIRS),
        masked_consensus = expand(config["samples_dir"] + "/{sample_dir}/{sample_dir}_masked_consensus.fa", sample_dir=SAMPLE_DIRS)
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
        bam_array=({input.sample_bams})   #Create array of all the sample bam files.

        for sample_bam in ${{bam_array[@]}}   #Loop through array of bam files.
        do

        newName=$(echo $sample_bam | rev | cut -d / -f 1 | cut -d _ -f 3- | rev)   #Get sample name from the file path. This will be used as the new filename and fasta header.
        cp $sample_bam {params.climb_dir}/$newName/$newName.bam   #Copy the sample bam file to the climb directory.
        awk -v newName="$newName" '{{if($1~/>/){{print ">"newName;}} else print}}' {params.sample_path}/$newName/$newName*_masked_consensus.fa > {params.climb_dir}/$newName/$newName.fa

        done 2>{log}
        """
