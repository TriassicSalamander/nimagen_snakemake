# nimagen_snakemake
A nimagen workflow implemented in snakemake.

## Install
```
git clone git@github.com:TriassicSalamander/nimagen_snakemake.git
cd nimagen_snakemake
mamba env create -f environment.yml
mamba create --name pangolin --channel conda-forge pangolin
conda activate nimagen_snakemake
```

## Creating Environments
It is recommended to use mamba (https://github.com/mamba-org/mamba) to create the environment, rather than conda.<br/>
Mamba can be installed by:
```
conda install mamba -n base -c conda-forge
```

## Processing Waste Water Samples
Currently, whether or not freyja is used in the pipeline is controlled by commenting/uncommenting the 'Freyja Outputs' in 'rule all' in the Snakefile. <br/>
See comments in Snakefile for more detail.


## Using autoNimagen
autoNimagen.sh is a wrapper which automatically executes the snakemake pipeline once the sequencing run is complete.
This is useful when you know the output path of a given sequencing run, but don't know exactly when the run will finish.<br/>
Usage:
```
bash autoNimagen.sh -i <path/to/run/directory> -c <number of cores>
```


## Pipeline Overview
### Read Alignment
The pipeline begins by demultiplexing the samples using bcl-convert, generating fastqs. <br/>
These fastqs are then quality trimmed using TrimGalore. <br/>
Prinseq is used to trim bases with quality <30 from the 3' end. <br/>
The trimmed reads are aligned to a reference using bwa mem. <br/>
The primers are trimmed from the aligned reads using ivar trim. <br/>
The consensus is generated from the aligned reads using ivar consensus. <br/>



## Config Notes
### Input Paths
-**sample_sheet**: Path to sample sheet generated as part of the illumina sequencing run. Used by bcl-convert for demultiplexing and in the Snakefile for getting sample names. Normally found within the run directory. <br/>
-**run_directory**: Path to the output directory generated from the sequencing run. Contains multiplexed basecall files. <br/>
-**fastq_directory**: Path to directory which will contain demulitplexed fastq files. Demultiplexing is done as the first step in the pipeline. <br/>

### Output Paths
-**output_path**: Output path which should contain other output directories. <br/>
-**batch_dir**: The name of your run. Recommended to be contained within output_path, e.g. output/batch_1 (assuming output_path is 'output'). <br/>
-**samples_dir**: The name of the directory containing your samples. Recommended to be contained within batch_dir, e.g. output/batch_1/Samples <br/>
-**summary_dir**: The name of the directory containing your summary and statistics files. Recommended to be contained within batch_dir, e.g. output/batch_1/Summary <br/>

### Other Parameters
-**threads**: Number of threads to use in rules where multithreading is possible. Deafult: 7. <br/>
-**bcl-convert_threads**: Number of threads to be used by bcl-convert for demultiplexing. Default 26. <br/>

### Resources
-**ref_genome**: Path to reference genome used for alignment. Default: resources/ref/MN908947.3.fa <br/>
-**primer_bed**: Path to primer bed file. Default: resources/primers/nCoV-2019-NimaGen-V402.bed (for primer version 4.02). <br/>
-**PE_primer_bed**: Path to paired end primer bed file. Default: resources/primers/nCoV-2019-NimaGen-PE-V402.bed (for primer version 4.02). <br/>
-**ambig_regions**: Path to mask file specifying which regions to mask in consensus sequences. Format should be same as primer_bed. See resources/mask/primer_V3_ambig_regions for example. Default: resources/mask/blank_mask. Blank mask is used as default as ambiguous nucleotides was mainly an issue while V3 primers were being used. <br/>

### Utilities
-**scripts**: Path to scripts used by pipeline. Default: utilities/scripts <br/>




## Rulegraph
![Alt text](documentation/rulegraph.svg "Rulegraph")
