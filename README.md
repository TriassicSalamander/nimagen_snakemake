# nimagen_snakemake
A nimagen workflow implemented in snakemake.

## Install
```
git clone git@github.com:TriassicSalamander/nimagen_snakemake.git
cd nimagen_snakemake
mamba env create -f environment.yml
conda activate nimagen_snakemake
```

## Creating Environment
It is recommended to use mamba (https://github.com/mamba-org/mamba) to create the environment, rather than conda.<br/>
Mamba can be installed by:
```
conda install mamba -n base -c conda-forge
```

## Using autoNimagen
autoNimagen.sh is a wrapper which automatically executes the snakemake pipeline once the sequencing run is complete.
This is useful when you know the output path of a given sequencing run, but don't know exactly when the run will finish.<br/>
Usage:
```
bash autoNimagen.sh -i <path/to/run/directory> -c <number of cores>
```


## Rulegraph
![Alt text](documentation/rulegraph.svg "Rulegraph")
