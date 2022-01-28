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
It is recommended to use mamba (https://github.com/mamba-org/mamba) to create the environment, rather than conda.
Mamba can be installed by:
```
conda install mamba -n base -c conda-forge
```

## Rulegraph
![Alt text](documentation/rulegraph.svg "Rulegraph")
