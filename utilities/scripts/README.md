# Scripts
For more information, check comments at head of scripts.

## ampliconDepth-NimaGen.sh


## getAlignmentStats


## getAmbPosCounts.py
Short for get ambiguous position counts. <br/> 
Takes csv of ambiguous nucleotide positions as input. <br/>
Writes a csv of the number of ambiguous nucleotides for each sample. <br/>

## getAmbiguousPositions.py
Takes aligned fasta file and name of the reference in alignment as input. <br/>
Writes a csv with columns: Sequence ID, Ambiguous Character and Position. <br/>
Writes a CSV with the number of Ns per sample. <br/>

## getAmpCoverageNG.R


## getCoverage.R


## maskAmbNucs.py
Script will replace ambiguous nucleotides with 'N' in regions of fasta sequences specified by regions file. <br/>

## process_freyja_demixed.py
Takes 'All-freyja-demixed.tsv' as input <br/>.
Writes a tsv with the lineages for each sample and what condition that lineage is speccific to. <br/>
Writes a csv with the abundance a given lineage has in a given sample and condition. Also states whether that lineage is also found in another condition for a given sample. <br/>

## removeRefAndGaps.py
Takes aligned fasta file and name of the reference in alignment as input. <br/>
Writes fasta file without reference sequence and gaps. <br/>

