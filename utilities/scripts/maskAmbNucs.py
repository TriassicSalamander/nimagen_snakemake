#Developed by Simon Daldry at University of Glasgow

from Bio import SeqIO
import sys


def main(regions, in_fasta, out_fasta):
    '''
    Description
    -----------
    Script will replace ambiguous nuceotides with 'N' in regions of fasta sequences specified by regions file.

    Parameters
    ----------
    regions : Tab-separated file.
        Should be similar in format to a bed file. Start and stop are expected to be in columns 2 and 3, respectively.
    in_fasta : Fasta file.
        Expects an alignment where sequences are aligned to reference specified in column 1 of regions file.
    out_fasta : Fasta file.
        Fasta file with ambiguous nucleotides replaced with N in specified regions.
    '''
    ambig_chars = ['B', 'D', 'H', 'K', 'M', 'R', 'S', 'V', 'W', 'Y']   #list of ambiguous nucleotides

    with open(regions, 'r') as reg_handle, open(out_fasta, 'w') as fasta_handle:   #open file containing coordinates for regions
        for region in reg_handle.readlines():   #loop through regions
            start = int(region.split('\t')[1])   #get region start
            end = int(region.split('\t')[2])   #get region end
            for record in SeqIO.parse(in_fasta, format='fasta'):   #loop through sequences in fasta file
                record.seq = record.seq.upper()   #makes sequence uppercase, if it wasn't already
                pre_reg = record.seq[:start]   #get sequence substring before region
                post_reg = record.seq[end:]
                for char in ambig_chars:   #loop through ambiguous nucleotides set
                    if char not in record.seq[start:end]:   #if char isn't in region, go to next char
                        continue
                    else:
                        new_reg = str(record.seq[start:end]).replace(char, 'N')   #replace ambig char with N within region
                        record.seq = pre_reg + new_reg + post_reg   #concatenate substrings
                SeqIO.write(record, fasta_handle, "fasta")   #write record to out_fasta


if __name__ == "__main__":
    reg_file = sys.argv[1]
    in_fasta = sys.argv[2]
    out_fasta = sys.argv[3]
    main(reg_file, in_fasta, out_fasta)
