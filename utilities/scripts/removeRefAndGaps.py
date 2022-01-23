#Developed by Simon Daldry at University of Glasgow

from Bio import SeqIO
from Bio.Seq import Seq
import sys

def main(in_fasta, ref_name, out_fasta):
    '''
    Description
    -----------
    Script which removes reference 'MN908947.3' from alignment and removes gaps.
    Intended to work with alignment of one sequence and a reference.

    Parameters
    ----------
    in_fasta : fasta
        Aligned fasta which contains a sequence and a reference.
    ref_name : string
        Fasta header for reference sequence(s) are aligned to.
    out_fasta : fasta
        Fasta to write.
    '''
    with open(out_fasta, 'w') as fasta_handle: 					#create fasta file to write to
        for record in SeqIO.parse(in_fasta, 'fasta'): 				#loop through records in fasta
            if record.id != str(ref_name): 					#if record isn't reference
                record.seq = Seq(str(record.seq).upper().replace('-', '')) 	#remove gaps
                SeqIO.write(record, fasta_handle, 'fasta') 			#write record to output
    fasta_handle.close()


if __name__ == "__main__":
    in_fasta = sys.argv[1]
    ref = sys.argv[2]
    out_fasta = sys.argv[3]
    main(in_fasta, ref, out_fasta)
