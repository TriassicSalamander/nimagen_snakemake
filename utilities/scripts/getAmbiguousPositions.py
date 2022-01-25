#!/usr/bin/python
#Script which expects a single file containing one or more fastas as input
#Writes a CSV with headers: Sequence ID, Ambiguous Character and Position
#Developed by Simon Daldry at University of Glasgow

from Bio import SeqIO
import sys


def get_index_positions(list_of_elems, element):
    '''
    Returns the indexes of all occurrences of give element in
    the list- listOfElements
    '''
    index_pos_list = []
    index_pos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            index_pos = str(list_of_elems).upper().index(element, index_pos)
            # Add the index position in list
            index_pos_list.append(index_pos + 1)
            index_pos += 1
        except ValueError:
            break
    return index_pos_list


def get_ref_seq_aligned_map(fasta, ref_name):
    '''
    Function which returns a mapping of where a base in an aligned reference
    corresponds to in the unaligned reference.
    '''
    ref_index = SeqIO.index(fasta, format='fasta') 						#parse aligned fasta file index
    ref_gappy_seq = ref_index[str(ref_name)].seq 						#get aligned reference sequence
    base_count = 0 										#init counter which will track bases in aligned sequence
    position_map = {} 										#init dict which will store position in aligned seq as keys and corresponding position in unaligned(gapless) seq as values

    for pos,char in enumerate(ref_gappy_seq): 							#loop through aligned seq
        if char != '-': 									#if character isn't a gap
            base_count += 1 									#increment counter
        position_map[(pos+1)] = base_count 							#store key:value pair to position map. pos+1 used because it's 0-indexed.

    return position_map


def main(fasta, ref_name, amb_tsv, N_counts):
    ambig_chars = ['B', 'D', 'H', 'K', 'M', 'R', 'S', 'V', 'W', 'Y'] 				#list of ambiguous nucleotide characters
    pos_map = get_ref_seq_aligned_map(fasta, ref_name)						#get map of where a base in an aligned seq corresponds to in the unaligned reference
    with open(amb_tsv,'w') as out_ambs, open(N_counts, 'w') as out_Ns: 				#create files to write to
        out_ambs.write('Sequence ID,Ambiguous Base,Position\n') 				#create header line
        out_Ns.write('Sequence ID,N Count\n') 							#create header line for N count file
        N_total = 0 										#init counter for total Ns across all samples
        for record in SeqIO.parse(fasta, format='fasta'): 					#loop through sequences in fasta file
            for char in ambig_chars: 								#loop through ambiguous nucleotides set
                positions = get_index_positions(record.seq, char) 				#get list of indeces for a given ambiguous nucelotide character in a given sequence
                for pos in positions: 								#loop through indeces
                    out_ambs.write(record.id + ',' + char + ',' + str(pos_map[pos]) + '\n') 	#write sequence id, ambiguous nucleotide character and index to output file. index is position in unaligned reference.
            Ns = str(record.seq).upper().count('N') 						#get number of Ns in current record
            N_total += Ns
            out_Ns.write(record.id + ',' + str(Ns) + '\n') 					#write number of Ns in record to file
        out_Ns.write('Total,' + str(N_total)) 							#write total Ns to file
    out_ambs.close()
    out_Ns.close()


if __name__ == "__main__":
    in_fasta = sys.argv[1]
    ref_genome = sys.argv[2]
    out_amb_tsv = sys.argv[3]
    out_N_counts = sys.argv[4]
    main(in_fasta, ref_genome, out_amb_tsv, out_N_counts)
