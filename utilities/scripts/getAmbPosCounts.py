#Developed by Simon Daldry at University of Glasgow

import pandas as pd
from collections import Counter
import sys


def main(amb_nuc_pos, amb_pos_count):
    '''
    Script which expects csv of ambiguous nucleotide positions as input (argument 1)
    and writes a csv (argument 2) of the number of positions for each sample.
    '''
    df = pd.read_csv(amb_nuc_pos, sep=',')                                     	#read in input csv
    pos_counts = Counter(df['Sequence ID'])                                     #make Counter object which counts the number of positions for each sample
    with open(amb_pos_count, 'w') as out_handle:                                #create out_csv to write to
        out_handle.write('Sequence ID,Amb Count\n')                             #write header
        for sample in df['Sequence ID'].unique():                               #loop through unique samples
            out_handle.write(sample + ',' + str(pos_counts[sample]) + '\n')     #write sample and position count to out_csv
        out_handle.write('Total,' + str(sum(pos_counts.values())))              #write total positions to out_csv
    out_handle.close()


if __name__ == "__main__":
    in_csv = sys.argv[1]
    out_csv = sys.argv[2]
    main(in_csv, out_csv)
