#Developed by Simon Daldry at the University of Glasgow

import pandas as pd
from functools import reduce
from collections import Counter
import sys


'''
Expects 'All-freyja-demixed.tsv' which is produced by the nimagen_snakemake pipeline when the freyja outputs are unncommented in the Snakefile.

Outputs lineage counts and lineage abundances, where:
-lineage_counts is a tsv with the columns 'sample', 'Intersect' and a column for each condition.
	Conditions are taken from the first column in input 'All-freyja-demixed.tsv' and are expected to be of the format 'SAMPLE-CONDITION_freyja_variants.tsv'.
	For a given sample in lineage_counts, each column will be a list of the lineages specific to that column name.
-lineage_abundances is a csv with the columns 'Sample', 'In Intersect', 'Condition' and 'Abundance'.
	For a given sample in lineage_abundances, it will have the sample name, one of the lineages found in that sample
	for this condition, a Yes/No value of whether that lineage is in found in other conditions for that sample,
	what condition this lineage and sample belong to and the abundance of this lineage in this sample for this condition.

To use:
	python process_freyja_demixed.py 'All-freyja-demixed.tsv' 'BATCHNAME_lineage_counts.tsv' 'BATCHNAME_lineage_abundances.csv'

Note: if the genomics team changes the format of the sample names to not use a hyphen, you may have to change how the string method 'split' is used to get the sample name and condition.
'''


def main(freyja_demixed, out_lin_counts, out_lin_abundances):
    #Parse input data
    demixed_df = pd.read_csv(freyja_demixed, sep="\t")
    demixed_df.rename(columns={'Unnamed: 0':'sample'}, inplace=True)
    demixed_df.dropna(inplace=True)

    #create mapping used for formatting sample names to conditions
    conditions_map = {sample:sample.split('_')[0].split('-')[1] for sample in demixed_df["sample"] if 'Undetermined' not in sample}

    #create columns for output dataframe
    conditions = list(set([sample.split('_')[0].split('-')[1] for sample in demixed_df["sample"] if 'Undetermined' not in sample])) #make list of unique conditions
    #conditions.remove('S0') #remove value from 'Undetermined' sample
    conditions.sort() #sort to make output column order consistent on subsequent script runs
    df_cols = [cond + '-specific' for cond in conditions] #append '-specific' to end of condition names
    df_cols.insert(0,'Intersect') #add other column names
    # =============================================================================
    # df_cols.insert(0,'Lineage Count %')
    # df_cols.insert(0,'Lineage Count')
    # =============================================================================

    #create output df to append to
    out_df = pd.DataFrame(columns=df_cols)

    #create sample list
    samples = list(set([sample.split('-')[0] for sample in demixed_df["sample"] if 'Undetermined' not in sample])) #Make unique list of samples
    #samples.remove('Undetermined')
    samples.sort() #Sort samples

    #loop through samples
    for sample in samples:
        sample_df = demixed_df[demixed_df["sample"].str.contains(sample)] #get dataframe subset with given sample

        lins_df = sample_df[["sample","lineages","abundances"]] #from sample df, select columns of interest
        lins_dict = lins_df.set_index("sample").to_dict() #create dict with sample,lineages as key,value pair. still has 'lineage' column name as outer key.
        lins_dict = lins_dict["lineages"] #remove outer key
        lins_dict = {conditions_map[key]:val.replace('\n','').strip("'[]").split("' '") for key,val in lins_dict.items()} #creates dict of condition:lineages for a given sample, with more appropriate formatting

        out_dict = {'Intersect':reduce(set.intersection, (set(val) for val in lins_dict.values()))} #Initialise dict with intersect of lineages between conditions

        for key in lins_dict.keys(): #loop through lins_dict to get relative complements
            other_lins = [lins_dict[key] for key in (lins_dict.keys() - {key})] #get lineages for conditions except for current key
            other_lins = [val for sublist in other_lins for val in sublist] #collapse nested lists
            out_dict[key+'-specific'] = set(lins_dict[key]) - set(other_lins) #add relative complement to dict
            out_dict = {key:list(val) for key,val in out_dict.items()} #typecast to prevent 'set()' appearing in final output

    # =============================================================================
    #     all_lins = [lins_dict[key] for key in lins_dict.keys()] #get all lineages for conditions
    #     all_lins = [val for sublist in all_lins for val in sublist] #collaps nested list
    #     out_dict['Lineage Count'] = dict(Counter(all_lins)) #add count of lineages to dict
    #
    #     out_dict['Lineage Count %'] = {lin:(count/(len(lins_dict.keys()))) for lin,count in out_dict['Lineage Count'].items()} #add lineage counts as % of total conditions for sample
    # =============================================================================

        out_df = out_df.append(out_dict, ignore_index=True)

    out_df['sample'] = samples #add samples to dataframe
    out_df.set_index(['sample'], inplace=True) #set samples to index

    out_df.to_csv(out_lin_counts, sep='\t') #write to tsv




    #THIS BLOCK CREATES A MAPPING OF SAMPLE,CONDITION AND LINEAGE TO ABUNDANCE
    #in retrospect, could have probably used defaultdicts to handle nonexistant keys
    lin_freq_map = {} #init dict

    for index,row in demixed_df.iterrows(): #loop through rows
        if 'Undetermined' in row["sample"]:
            continue
        else:
            sample_name = row["sample"].split('-')[0]
            cond = row["sample"].split('_')[0].split('-')[1]

            if sample_name not in lin_freq_map.keys(): #create subdict for sample name, if it doesn't exist
                lin_freq_map[sample_name] = {}

            if cond not in lin_freq_map[sample_name].keys(): #create sub-subdict for condition
                lin_freq_map[sample_name][cond] = {}

            for lin_ind,lin in enumerate(row["lineages"].replace('\n','').strip("'[]").replace("'","").split()): #loop through lineages in row
                if lin not in lin_freq_map[sample_name][cond].keys(): #create sub-sub-subdict for lineage
                    lin_freq_map[sample_name][cond][lin] = {}

                lin_freq_map[sample_name][cond][lin] = row["abundances"].replace('\n','').strip("'[]").split()[lin_ind] #add lineage abundance as value for given sample, condition and lineage




    #THIS BLOCK CREATES A DATAFRAME WITH THE BELOW COLUMNS
    abund_cols = ['Sample', 'Lineage', 'In Intersect', 'Condition', 'Abundance']

    abunds_df = pd.DataFrame(columns=abund_cols)

    for sample in lin_freq_map.keys():
        if sample=="Undetermined":
            continue
        for cond in lin_freq_map[sample].keys():
            for lin in lin_freq_map[sample][cond].keys():
                if lin in out_df.loc[sample, 'Intersect']:
                    in_intersect = 'Yes'
                else:
                    in_intersect = 'No'

                abund = lin_freq_map[sample][cond][lin]

                new_row = {'Sample':sample, 'Lineage':lin, 'In Intersect':in_intersect, 'Condition':cond, 'Abundance':abund}

                abunds_df = abunds_df.append(new_row, ignore_index=True)

    abunds_df.sort_values(by=['Sample'], inplace=True)

    abunds_df.set_index(['Sample'], inplace=True)
    abunds_df.to_csv(out_lin_abundances)




if __name__ == "__main__":
    in_freyja = sys.argv[1]
    out_counts = sys.argv[2]
    out_abundances = sys.argv[3]
    main(in_freyja, out_counts, out_abundances)
