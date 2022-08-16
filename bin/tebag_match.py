import os
import sys 
import csv
import glob
import pickle 
import numpy as np
import pandas as pd

with open(sys.argv[1], 'rb') as f:
    TE_DB = pickle.load(f)

#quant_file = pd.read_csv (sys.argv[2], sep='\t', header=0)
quant_file = pd.read_csv(sys.argv[2], sep = '\t', header = 0)
#from IPython import embed; embed()
# removing non-TE rows
TE_quants = quant_file[quant_file["insertion"].str.contains("ENST")==False]
TE_quants = TE_quants[TE_quants["insertion"].str.contains("\)n")==False]
TE_quants = TE_quants.reset_index(drop = True)

TE_quants = TE_quants.set_index(TE_quants['insertion'])
TE_quants = TE_quants.drop(['insertion'], axis = 1)

# remove rows where all values are 0 
TE_quants = TE_quants[(TE_quants > 0).any(1)]
TE_quants.reset_index(inplace=True)
TE_quants_dict = TE_quants.to_dict()
# set TE TPM cutoff to 50
#TE_quants = TE_quants[TE_quants["TPM"] > 0
#TE_quants = TE_quants.reset_index(drop = True)

species = list(TE_DB.columns)[:-1]
species.insert(0,'Human')
species_TE_df = pd.DataFrame(index=range(0),columns=species)

# def append_to_new_df(index, lookup_table, append_df):
#     species_lookup = check_if_in_species(index, lookup_table)
#     return species_lookup
#     #print(species_lookup)
#     #append_df = pd.concat([append_df,species_lookup])
#     #return append_df

TE_insertions = TE_quants['insertion']
TE_insertions = TE_insertions.to_frame()
TE_insertions = TE_insertions.merge(TE_DB, left_on = 'insertion', right_on = 'element')

# def check_if_in_species(row):
#     new_row = [row['Chimp'],row['Gorilla'],row['Orangutan'],row['Bonobo']]
#     return new_row

#TE_insertions['insertion_status'] = TE_insertions.apply(lambda row: check_if_in_species(row), axis=1)
TE_quants = TE_quants.merge(TE_insertions, left_on='insertion', right_on = 'element')
TE_quants = TE_quants.drop(columns = ['insertion_x','insertion_y'])
# LI = check_if_in_species(TE_quants.loc['AluSc_range=chrY:57202568-57202876_strand=+'], TE_DB)
# #species_TE_df = pd.concat([species_TE_df,LI])

# tqdm.pandas()
# species_TE_list = [i for i in TE_quants.progress_apply(append_to_new_df, lookup_table = TE_DB, append_df = species_TE_df,  axis = 1)]
# species_TE_df = pd.concat(species_TE_list)
# species_TE_df = species_TE_df.reset_index(drop = True)

# merged_species_TE_df = species_TE_df.merge(TE_quants, left_on = 'Human', right_on = 'Name')

filename =  "merged_species_TE_df.TEBAG.csv"

TE_quants.to_csv(filename, index = False)

#from IPython import embed; embed()

#melted = pd.melt(merged_species_TE_df)

#ggplot(merged_species_TE_df, aes(x='Human', y='TPM')) + geom_point()
