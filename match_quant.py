import os
from re import S
import sys 
import csv
import glob
import pickle 
import numpy as np
import pandas as pd
from tqdm import tqdm
from plotnine import *

with open('/Users/vikas/Documents/UCSC/lab/kim_lab/TEBAG/TEbag_DB.pkl', 'rb') as f:
    TE_DB = pickle.load(f)

quant_file = pd.read_csv (sys.argv[2], sep='\t', header=0)

# removing non-TE rows
TE_quants = quant_file[quant_file["Name"].str.contains("ENST")==False]
TE_quants = TE_quants[TE_quants["Name"].str.contains("\)n")==False]
TE_quants = TE_quants.reset_index(drop = True)

# set TE TPM cutoff to 50
TE_quants = TE_quants[TE_quants["TPM"] > 50]
TE_quants = TE_quants.reset_index(drop = True)


species = ['Human','Chimp','Gorilla', 'Orangutan','Bonobo']
species_TE_df = pd.DataFrame(index=range(0),columns=species)


def check_if_in_species(element, lookup_table):
    TE = element[0]
    lookup_index = TE_DB.loc[TE_DB['Human'] == TE]
    return(lookup_index)

def append_to_new_df(index, lookup_table, append_df):
    species_lookup = check_if_in_species(index, lookup_table)
    #print(species_lookup)
    return species_lookup
    #append_df = pd.concat([append_df,species_lookup])
    #return append_df

    
#LI = check_if_in_species(TE_quants.loc[1], TE_DB)
#species_TE_df = pd.concat([species_TE_df,LI])

tqdm.pandas()
species_TE_list = [i for i in TE_quants.progress_apply(append_to_new_df, lookup_table = TE_DB, append_df = species_TE_df,  axis = 1)]
species_TE_df = pd.concat(species_TE_list)
species_TE_df = species_TE_df.reset_index(drop = True)

merged_species_TE_df = species_TE_df.merge(TE_quants, left_on = 'Human', right_on = 'Name')

base = sys.argv[1]
filename = base + ".merged_species_TE_df.TEBAG.csv"

merged_species_TE_df.to_csv(filename)

#from IPython import embed; embed()

#melted = pd.melt(merged_species_TE_df)

#ggplot(merged_species_TE_df, aes(x='Human', y='TPM')) + geom_point()
