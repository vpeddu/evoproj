import os
import sys 
import csv
import glob
import pickle 
import pandas as pd
from tqdm import tqdm

with open('/Users/vikas/Documents/UCSC/lab/kim_lab/TEBAG/TEbag_DB.pkl', 'rb') as f:
    TE_DB = pickle.load(f)

quant_file = pd.read_csv ('/Volumes/metagenomics_drive/liftover_project/test_files/te_aware_quant_for_vikas.tsv', sep='\t', header=0)

# removing non-TE rows
TE_quants = quant_file[quant_file["Name"].str.contains("ENST")==False]
# set TE TPM cutoff to 50
TE_quants = TE_quants[quant_file["TPM"] > 50]
TE_quants = TE_quants.reset_index(drop = True)


species = ['Human','Chimp','Gorilla', 'Orangutan','Bonobo']
species_TE_df = pd.DataFrame(index=range(len(TE_quants)),columns=species)


def check_if_in_species (element, lookup_table):
    TE = element[0]
    lookup_index = TE_DB.loc[TE_DB['Human'] == TE]
    return(lookup_index)

from IPython import embed; embed()

#species_TE_df = pd.concat([species_TE_df,lookup_index])