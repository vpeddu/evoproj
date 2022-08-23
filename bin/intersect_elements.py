import os
import sys 
import csv
import glob
import pickle 
import pandas as pd

class repeat_element():
    def __init__(self, chromosome, start, stop, name, score, strand):
        self.chromosome = chromosome
        self.start = start
        self.stop = stop
        self.name = name
        self.score = score
        self.strand = strand

def read_bed(bedpath): 
    print('reading bed ' + str(bedpath))
    element_dict = {}
    with open(bedpath) as openbed:
        lines = [line for line in openbed]
    for element in lines: 

        ename = element.split()[3]
        echr = element.split()[0]
        estart = element.split()[1]
        estop = element.split()[2]
        escore = element.split()[4]
        estrand = element.split()[5]
        
        try: 
            out_repeat_element.chromosome.append(echr)
            out_repeat_element.start.append(estart)
            out_repeat_element.stop.append(estop)
            out_repeat_element.name.append(ename)
            out_repeat_element.score.append(escore)
            out_repeat_element.strand.append(estrand)
        except:
            out_repeat_element = repeat_element([echr], [estart], [estop], [ename], [escore], [estrand])
    return out_repeat_element

            
    #tmpname = repeat_element(echr, estart, estop, ename, escore, estrand)

def find_shared(new_species, old_species):
    new = set(new_species.name)
    old = set(old_species.name)
    new_old_intersect = new.intersection(old)
    return(new_old_intersect)
    
def find_new_unique(new_species, old_species):
    new = set(new_species.name)
    old = set(old_species.name)
    new_unique = new - old
    return(new_unique)
    
def find_old_unique(new_species, old_species):
    new = set(new_species.name)
    old = set(old_species.name)
    old_unique = old - new
    return old_unique
    
def create_intersections(ns, os):
    s = find_shared(ns,os)
    n = find_new_unique(ns,os)
    o = find_old_unique(ns,os)
    return [s, n, o]
    
def check_if_intersects(df_index, comp_object): 
    if df_index['element'] in comp_object[0]:
        return True
    else:
        return False

human_bed_path = sys.argv[1]
lifted_files = glob.glob("*.check.lifted.bed")
species_comparisons = {}
human_bed = read_bed(human_bed_path)

for file in lifted_files: 
    species_name = file.split('.')[0]
    current_bed = read_bed(file)
    species_comp = create_intersections(human_bed, current_bed)
    species_comparisons[species_name] = species_comp
    
species = species_comparisons.keys()
intersection_table = pd.DataFrame(index=range(len(human_bed.name)),columns=species)
intersection_table['element'] = human_bed.name

for comp in species_comparisons.keys():
    intersection_table[comp] = intersection_table.apply(check_if_intersects, comp_object = species_comparisons[comp], axis = 1)
    
appended_column_names = list(intersection_table.columns)
for i in range((len(appended_column_names) - 1)):
    appended_column_names[i] = ('species_' + str(appended_column_names[i]))

intersection_table.columns = appended_column_names

intersection_table.to_csv("intersection_table.csv")

TE_DB = open('TEbag_DB.pkl', 'wb') 
pickle.dump(intersection_table, TE_DB)

#from IPython import embed; embed()
