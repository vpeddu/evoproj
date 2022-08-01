import os
import sys 
import csv
import glob
import pickle 
import pandas as pd
from tqdm import tqdm


class repeat_element():
    def __init__(self, chromosome, start, stop, name, score, strand):
        self.chromosome = chromosome
        self.start = start
        self.stop = stop
        self.name = name
        self.score = score
        self.strand = strand

def read_bed(bedpath): 
    print('reading bed')
    element_dict = {}
    with open(bedpath) as openbed:
        lines = [line for line in openbed]
    for element in tqdm(lines, total = len(lines)): 

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
    


Human_bed = read_bed('/private/groups/kimlab/vikas/liftover_tool/lnc/g39_lnc.bed')
Human_Chimp_lifted_bed = read_bed('/private/groups/kimlab/vikas/liftover_tool/lnc/Human.Chimp.check.lnc.lifted.bed')
Human_Gorilla_lifted_bed = read_bed('/private/groups/kimlab/vikas/liftover_tool/lnc/Human.Gorilla.check.lnc.lifted.bed')
Human_Orangutan_lifted_bed = read_bed('/private/groups/kimlab/vikas/liftover_tool/lnc/Human.Orangutan.check.lnc.lifted.bed')
Human_Bonobo_lifted_bed = read_bed('/private/groups/kimlab/vikas/liftover_tool/lnc/Human.Bonobo.check.lnc.lifted.bed')

HC_comp = create_intersections(Human_bed, Human_Chimp_lifted_bed)
HG_comp = create_intersections(Human_bed, Human_Gorilla_lifted_bed)
HO_comp = create_intersections(Human_bed, Human_Orangutan_lifted_bed)
HB_comp = create_intersections(Human_bed, Human_Bonobo_lifted_bed)

species = ['Human','Chimp','Gorilla', 'Orangutan','Bonobo']
intersection_table = pd.DataFrame(index=range(len(Human_bed.name)),columns=species)
intersection_table['element'] = Human_bed.name

def check_if_intersects(df_index, comp_object): 
    if df_index[0] in comp_object[0]:
        return True
    else:
        return False

tqdm.pandas()
intersection_table['Chimp'] = intersection_table.progress_apply(check_if_intersects, comp_object = HC_comp, axis = 1)
intersection_table['Gorilla'] = intersection_table.progress_apply(check_if_intersects, comp_object = HG_comp, axis = 1)
intersection_table['Orangutan'] = intersection_table.progress_apply(check_if_intersects, comp_object = HO_comp, axis = 1)
intersection_table['Bonobo'] = intersection_table.progress_apply(check_if_intersects, comp_object = HB_comp, axis = 1)
intersection_table['Human'] = True

intersection_table.to_csv("intersection_table.csv")
# for ind in tqdm.tqdm(intersection_table.index, total = len(intersection_table.index)):
#     intersection_table['Chimp'] = check_if_intersects(intersection_table, ind, HC_comp)
#     intersection_table['Gorilla'] = check_if_intersects(intersection_table, ind, HG_comp)
#     intersection_table['Orangutan'] = check_if_intersects(intersection_table, ind, HO_comp)
#     intersection_table['Bonobo'] = check_if_intersects(intersection_table, ind, HB_comp)


# need to make csvwriter
#cw = csv.writer(open("chimp_gorilla_intersection.csv",'w'))
#cw.writerows(list(chimp_gorilla_intersections[0]))


TE_DB = open('lnc_TEbag_DB.pkl', 'wb') 
pickle.dump(intersection_table, TE_DB)

#from IPython import embed; embed()

