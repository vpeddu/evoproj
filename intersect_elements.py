from imp import cache_from_source
from multiprocessing import sharedctypes
from operator import ne
import os
import sys 
import glob
from cv2 import absdiff
from tqdm import tqdm

shared_elements = glob.glob('/Volumes/metagenomics_drive/liftover_project/lifted_beds/*.lifted.bed')
unique_elements = glob.glob('/Volumes/metagenomics_drive/liftover_project/lifted_beds/first.unlifted.bed')

#species = [set(i.split('.')[0] for i in lifted_beds)]
species = ['Human','Chimp','Gorilla']

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
    print(element_dict)
    

testhead = read_bed('test.bed')
testheadtail = read_bed('test2.bed')

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
    

print(shared_elements)

Human_Chimp_lifted_bed = read_bed('/Volumes/metagenomics_drive/liftover_project/lifted_beds/Human.Chimp.check.lifted.bed')
Human_Gorilla_lifted_bed = read_bed('/Volumes/metagenomics_drive/liftover_project/lifted_beds/Human.Gorilla.check.lifted.bed')

chimp_gorilla_intersections = create_intersections(Human_Chimp_lifted_bed, Human_Gorilla_lifted_bed)

from IPython import embed; embed()
