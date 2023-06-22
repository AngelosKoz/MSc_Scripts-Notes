import itertools, re, time
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

start_time = time.time()

def get_lines(fl):
    with open(f'{fl}') as get_lines:
        line_number = list(enumerate(get_lines))[-1][0]+1
    return line_number

whole_seq = ''
with open('e_coli.fa') as e_col:
    for i in range(get_lines('e_coli.fa')):
        one_line = e_col.readline()
        whole_seq = ''.join([whole_seq,one_line])
whole_seq = re.sub(r'\n','',whole_seq)
eightmers = itertools.product(*itertools.repeat(['A','T','G','C'], 8))
dict_mers = {}
for rnd in eightmers:
    dict_mers[''.join(list(rnd))] = 0 #create a dictionary with all 65536 possible 8mers

for j in range(len(whole_seq)-8+1): #we assign each time we find an 8mer to the dictionary, moving at a pace of 1 nucleotide from left to right at a time
    motif = j+8
    e_mer = whole_seq[j:motif]
    dict_mers[e_mer] += 1

end_time = time.time()
print(end_time - start_time)
#filtered_dict_mers = {k:v for k,v in dict_mers.items() if v != 0} #65360 ... Only 176 not found! # We would use this to filter, but no need to here
all_kmer_occurences = list(dict_mers.values())
occurence_counter = dict(Counter(all_kmer_occurences))






# %% ################--------- Project 10 ---------##################

def generate_seq():
    s = ''
    with gzip.open(f'./project_10_aggk/chr22.fa.gz', 'rt') as f:
        for i, l in enumerate(f):
            if i==0:
                continue
            if i%50_000 == 0:
                print (f'{i}')
            s += l.strip().replace('N','')
    return s


def get_motif():
    all_seq = generate_seq()
    for j in range(3,100): #I chose 3 because there is no way 2-mers wont exist (not that 3 changes it but still)
        kmers = itertools.product(*itertools.repeat(['A','T','G','C'], j))
        dict_mers = {}
        for rnd in kmers:
            dict_mers[''.join(list(rnd))] = 0 
        for n in dict_mers.keys():
            if n in all_seq:
                continue
            else:
                got_it = f'{n} is the smallest motif not found'
                return got_it
                
print(get_motif())