# Reading ecoli genome
#Edw den pairnw unique alla pairnw olo to set

file = open('files/ecoli.fa', 'r')
ecoli = ''
count = 0
k = 10
for line in file: #auth h loypa petaei th prwth grammh, thn antikathista kai to vazei olo se mia grammh
    count += 1
    if (count > 1): # the first line contains the non-sequence header so we discard it 
        ecoli += line.replace("\n", "") # we string the newline character from the end of each line

# Creating a sort list of all k-mers in the genome
kmers = [ecoli[i:i+k] for i in range(len(ecoli)-k+1)]
kmers.sort() #Edw exw ola ta k-mers (oxi unique) me alfabhtikh seira



## Dichotomous Search for k-mers

# Using time to measure time of execution
import time
start_time = time.time()

# Pattern search

pattern = 'AGTTAGGCCT' #pattern poy den yparxei epithdes
#pattern = 'TCGGCATCAG'
matches = 0

iter = 0
min = 1
max = len(kmers)

midpoint = int((max+min)/2)

import math 

while iter <= math.log2(len(kmers)):
    iter += 1
    if (pattern == kmers[midpoint]):
        matches = kmers.count(kmers[midpoint])
        print("Pattern matched ", matches, " times")
        break
    if (pattern > kmers[midpoint]):
        min = midpoint
        midpoint = int((max+min)/2)
    if (pattern < kmers[midpoint]):
        max = midpoint
        midpoint = int((max+min)/2)
if (matches == 0):
    print("No matches found")

print("--- %s seconds ---" % (time.time() - start_time))
