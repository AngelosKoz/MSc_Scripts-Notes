# %%
# %% Reading a file fast
import itertools, time

start_time = time.time()

file = open('e_coli.fa', 'r')
ecoli = ''
count = 0
k = 10
for line in file:

    count += 1
    if (count > 1): # the first line contains the non-sequence header so we discard it
        ecoli += line.replace("\n", "") # we string the newline
kmers = [ecoli[i:i+k] for i in range(len(ecoli)-k+1)]
kmers.sort()
setkmers = set(kmers)

tenmers = set([''.join(list(motf)) for motf in itertools.product(*itertools.repeat(['A','T','G','C'], 10))])
non_existing = list(tenmers - setkmers)

print(len(non_existing))

end_time = time.time()
print(end_time - start_time)
# %% Recursion euclid

def rec_euclid(a,b):
    if a % b == 0:
        return b
    else:
        print(b, a%b, end = 'hey\n')
        return rec_euclid(b, a%b)
x = 1920
y = 1080
rec_euclid(x, y)



# %% Simple Sorting

def SimpleSort(N):
    S=[]
    i=0
    minind=0
    
    while (i < len(N)):
        minimum = N[i]
        print(minimum)
        j = i

        while (j < len(N)):
            
            if minimum >= N[j]:
                
                minimum = N[j]
                print(N[j], end = '  minimum \n')
                
                minind = j
                print(minind, end = ' minind \n')
                
                j += 1
                print(j, end = ' if j \n')
            
            else:
                
                j += 1
                print(j, end = ' else \n')
        
        S.append(N[minind])
        print(S)
        N.remove(N[minind])
    return(S)
numbers = [2,7,14,3,12,9,11,6]
SimpleSort(numbers)


# %% Recursion merge-sort
def merge(a,b):
# Function to merge two arrays
    c = []
    while len(a) != 0 and len(b) != 0:
        if a[0] < b[0]:
            c.append(a[0])
            a.remove(a[0])
        else:
            c.append(b[0])
            b.remove(b[0])
    if len(a) == 0:
        c += b
    else:
        c += a
    return c
def mergesort(x):
# Function to sort an array using merge sort algorithm
    if len(x) == 0 or len(x) == 1:
        return x
    else:
    # print(len(x)) # use this line to get an idea how this works
        middle = int(len(x)/2)
        a = mergesort(x[:middle])
        b = mergesort(x[middle:])
        return merge(a,b)

numbers = [2,7,14,3,12,9,11,6]
sortedNs = mergesort(numbers)

print(sortedNs)




# %%
# Naive GC content
def naiveGC(genomefile):
    import re
    file = open(genomefile, 'r')

    seq = ""
    window = 1000
    total = 0
    nG = nC = 0
    GCCont = 0
    times = 0
    count = 0
    for line in file:
        count += 1
        if (count > 1): # the first line contains the non-sequence header so we discard it 
            length=len(line)
            total=total+length
            seq=seq+line[0:length-1]

    for k in range(len(seq)):
        if(seq[k]=="G"):
            nG+=1
        elif(seq[k]=="C"):
            nC+=1
    GCContent=(nG+nC)/len(seq)
    return(GCContent)
EcoliGC = naiveGC('files/ecoli.fa')
StaurGC = naiveGC('files/Staaur.fa')
print(EcoliGC, StaurGC)



# %%
def fastGC(genomefile):

    import regex as re

    file = open(genomefile, 'r')

    seq = ""
    total = 0
    nG = nC = 0
    GCCont = 0
    times = 0
    count = 0
    for line in file:
        count += 1
        if (count > 1):
            length = len(line)
            total = total+length
            seq = seq+line[0:length-1]
    file.close()
    
    nC = seq.count("C")
    nG = seq.count("G")
    GCContent = (nG+nC)/len(seq);
    return(GCContent)
EcoliGC = fastGC('files/ecoli.fa')
StaurGC = fastGC('files/Staaur.fa')
print(EcoliGC, StaurGC)




# %%

f=open('files/GCContent.tsv', 'r')

i=0
GCC={}
for line in f:
    i=i+1
    if(i>1):
        species=line.split()[0]
        GC=line.split()[1]
        GCC[species]=float(GC)

gcdistances={}
for genome1 in GCC.keys():#Gia kathe gonidiwma sto dict pairnei to epomeno kai ypologizei distance memtaksy twn 2
    for genome2 in GCC.keys():
        pair=genome1+":"+genome2
        gcdistances[pair]=abs(float(GCC[genome1])-float(GCC[genome2]))
        gcdistances[pair]=round(gcdistances[pair],2)

sorted(gcdistances.items(), key=lambda x: x[1])


# %%
## Clustering of a dataset 

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

# Load the dataframe and assign values/labels
df = pd.read_csv('files/GCContent_simple.csv')
dvalues = df['GCContent'].values.reshape(-1,1)
dlabels = list(df['Genome'])

# Calculate the distances
distances = pdist(dvalues)

# Convert the pairwise distances into a square distance matrix
distance_matrix = squareform(distances)

# Calculate the linkage matrix using Ward's method
linkage_matrix = linkage(distance_matrix, method='ward')

# Plot the dendrogram
sns.set_style('white')
dendrogram(linkage_matrix, labels=dlabels, color_threshold=0, orientation='left')

# Show the plot
plt.show()




# %%# Prepei na ypologisw ena parathyro sto opoio tha upologizw to GC kathws an exei ginei metafora tha doyme kapoioy eidoys peak ___|-|__
def windowGC(genomefile, window):

    import regex as re
    f=open(genomefile, 'r')

    seq = ""
    nG = nC = 0
    total = 0
    window = int(window)

    for line in f:
        x=re.match(">", line)
        if x == None:
            length=len(line)
            total=total+length
            seq=seq+line[0:length-1]
    f.close()

    step=int(window/10) # we use a 10% sliding overlap between windows
    times=int(len(seq)/step);

    GCwin = {}
    for i in range(times): 
        DNA=seq[i*step:i*step+window]
        nC=DNA.count("C")
        nG=DNA.count("G")
        GCwin[i*step] = (nG+nC)/window
    
    return(GCwin)
EcoliWGC = windowGC("files/ecoli.fa", 10000)
print(list(EcoliWGC.values())[1:10])


# %%

def Z_GC(genomefile, window, threshold):
    
    import regex as re
    import numpy as np
    f=open(genomefile, 'r')

    seq = ""
    nG = nC = 0
    total = 0
    window = int(window)

    for line in f:
        x=re.match(">", line)
        if x == None:
            length=len(line)
            total=total+length
            seq=seq+line[0:length-1]
    f.close()

    step=int(window/10) # we use a 10% sliding overlap between windows
    times=int(len(seq)/step);

    GCwin = {}
    for i in range(times): 
        DNA=seq[i*step:i*step+window]
        nC=DNA.count("C")
        nG=DNA.count("G")
        GCwin[i*step] = (nG+nC)/window
    
    GCcont =  list(GCwin.values())
    mGC=np.mean(GCcont)
    sdGC=np.std(GCcont)
    zGC=(GCcont-mGC)/sdGC #gc - mesh timh / tupikh
    for i in range(len(zGC)):
        if abs(zGC[i]) >= threshold:#elegxw an kathe timh ths zgc einai megaluterh apo thresh, kai an nai tupwse
            print(i*step, zGC[i])
Z_GC('files/Staaur.fa', 1000, 5) 
#anti na anaferontai ta parathura ena ena, na elegxw an einai sunexomena, an nai tote na katagrafw to euros. Px apo 608000 ews 614000 pws tha ebgaza mesh timh

# %% lathos tropos

def kmers(genomefile, k):
    import regex as re

    file = open(genomefile, 'r')

    seq = ""
    kmertable = {} 

    count = 0
    for line in file:
        count +=1
        if (count > 1) :
            length=len(line)
            seq=seq+line[0:length-1]
            
    file.close()

    seq = re.sub("[^AGCT]", "", seq)

    for i in range(len(seq)-k):
        DNA=seq[i:i+k]
        if DNA not in kmertable.keys():
            kmertable[DNA]=1
        else:
            kmertable[DNA]+=1

    kmertable = {k: float(v) / len(seq) for k, v in kmertable.items()}
    return(kmertable)


# %%


# %% Exercise 1 : Finding palindromes in a sequence using recursion

#In order to find palindromes we can aproach it by spliting a string in half and see if the 1st counterpart is the same as the reverse of the second counterpart.
#Of course we have to be careful since the aproach for even and odd subsequence length has to be different


# %% Recursion merge-sort
def merge(a,b):
# Function to merge two arrays
    c = []
    print('i am here\n')
    while len(a) != 0 and len(b) != 0:
        if a[0] < b[0]:
            c.append(a[0])
            print(c, end = ' if \n')
            a.remove(a[0])
        else:
            c.append(b[0])
            print(c, end = ' else\n')
            b.remove(b[0])
    if len(a) == 0:
        print(c, end = ' before b\n')
        c += b
        print(b)
        print(c, end = ' c from b \n')
    else:
        print(c, end = ' before a\n')
        c += a
        print(a)
        print(c, end = ' c from a \n')
    return c
def mergesort(x):
# Function to sort an array using merge sort algorithm
    if len(x) == 0 or len(x) == 1:
        print(f' {x} returnreturnreturnreturn\n')
        return x
    else:
        #print(len(x), end = ' len of x\n') # use this line to get an idea how this works
        middle = int(len(x)/2)
        print(middle, end = ' middle\n')
        print(x[:middle], end = 'aaaaaaaaaaaaaaaaaaaaaaaaaeskdjfksjdfksjdkfjskdjfksjdkf\n')
        a = mergesort(x[:middle])
        print(a, end = ' aaaa sfkdjfksjdkjf \n')
        print(x[:middle], end = 'bbbbbbbbbbbbbbbbbbbbbbbbbbbeskdjfksjdfksjdkfjskdjfksjdkf\n')
        b = mergesort(x[middle:])
        print(b, end = ' bbbb sfkdjfksjdkjf \n')
        return merge(a,b)

numbers = [14,7,3,12,9,11,6,2]
sortedNs = mergesort(numbers)

print(sortedNs)
# %%
