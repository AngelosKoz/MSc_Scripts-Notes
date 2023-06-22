# %%
#Apostash hamming: How much 2 set of same length strings differ. If 3 different nucleotides, d=3.
#Problem is it does not take into account the meaning that different positions have in the motif.
#For example GGG[AG][AC]TT[TC]CC, not all positions matter the same since [] means it can have more than 1 option

#PWM (Position wight Matrix) is a good way to evaluate each position in the motif, and that depends on 
# the probability of a nucleotide to appear in specific position. This is done for example, if we 
# have n=104 sequences of length 10, we check what nucleotide appears in each position, then divide 
# by n. So if first position has 103 times appearance of G then it will be 101/104 = 0.97. 
# The sum of each column should be 1

#- Steps -#
#1) Δηλώνω πίνακα n σημείων πρόσδεσης μήκους l, TFBS[n,l]
#2) Δηλωση πίνακα Ν τεσσάρων νουκλεοτιδίων (A,G,C,T)
#3) Δήλωση πίνακα Ρ πιθανοτήτων νουκλεοτιδίων (A,G,C,T)
#4) Δήλωση πίνακα PWM[4,l]

#Απαρίθμηση 1: 
#για θέση i=1 έως i=l ανά 1;
#   Δημιούργησε τη σειρά C=TFBS[1:n,i];
#   Απαρίθμηση 2: 
#   για θέση j=1 έως j=n ανά 1;
#       διάβασε s=C[j];
#       αύξησε το πλήθος πίνακα νουκλεοτιδίων N[s]++
#   Τέλος Απαρίθμηση 2

#   Συχνότητα:
#   Για κάθε νουκλεοτίδιο s;
#       Υπολόγισε τη συνότητα P[s]=N[s]/n
#       Απόδοση στον PWM[i,s]=P[s]
#   N=0; P=0; $αρχικοποίηση πινάκων πλήθους συχνοτήτων
#Τέλος αρίθμηση 1
#Απόδωσε αποτέλεσμα: Πίνακας PWM

#Read the file into a list. Since it has N, we can do a preliminary and remove all the N and read it after
import re
f=open("e_coli_fixed.fa", 'r')
seqs = []
for line in f:
    x=re.match(">", line)
    if x == None:
        #seqs.append(line.rstrip().replace('N',''))
        seqs.append(line.rstrip())
        

# %%
def readmultifasta(file):
    import re
    f = open(file, 'r')
    seqs = []
    for line in f:
        x=re.match(">", line)
        if x == None:
            seqs.append(line.rstrip())
    return(seqs)

#seqs = readmultifasta("e_coli_fixed.fa")
#seqs = readmultifasta("staaur.fa")
seqs = readmultifasta("ex2_motifs.fa")


def pwm(sequences):
    nuc = ['A', 'C', 'G', 'T']
    profile=[[0 for i in range(len(sequences[0]))] for j in range(len(nuc))]
    for instance in sequences:
        for j in range(len(instance)):
            residue=instance[j]
            if residue == 'N':
                continue
            profile[nuc.index(residue)][j]+=1
            profile[nuc.index(residue)][j]=float(profile[nuc.index(residue)][j])
    import numpy as np
    pwm = np.array(profile)
    pwm = pwm/len(sequences)
    return(pwm)

mypwm=pwm(seqs)
print(mypwm)
# %%
#PSSM (Position-specific Scoring Matrix): They depends on the PWM matrix and the nucleotide sequence 
# of a background (κατάλοιπα υποβάθρου) motif. We get a matrix with a score. Negative if it doesnt concur
#or positive if it is strongly seen in that position. It compares our motif PWM and the background motif
#R = log2(P[i,j]/Q[i,j])
#1) Δήλωση αλληλουχίας S μήκος n 
#2) Δήλωση πίνακα PWM[4,l]
#3) Δήλωση πίνακα Score[n-l+1]

#Απαρίθμηση 1:
#για θέση i=1 έως i=n-l+1 ανά 1;
#   Δημιούργησε την υποαλληλουχία s=S[i:i+l-1] #μήκους l
#   Απαρίθμηση 2:
#   Για κάθε θέση j=1 έως j=l ανά 1;
#       Score[i]=Score[i]+PWM[s[j],j]
#   Τέλος απαρίθμηση 2
#Τέλος απαρίθμηση 1
#Απόδωσε αποτέλεσμα: Πίνακας Score

#We can transform this table by applying a normalization and log-transformation against nucleotide 
# occurrences from a given sequence.

#We shall first write a function to read the sequence in fasta format
import re
def readfasta(fastafile):
    
    f=open(fastafile, 'r')
    seq = ""
    total = 0
    for line in f:
        x=re.match(">", line)
        if x == None:
            length=len(line)
            total=total+length
            seq=seq+line[0:length-1]
    seq=seq.replace('N','')
    f.close()
    return(seq)

#targetsequence=readfasta("e_coli.fa")
#targetsequence=readfasta("staaur.fa")
targetsequence=readfasta("ex2_motifs.fa")
print((targetsequence))

# %%
#And then write another function to calculate nucleotide frequencies
def nuccomp(sequence):
    import numpy as np
    nucfreq = [0, 0, 0, 0]
    nuc = ['A', 'C', 'G', 'T']
    for i in range(len(nuc)):
        nucfreq[i]=sequence.count(nuc[i])
    nucfreq=np.array(nucfreq)/len(sequence)
    return(nucfreq)


# %%  --------- ALTERNATIVE TO NUCCOMP ---------- # 
def kmers(genomefile, k):
    import re

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
    kmertable = dict(sorted(kmertable.items()))
    return(list(kmertable.values()))
#nucfreqs = kmers("staaur.fa", 1)
nucfreqs = kmers("ex2_motifs.fa", 1)
print(nucfreqs)


# %%
#Τhe next obvious step is to combine the PWM with the Array of the nucleotide composition of the target 
#sequnce and log-transform the resulting table into a PSSM
# Creation of a PSSM
def pssm(pwm, nucfreqs):
    import numpy as np
    import math
    pseudocount=0.01
    pssm=[[0 for i in range(len(pwm[0]))] for j in range(len(nucfreqs))]
    for i in range(len(nucfreqs)):
        pssm[i]=(np.array(pwm[i])+pseudocount)/nucfreqs[i]
    for i in range(len(pssm)):
        for k in range(len(pssm[0])):
            pssm[i][k]=math.log(pssm[i][k])/math.log(2)
    return(np.array(pssm))

mypssm=pssm(mypwm, nucfreqs)
print(mypssm)

# %%
#c. Searching a sequence with a PWM/PSSM or with Hamming Distance 

def pssmSearch(pssm, sequence, threshold):
    nuc = ['A', 'C', 'G', 'T']
    hits = []
    instances = []
    for i in range(len(sequence)-len(pssm[0])):
        instance=sequence[i:i+len(pssm[0])]
        score=0
        for l in range(len(instance)):
            score=score+pssm[nuc.index(instance[l])][l]
        if (score > threshold):
            hits.append(i)
            instances.append(instance) 
    return(hits, instances)

out=pssmSearch(mypssm, targetsequence, 9)
print(out)

# %%
def pssmSearch(pssm, sequence, threshold):

    nuc = ['A', 'C', 'G', 'T']
    hits = []
    instances = []
    x = []
    allscores = [] # for plotting reasons
    # Step 1: Calculation of maximum possible PSSM score
    maxPssm = 0
    for j in range(len(pssm[0])):
        
        maxPssm = maxPssm + max(pssm[:,j])
    print(maxPssm)
    # Step 2: Search
    for i in range(len(sequence)-len(pssm[0])):
        x.append(i)
        instance=sequence[i:i+len(pssm[0])]
        score=0
        for l in range(len(instance)):
            score=score+pssm[nuc.index(instance[l])][l]
        #if (score > threshold*maxPssm): #This doesnt need maxpssm?
        if (score > threshold):
            hits.append(i)
            instances.append(instance)
            allscores.append(score)
    
    return(hits, instances, allscores)

out = pssmSearch(mypssm, targetsequence, 0.95)

print(out)

# %%
#And we can apply Entropy and Information Content calculations on the resulting hits/matches

#Entropia(systhmatos X) = arnhtiko athroisma twn pi*log2(pi), opoy pi = pithanothta kathe endexomenoy toy systhmatos

#-Xamhlh entropia --> Poly plhroforia  (antistrofos analoga san ennoies kai oroi)
#mikros arithmos endexomenwn, h entropia mporei na metablithei sto systhma, alla to megisto symbainei otan
#exoyn thn idia pithanothta (px nomismata)

#H plhroforia orizetai ws h metatroph ths entropias I = Hmeta - Hprin (arnhtiko proshmo). H metatroph
#ths entropias meta apo mia diadikasia

#Gia mia thesh p.x an parathrw mono A, tote Hfin = 1*log(1) + 0 + 0 + 0 = 0, afou pithanothta toy A = 1
#thewreitai pws 0*log0 = 0. Ara h telikh entropia einai 0, ara megisth plhroforia

#H megisth entropia einai Hinit = -4*(0.25log(0.25)) = -2 <------- PANTA H ARXIKH ENTROPIA GIA YPOLOGISMO. 
# Edw einai avevaiothta
# I = Hfin - Hinit = 0 - (-2) = 2.  Edw einai plhroforia. Apo -2 ews 2?.
#--- Pairnw th mesh timh ths plhroforias kai vriskw thn sunolikh plhroforia toy motivoy.
#Motivo me 10 theseis tha exei anamenomena megalyterh se athroisma plhroforia apo motivo me 5 theseis

#Oi times autes --> SeqLogo 
#1) PSSM --> Perissoterh plhroforia (I). Einai pio eidikh
#2) PWM --> ligoterh

def EntropyInformation(pwm):
    
    import numpy as np
    k = pwm.shape[1]

    information = np.zeros([1,k]) #computing the information of each position
    for i in range(k):
        information[0,i] = 2-abs(sum([elem*np.log2(elem) for elem in pwm[:,i] if elem > 0]))
    
    sumInfo = np.sum(information)
    scaledSumInfo = sumInfo/k
    
    return(information, sumInfo, scaledSumInfo)

# %%
pwmsearches = pssmSearch(mypwm, targetsequence, 0.90)[1]
# %%

EntropyInformation(pwm(pwmsearches))
#Prwto meta to array deixnei megisth plhroforia kai to katw einai h mesh plhroforia (9.33.../6)
#Oso dhladh o arithmos apo len(array)
# %%
pssmsearches = pssmSearch(mypssm, targetsequence, 0.90)[1]
print(pssmsearches)

# %%
EntropyInformation(pwm(pssmsearches))
#O pssm deinei kalytera apotelesmata kathws vveltiwnei kai th teleytaia thesh. Megalyterh megisth
#Pio eidikh anazhthsh. Pio plhroforiaka plousia instances



# %%
#Denovo anazhthsh motivou.

#Pairnw to set k-merwn kai metraw poses fores prokuptei to kathe ena.

#Vlepw poses fores emfanizetai kai thelw na dw to anamenomeno poio einai
# (isws edw athroisma olwn kai diairesh dia arithmo motivwn )
# Vriskw pio uperekprosophmena MOTIVA oxi sequence

#Pws tha epileksw apo kathe allhloyxia (s) ena instance (kmer), na ta valw ola se mia syllogh kai na
#ftiaksw sth synexeia to pwm. To kmer tha prepei na einai auto poy tha dinei to kalytero pwm
#estw dhladh oti emfanizetai mono mia fora
#Prepei na doylepsw ana allhloyxia kai na kanw syndiasmoys



#--> Greedy prosegish, na vrw to pio over-represented se kathe Allhloyxia (statistika)
#ta krataw kai kanw PWM alla den eksasfalizei swsto PWM
#Edw problhma emafnizoyn ta repeats allhloyxiwn, den eksetazei th sxesh metaksy allhloyxias

#--> Brute force: Ola ta kmers, apo oles tis allhloyxies kai na kanw olous tous syndiasmoys,
#na kanw meta PWM gia to kathena apo auta kai na krathsw auto me to kalutero. 
# SIGOURH lush alla pairnei poly wra.




# --> GIBBS sampling: Randomized algorithmos. 

#1)We start with a collections of s k-mers, one from each sequence
#2)Build a profile (pwm)
#3)Scan sequences with that profile
#4)Update profile and repeat until the k-mer set is good enough match for the updated profile

#pseudocode

print(len(seqs))
# %%
for seq in seqs:  
    nuc = ['A', 'C', 'G', 'T']
    profile=[[0 for i in range(len(seqs[0]))] for j in range(len(nuc))]
    profile[seq]<-random(k, seq)  
    while distance(profile, sequences)>threshold  
        for seq in sequences:  
        profile[seq]<-max(k, profile, seq)  


1)Dialegw 1 tyxaio apo kathe allhloyxia kai ftiaxnw ena pwm 
2)Xrhsimopoiw pwm gia na sarrwsw kathe allhloyxia apo thn arxh
3)Tha entopisw pws auto poy dialeksa eksarxhs den einai auo poy dinei to kalytero score


#Stamataw otan exw ftasei ena information threshold, h otan paramenei 
# statheros o pwm meta apo merika repeats

#taksh twn xiliwn. Orio sto meso I. Otan to ftanw stamataw


I = 1-1.6



# %% Read for the PWM

import numpy as np
import random, re, math
from collections import Counter

# We start of by the given function to read our fasta file. In this case the 50 (l=100) sequences.
def readmultifasta(file):
    import re
    f = open(file, 'r')
    seqs = []
    for line in f:
        x=re.match(">", line)
        if x == None:
            seqs.append(line.rstrip())
    return(seqs)

# Secondly, let's create another function to take random k-mers from each sequence.
def randomizer(all_seq,len_kmer):
# We use 2 variables, 1 for the fasta (or sequences) and 1 for the length of k-mers we are interested    
    random_seq = [] #Initialize an empty list to add all our kmers
    for seq in all_seq:
        k_index = random.randint(0,len(max(all_seq))-len_kmer) #This way it is provided we will take a different random kmer each time from each sequence
        random_seq.append(seq[k_index:k_index+len_kmer])
      
    return random_seq

# Now let's create a function to generate the position weight matrix, as seen from the course.
def pwm_get(sequences): #The input here is sequences or kmers
    nuc = ['A', 'C', 'G', 'T'] #initialize the order of the nucleotides
    profile=[[0 for i in range(len(sequences[0]))] for j in range(len(nuc))]# and create an empty initial profile
    for instance in sequences:
        
        for j in range(len(instance)):
            residue=instance[j]
            if residue == 'N': #I also added this in case we are dealing with unidentified nucleotides and are not using the single string read method
                continue
            profile[nuc.index(residue)][j]+=1
            profile[nuc.index(residue)][j]=float(profile[nuc.index(residue)][j]) #we now fill in the profile depending in the appearance of a nucleotide
    import numpy as np
    pwm = np.array(profile)
    pwm = pwm/len(sequences)
    return(np.array(pwm))

#  With this function we can now calculate tha information content of a given pwm and consider if we are willing to keep it or not
def EntropyInformation(pwm):
    
    import numpy as np
    k = pwm.shape[1]
    information = np.zeros([1,k]) #computing the information of each position
    for i in range(k):
        information[0,i] = 2-abs(sum([elem*np.log2(elem) for elem in pwm[:,i] if elem > 0]))
    sumInfo = np.sum(information)
    scaledSumInfo = sumInfo/k
    
    return(information, sumInfo, scaledSumInfo)



# %%     @$^@$%^#$%^ TESTING @%##$^ $%^$^ -- With USE of max

def gibbs(sequences, kmer, threshold, conf):
    
    nuc = ['A', 'C', 'G', 'T']
    rand_kmers = randomizer(seqs,kmer)
    pwm = pwm_get(rand_kmers)
    counts = 0
    while True:
        
        info, sum_info, mean_info = EntropyInformation(pwm)
        if mean_info >= threshold and counts > 2000:
            return (rand_kmers, info, sum_info, mean_info, pwm)
        
        #print(mean_info, end = ' ,, counts \n')
        seq_dict = {}
        all_scores = []
        init_index = 0
        maxPssm = 0
 
        
        for j in range(len(pwm[0])):
            maxPssm = maxPssm + max(pwm[:,j])

        for k in sequences:
            score_init = 0    
            
            for l_init in range(len(rand_kmers[score_init])):
                score_init += pwm[nuc.index(rand_kmers[init_index][l_init])][l_init]

            seq_dict[k] = score_init

            for i in range(len(k)-kmer+1):
                instance = k[i:i+kmer]
                score = 0
            
                for l in range(len(instance)):
                    score += pwm[nuc.index(instance[l])][l]
                
                inst_score = (instance,score)
                all_scores.append(inst_score) 
            
            max_instance, max_score = max(all_scores, key = lambda x: x[1])
                
            if max_score >= (seq_dict[k] and maxPssm*conf):
                seq_dict[k] = max_score
                del rand_kmers[init_index]
                rand_kmers.insert(init_index, max_instance)
            
            init_index += 1
            counts += 1
            
            if counts >= 100_000:
                print(f'Stopped after {counts} counts')
                return (rand_kmers, info, sum_info, mean_info, pwm)
        
        pwm = pwm_get(rand_kmers)
    #return (rand_kmers, info, sum_info, mean_info, pwm)

seqs = readmultifasta("ex2_motifs.fa")
motifs = list(range(3,8))


# %%
k_dict = {}
for km in motifs:
    k_dict[km] = gibbs(seqs, km, 1.85, 0.95)
    
    print(f'These are the motifs found: {sorted(dict(Counter(k_dict[km][0])).items(), key = lambda x:x[1],reverse = True)}')
    print(f'The mean information when kmer length is {km}: {k_dict[km][3]} and the PWM is: \n\n {k_dict[km][4]} \n')


# %%
for km in motifs:
    print(km, end = ' <-- kmer repeat \n')
    
    for i in range(10):
        k_dict = {}
        k_dict[km] = gibbs(seqs, km, 1.85, 0.95)
        print(f'These are the motifs found: {sorted(dict(Counter(k_dict[km][0])).items(), key = lambda x:x[1],reverse = True)}')
        print(f'The mean information when kmer length is {km}: {k_dict[km][3]} and the PWM is: \n\n {k_dict[km][4]} \n')
    print('\n\n\n')
# %%                    WITH PSSM
import re
import numpy as np
import math

def readfasta(fastafile):
    
    f=open(fastafile, 'r')
    seq = ""
    total = 0
    for line in f:
        x=re.match(">", line)
        if x == None:
            length=len(line)
            total=total+length
            seq=seq+line[0:length-1]
    seq=seq.replace('N','')
    f.close()
    return(seq)




def nuccomp(sequence):
    nucfreq = [0, 0, 0, 0]
    nuc = ['A', 'C', 'G', 'T']
    for i in range(len(nuc)):
        nucfreq[i]=sequence.count(nuc[i])
    nucfreq=np.array(nucfreq)/len(sequence)
    return(nucfreq)


def pssm_get(pwm, nucfreqs):
    import numpy as np
    import math
    pseudocount=0.01
    pssm=[[0 for i in range(len(pwm[0]))] for j in range(len(nucfreqs))]
    for i in range(len(nucfreqs)):
        pssm[i]=(np.array(pwm[i])+pseudocount)/nucfreqs[i]
    for i in range(len(pssm)):
        for k in range(len(pssm[0])):
            pssm[i][k]=math.log(pssm[i][k])/math.log(2)
    return(np.array(pssm))




def gibbs(sequences, kmer, threshold, nuc_frequency):
    
    nuc = ['A', 'C', 'G', 'T']
    rand_kmers = randomizer(seqs,kmer)
    pwm = pwm_get(rand_kmers)
    pssm = pssm_get(pwm, nuc_frequency)
    counts = 0
    while True:

        info, sum_info, mean_info = EntropyInformation(pssm)
        if mean_info >= threshold and counts > 2000:
            return (rand_kmers, info, sum_info, mean_info, pwm)
        
        #print(mean_info, end = ' ,, counts \n')
        
        seq_dict = {}
        all_scores = []
        init_index = 0
        


        for k in sequences:
            score_init = 0    
            
            for l_init in range(len(rand_kmers[score_init])):
                score_init += pssm[nuc.index(rand_kmers[init_index][l_init])][l_init]

            seq_dict[k] = score_init

            for i in range(len(k)-kmer+1):
                instance = k[i:i+kmer]
                score = 0
            
                for l in range(len(instance)):
                    score += pssm[nuc.index(instance[l])][l]
                
                inst_score = (instance,score)
                all_scores.append(inst_score) 
            
            max_instance, max_score = max(all_scores, key = lambda x: x[1])
                
            if max_score > seq_dict[k]:
                print(seq_dict[k], max_score)
                seq_dict[k] = max_score
                del rand_kmers[init_index]
                rand_kmers.insert(init_index, max_instance)
            
            init_index += 1
            counts += 1
            
            if counts >= 100_000:
                print(f'Stopped after {counts} counts')
                return (rand_kmers, info, sum_info, mean_info, pwm)
        
        pwm = pwm_get(rand_kmers)
        pssm = pssm_get(pwm, nuc_frequency)
        print(pssm)
    #return (rand_kmers, info, sum_info, mean_info, pwm)

seq_string=readfasta("ex2_motifs.fa")
nucs = nuccomp(seq_string)

seqs = readmultifasta("ex2_motifs.fa")
motifs = list(range(4,6))


k_dict = {}
for km in motifs:
    k_dict[km] = gibbs(seqs, km, 1.85, nucs)
    
    print(f'These are the motifs found: {sorted(dict(Counter(k_dict[km][0])).items(), key = lambda x:x[1],reverse = True)}')
    print(f'The mean information when kmer length is {km}: {k_dict[km][3]} and the PWM is: \n\n {k_dict[km][4]} \n')


# %% by replacement and not by max NOT WORKING

def gibbs(sequences, kmer, threshold):
    
    nuc = ['A', 'C', 'G', 'T']
    rand_kmers = randomizer(seqs,kmer)
    pwm = pwm_get(rand_kmers)
    counts = 0
    
    while counts < 50_000:
        info, sum_info, mean_info = EntropyInformation(pwm)

        if mean_info >= threshold:
            print(mean_info, end = ' i am bigger than life \n')
            return (rand_kmers, info, sum_info, mean_info, pwm)
        #print(mean_info, counts, end = ' ,, counts \n')
        seq_dict = {}
        init_index = 0

        for k in sequences:
            
            score_init = 0    
            for l_init in range(len(rand_kmers[score_init])):
                score_init += pwm[nuc.index(rand_kmers[init_index][l_init])][l_init]
            seq_dict[k] = score_init

            for i in range(len(k)-kmer+1):
                instance = k[i:i+kmer]
                score = 0
            
                for l in range(len(instance)):
                    score += pwm[nuc.index(instance[l])][l]
            
                if score >= seq_dict[k]:
                    seq_dict[k] = score
                    del rand_kmers[init_index]
                    rand_kmers.insert(init_index, instance)
            
            init_index += 1
            counts += 1
        print(rand_kmers)
        pwm = pwm_get(rand_kmers)
        
    #return (rand_kmers, info, sum_info, mean_info, pwm)

seqs = readmultifasta("ex2_motifs.fa")
motifs = list(range(3,8))


for i in range(10):
    print(i, end = ' <-- repeat \n')
    k_dict = {}
    for km in motifs:
        k_dict[km] = gibbs(seqs, km, 2)
        print(Counter(k_dict[km][0]))
        print(f'The mean information when kmer length is {km}: {k_dict[km][3]} and the PWM is: \n\n {k_dict[km][4]} \n')
    print('\n\n\n')

