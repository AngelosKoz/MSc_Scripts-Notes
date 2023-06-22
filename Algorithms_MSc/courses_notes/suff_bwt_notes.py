# %%

# Read from NGS --> 10^6, with each beeing > 100bp. 
#We compare those to genomes with more than 10^9.

#Online = Works on patterns. Takes the shorter and searches it in the longer
#Problem here occurs when we have alot of those "query" like from NGS

#Offline = We work on the target sequence (genome) instead of query.
#Idea here is to transform/index the text (sequence-genome) to work faster on it.




## Suffix Trees/Tries
S = 'GAGTAAGTCA'
print(len(S))
# We work with sub-sequences of the target sequence (suffixes - endings). We try to map all suffixes of the target sequence
# We use '$' to show the end of the sequence. After finding them, we sort them alphabetically. We only need to show the importance of A$ vs AA$ (lets sort low->high)
#Then we create a tree from the suffixes with the root on top. We then start adding suffixes with a '$' at the end
# If the prefix of a suffix is already included there we continue from it (removing the prefix of the longer suffix)
#and adding the rest with again a '$' at the end. By the end we will have 4 starting "clades" with A,G,T,C 

#This aproach takes a lot of time and memory, but since we will usually use that item (genome) alot of time then we
#should take a time to create this.


## Suffix Arrays
#To know the position where the sequence exists. We do this by adding the position of the suffix when creating suffixes
#So this array is a list with the positions of suffixes (in the order they sorted them alphabetically).
#By this we can know whenever we get a "hit" in the target sequence when "scanning" with the query


#Basis: Every substring of S is a prefix of some suffix of S
# -Decide if a pattern exists in a sequence
# -Find a pattern in a sequence
# -Find how many times a pattern occurs in a sequence
# -Find the position(s) of a patten in a sequence
# -Find the longest repeat in a sequence
# -Find the longest common pattern between Try this: https://visualgo.net/en/suffixtree


S = 'GAGTAAGTCA'
suffixArray = [9,4,1,5,8,0,2,6,3,7]
#For example 9 = A, 4 = AAGT...., ....
#We can use the suffix array using Binary Search. We take the array and index it. 
# len(suffix_array) = len(sequence)
#Lets say we want to look for AGT
# Because the array is alphabetically sorted, we go in the middle of our sequence and look at what that suffix is. 
# We then know if our motif is before or after. Then we do the same thing again. We need log2(S) steps


# %%
## Burrows Wheeler Transform --- Pop end to front
# Logic of suffix tree (we work on genome), uses endings as suffix arrays, to create a same or smaller length seq 
# that can be used to create smaller genomes for example.

# 2 characters, 1 for start and 1 for end (or just 1 of them) '^' for start '$' or '|' for end

# We create all the rotations of the genome. We take the last letter of the sequence and we move the sequence at the end
#K_kmer[motif] += [K_kmer[motif].pop()]
#After we create all the rotations, we sort them alphabetically and take ONLY the last column. (last element of each)
#Last column (output) is the same length as the original sequence.



#We must know
# 1)How to go from the BWT to the original sequence (BWT inversion)
# We only need to create 1 sequence and we then know the original.
# -The first column is created and enumerated by ordering alphabetically the BWT, and last column is enumerated
# -We know that the sequence will end with '$'. This means in a rotation, that the letter "opposite" of '$' is the first
# -After that, we go to the BWT and find where that letter was found and do the same as above (index 0 index -1)
# -For more than 1 instances, we know that the residues (index [0]) will also be the order they appear in the sequence
#(first last rule) So we enumerate the last column of the BWT

def BWT(seq):
    counter = 0
    bwt_ext = {}
    
    #This will account for sequences that implement their own terminals and/or suffixes
    if seq.isalpha():
        seq += '$'
    
    #A simple for loop to get all the cyclic rotations. The idea is we just displace the last index at the beginning.
    for rotation in range(len(seq)):
        seq = seq[-1] + seq[:-1] #And this is the displacement
        
        #pref_suf = (seq[0], seq[-1])
        #bwt_ext[counter] = pref_suf
        
        bwt_ext[counter] = seq #We use a dummy counter to add all the possible strings
        counter += 1

    bwt_ext = sorted(bwt_ext.items(), key = lambda x: x[1]) #Here we sort based on the string
    
    return [s[-1] for c,s in bwt_ext] #This will return only the last letter of the string

print(BWT('^BANANA|'))
S = 'GAGTAAGTCA'
bwt = (BWT(S))

# %%

def BWTinverse(bwt):
    bwt_inv = {}
    bwt_enu = list(enumerate(bwt)) #Use to retain original position
    inv_sort = sorted(bwt)
    counter = 0
    increment = '1'
    
    print(bwt_enu)
    for i in inv_sort:
    #for i in list(enumerat)
        if not i.isalpha():
            counter += 1
            continue
        
        if counter == len(bwt)-1:
            
            if inv_sort[counter] == inv_sort[counter-1][0]:
                inv_sort[counter] += increment
                break
            
            #This accounts for the last digit beeing different (maybe for words rather than sequences)
            inv_sort[counter] += '1'
            break
        
        if inv_sort[counter] == inv_sort[counter+1]:
            inv_sort[counter] += increment
            increment = str(int(increment)+1)
            counter += 1
            continue 
        inv_sort[counter] += increment
        increment = '1'
        counter += 1
    
    
    print(inv_sort)
    #d['k10'] = d.pop('k1')

#S = '^BANANA|'
S = 'GAGTAAGTCA'

print(BWTinverse(BWT(S)))


# %%
# 2)How to search for substrings (pattern search) within BWT.



a[10]