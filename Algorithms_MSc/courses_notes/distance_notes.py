# Edit Distance - Humming Distance : Epitrepei mono antikatastaseis (dhladh posa diaferoyn)
#Comparison is done for each of the nucleotides without allowing for re-arrangement, deletion or insertion

#+1 match // -2 mismatch // gap = -2
#Aristera --> Liga kena me polla mismatches (diagwnies metatopiseis)
##Deksia perissotera kena (orizonties kai kathetes metatopiseis)


#Dynamic Programming
def SimpleFibonacci(N):
	fib=[]
	fib.append(0)
	fib.append(1)
	for i in range(2,N+1):
		fib.append(fib[i-1]+fib[i-2])
	return fib[i]

#Recursion in Fibonacci
def RecursiveFibonacci(N):
    fib = {}
    if N in {0,1}:
        fib[N] = N
    if N > 1:
        fib[N] = RecursiveFibonacci(N-2) + RecursiveFibonacci(N-1)
    return fib[N]

#We can deploy some sort of dynamic programming approach to speed up the recursion. 
#We can basically check if the value to be calculated is already calculated and then skip it calculation step.

def DPFibonacci(N):
    fib = {}
    if N in {0,1}:
        fib[N] = N
    if N > 1:
        if (N not in fib):
            fib[N] = RecursiveFibonacci(N-2) + RecursiveFibonacci(N-1)
    return fib[N]
import time
start_time = time.time()
myFib = SimpleFibonacci(10)
print(myFib)
print("--- %s seconds ---" % (time.time() - start_time))
#
start_time = time.time()
myFib = RecursiveFibonacci(10)
print(myFib)
print("--- %s seconds ---" % (time.time() - start_time))
#
start_time = time.time()
myFib = DPFibonacci(10)
print(myFib)
print("--- %s seconds ---" % (time.time() - start_time))


#NGS - Πολλες μικρες αλληλουχιες εναντι μιας πολυ μεγαλυτερης αλληλουχιες (ολοκληρο χρωμοσωμα ή γονιδιωμα)
#Αλληλουχια S, Μικρότερη P.

#Allignment Needleman-Wunsch : Καλός για μικρες αλληλουχίες, χάνει για μεγάλες. 
#-Σύγκριση 2 αλληλουχιων
#-1 εναντι με 10
#-ΟΧΙ να βρουμε μια μικρη σε μια μεγαλυτερη

#We use pattern matching. Trying to find for example where P is found on the bigger S sequence (though
# we can still find it many times and not just one).
#In pattern matching we usually allow for mismatches. In this case we are talking about IDENTICAL 
# pattern matching.
#The difficulty here comes from the length of the seqs



#1)Start from sequence[0] and loop one residue at a time (i=0)
#2)Take a substring from sequence equal to pattern
#3)Start from pattern[0] and compare to sequence[i] (j=0)
#4)If there is a mismatch, exit and go to 2, take next i
#5)If there are no mismatches and j reaches the size of length of pattern, report a full match, go to 2 and take next i

# %%
# Naive Pattern Searching algorithm
def naivePatternSearch(pattern, sequence):
    p = len(pattern)
    s = len(sequence)
    print(s)
    no = 0
 
    # We slide pattern one residue/character at a time
    for i in range(s - p + 1):
        no += 1
        j = 0
        
        # for each pairing of pattern to sequence we check characters starting from the beginning
        while(j < p):
            if (sequence[i + j] != pattern[j]): 
                break
            j += 1
 
        if (j == p):
            print("Pattern found at index ", i)
    return("Number of steps taken=", no)
 
naivePatternSearch('abba', 'IhateabbaandIdontlikealibabba')


# %%

#There are actually two points that we should consider. 
# One is that in a case of a long pattern that matches the fist p-1 positions but not in the final 
# one, we may be doing a lot of comparisons we can do away with. 
# The other is that we may be able to speed up the process provided that the pattern has some 
# particular properties, by sliding not 1 character at a time, but more.
# The condition is for the pattern to be variable because when all characters of the pattern are 
# different, we can slide the pattern by more than 1. This is because when a mismatch occurs 
# after j matches, we know that the first character of pattern will not match the j matched 
# characters because all characters of pattern are different. So we can always slide the pattern 
# not by 1 but by j without missing any valid shifts.


# Optimized Naive Search 
def optNaivePatternSearch(pattern, sequence):
    p = len(pattern)
    s = len(sequence)
    i = 0
    no = 0
  
    while i <= s-p:
        no+=1
        # For current index i, check for pattern match
        for j in range(p):
            if sequence[i+j] != pattern[j]:
                break
            j += 1
        if j==p:    # if pat[0...M-1] = txt[i,i+1,...i+M-1]
            print("Pattern found at index " + str(i))
            i += p
        elif j==0:
            i += 1
        else:
            i += j    # slide the pattern by j
    return("Number of steps taken=", no)

optNaivePatternSearch('abba', 'IhateabbaandIdontlikealibabba')


# %%

# KMP (Knuth-Morris-Pratt Algorithm  -- Pattern Sliding Algorithm). This algorithm, as well as BM,RK..
# firstly do pattern transformation. It means they check on the pattern we are looking for

# Here we create a prefix array --> LPS (Longest-Prefix-Suffix).
# Prefix = Sub-sequence that the pattern starts with
# Suffix = Sub-sequence that the pattern ends with
# Its an array with len(LPS) == len(pattern)
# Each number in the array is the length of the biggest prefix that is also a suffix (both start + end)
# When looking at prefixes, we dont want to include ourself, thus the first number of the array is 0.
# For example in ACACA, the first prefix is A. But we dont want that. The next is AC. Pref = A, Suff = C
# Because A != C, again we input 0. 
# Next is ACA. Pref = A, AC // Suf = A, CA. Here A and A match, so we input the length (which is 1)
# Next is ACAC. Pref = A, AC, ACA // Suf = C, AC, CAC. Here AC match, length = 2, we input 2.
# And so on.

#Consider a prefix r such as that r=ACACA

#All proper(=not containing itself) prefixes of r are:
#A, AC, ACA and ACAC
#All proper suffixes of r are: CACA, ACA, CA, A
#The longest substring that is both a prefix and a suffix of r is ACA. It has a length=3.
#We assign 3 at the position where ACACA terminates in the prefix array.

def kmpPrefixArray(pattern):
    kmp=[]
    for i in range(len(pattern)+1):
        subpattern=pattern[:i]
        #print(subpattern, end = ' subpat<<< \n')
        maxlen=0
        for j in range(len(subpattern)):
            print(f' {subpattern[:j]} and {subpattern[-j:]}')
            #print(len(subpattern[j:]))
            if (subpattern[:j]==subpattern[-j:]) and (len(subpattern[j:])>maxlen):
                maxlen=len(subpattern[:j])
        kmp.append(maxlen)
    kmp.pop(0)
    return(kmp)

print(kmpPrefixArray('ACACA'))

#print(kmpPrefixArray('GAGCCGAGCCGAGTCTG'))

#print(kmpPrefixArray('GATCGCGACGTTCAGCT'))



# %%

#Now we will use the prefix array to search the sequence

#The idea is to use the LPS array to speed up the process by directing us at which characters 
# between string and pattern to match every time.

#>The algorithm may be outligned as follows:<#

#1)Start from the beginning i of sequence and j of pattern.
#2)For as long as there is a match continue matching.
#3)Upon hitting a mismatch, check if j=0. If yes, move one character in sequence.
#4)In j is not equal to 0 (which means we are in the middle of the pattern), look at the value of LPS[j-1].
#5)Consider j (the position in the pattern) to be equal to the value of LPS[j-1] (e.g. if LPS[j-1]==2 
# then move to the third position of the pattern)
#6)If j reaches the length of the pattern, then a complete match has been reached. Return exact match.
#7)After a match is found move the pattern as many places as the last value of LPS indicates.


#>The key is this: In every case (either a match or a mismatch):<#

#1)we take the last matching residue
#2)we look up its LPS value
#3)we move the pattern whose index matches that LPS value and place it at the position of the mismatch

# In case of matches we move as many positions as the last match at the LPS tells us. In case of mismatch 
# we move the pattern according to the index-value of the previous index of the LPS of the mismatch.
# This means we will place that given index, below the mismatch

def kmpSearch(pattern, sequence):
    successes=[]
    kmppattern=kmpPrefixArray(pattern)
    s=len(sequence)
    p=len(pattern)
    i=0 
    j=0
    while i < s:
        if (sequence[i] == pattern[j]):
            i += 1
            j += 1
        else: 
            if j != 0:
                j = kmppattern[j-1]
            else:
                i += 1
        if j == p:
            successes.append(i-j)
            j = kmppattern[j-1]
    return(successes)

kmpSearch('abba', 'IhateabbaandIdontlikealibabba')

#Other approaches : Boyer - Moore, ..


# %%

# FASTA Algorithm.

# The logic is :
#1) Get two seqeunces (Q, S)
#2) Choose a length k for kmers
#3) Map all kmers in Q, S

#We look for identical k-mers and see if we can extend our search around those.

# We create a table for L = longer-pattern, K = smaller-pattern, and length-of-kmers = n

#Identity/Dot-matrix. We have the longer pattern as axis_X and the smaller as axis_y. Whenever we get
# an identical residue, we place a dot. The fasta algorithm then looks at diagonals (l+k-1 diagonals) 
# and as we move diagonals, we score the number of CONTINUOUS dots in each. 

#L = GCATCGCCATCG
#K = GCATCGCGC
#1) We start from the smaller sequence and take the first k-mer (i=1), K[i] = GC. We start with i=1 not 0
#2) We search K[i] in L (since we have the positions in the array). For example position 6 (p=6)
#3) Index d is i-p. The first position we started on K - the position we found it on L. 
#4) We go to the score of those diagonals (S dict) and add 1. S[d] = S[d] + 1.
#5) Repeat for each 2-mer of K.


# We start the search with a k depending on the lengths of sequences. 
#After that we look for the optimal length of matching (small with big). We have scores where again,
# we take another parameter as a type of threshold (s > l1) that whichever score is lower than that we dont 
# take it into account. (We try to remove low scores). So S >= l1
# Lastly we want to add a gap score, which we want to be lower than a gap threshold. (g < l2). 
# The "g" is the distance between the diagonal (S) indexes.


# %% Index the kmers of a sequence
def indexKmers(sequence, k):
    import itertools
    import re
    alphabet = ['A', 'C' , 'G', 'T']
    kmers = [''.join(p) for p in list(itertools.product(alphabet, repeat=k))]
    positions = {}
    for kmer in kmers:
        regex = '(?='+kmer+')'
        positions[kmer] = [m.start() for m in re.finditer(regex, sequence)]
    return(positions)

#kmersTEST = indexKmers('ACACACCGTCACACACGTTACACA', 2)
#print(kmersTEST)

# %% Exercise
    
import time

def fasta(L, K, kmer, lthr, gthr, final_score):
    K_kmer = indexKmers(K,kmer) #We start by getting the index of each kmer in our smaller sequence
    successes = {}
    clean_succ = {}
    gene_seqs = {}

    
    for i in L: #Read each line of the file
        
        if ">" in str(i): #Skip lines without sequences
            curr_gene = str(i.rstrip('\n'))
            continue
        
        diag = {}
        start_points_dummy = {}
        
        L_kmer = indexKmers(i.rstrip('\n'),kmer) #Get the index of each kmer in each line (with sequences) of our file
                
        for j in range(len(K) - kmer + 1): #A simple sliding loop with a window of 1
            
            motif = K[j:j+kmer] #Assign each k-mer
                    
            L_match = L_kmer[motif] #And use it to get the indexes of that k-mer in our bigger sequence
            if not L_match: #We move on to the next motif if it is not found in the bigger sequence
                continue
            
            K_match = K_kmer[motif] #If we find it then we get the indexes of that k-mer
    
            pos = int(K_match[0]) #We will use the first element of our list since we will displace it later on
            count = 0
            
            while count != len(L_match): #We will repeat our loop for every element in the bigger sequence index list
                
                d = pos - int(L_match[count]) #We calculate the diagonal index
                if d not in diag.keys():
                    diag[d] = 1 #initialize it if it does not exist
                    start_points_dummy[d] = pos
                    count += 1 
                    continue
                                
                if pos < start_points_dummy[d]:
                    start_points_dummy[d] = pos
                
                diag[d] +=1 #If it exists, we just add 1
                count += 1    

     
            K_kmer[motif] += [K_kmer[motif].pop(0)]#Since we are using the first element always as the index, we can just move the element at the end, so we can keep the order by the end.
        
        #Here we only need to check if the max value is above the threshold
        if max(diag.values()) >= lthr:
            #Then we keep only the values that are above the score threshold
            diag_clean = dict(sorted({k:v for k,v in diag.items() if v >= lthr}.items(), key = lambda x : x[0], reverse=True))
            #This line can be used if we are also interested in the relative position of score assignment inside the loop (enumerate)
            #diag_clean = dict({k:v for k,v in enumerate(diag.items()) if v[1] >= lthr}.items(), key = lambda x : x[0])
            
            
            #Here are some dummy dictionaries to check if my chaotic code below works for different cases. Weirdly enough it seems to work.
            #diag_clean = {-2: 43, 3: 40, 4: 29, 17: 1685, 20:15, 50:100}
            #diag_clean = {2: 43, 3: 40, 4: 29, 17: 1685, 20:15, 50:100}
            #diag_clean = {2: 43, 3: 40, 4: 29, 17: 1685, 20:15}
            
            #print(diag_clean, end = f' for {curr_gene} \n') #Remove comment to check "real-time" values that exceed the threshold and for which gene
            
            
            dist = 1 
            dist_dum = 0
            score_sum = 0
            successes[curr_gene] = {} #This dictionary will hold all of the hits
            clean_succ[curr_gene] = {} #While this will only keep those that pass our final score threshold
            
            
            #And now the madness begins. What i tried to do here, is to check consecutive elements for the gap threshold.
            #After alot of experimenting i got it to work and so far seems to account for most cases i checked.
            #Of course, this part is entirely optional and can be removed, since this is only done to try and keep all of the successes
            # and not only the scores that exceed 50% of the length of the query
            while dist != len(diag_clean):
                
                #First we check if the consecutive elements have an index gap < threshold and save the score
                #Every time we check the next element (since we start with dist 1) and compare to previous
                if abs(list(diag_clean.keys())[dist] - list(diag_clean.keys())[dist-1]) <= gthr:
                    if score_sum == 0:
                        score_sum = list(diag_clean.values())[dist] + list(diag_clean.values())[dist-1]
                    
                    else:
                        score_sum += list(diag_clean.values())[dist]
                    
                    dist += 1
                    continue
                
                #This accounts for cases where gap threshold is not met and takes the previous element (since we start with position and not python indexing)
                #The dist_dum helps keep count of the consecutive elements we will add to the dictionary so this way we will have an overview of indexes and scores
                if score_sum == 0:    
                    score_sum = list(diag_clean.values())[dist-1]
                    successes[curr_gene][score_sum] = list(diag_clean.items())[dist_dum:dist]
                    dist_dum = dist
                    score_sum = 0
                    dist += 1
                    continue  
            
                successes[curr_gene][score_sum] = list(diag_clean.items())[dist_dum:dist]
                
                dist_dum = dist
                score_sum = 0
                dist += 1
                
            #To account for last element(s) from experimenting i found that this was necessary.
            
            #Since the loop will stop when reaching the last element, if no score was assigned we add it here
            if not score_sum:
                score_sum = list(diag_clean.values())[dist-1]
                successes[curr_gene][score_sum] = list(diag_clean.items())[dist_dum:dist]
            
            #And this part is necessary in case the last 2 elements had a gap <= gap_threshold
            if score_sum:
                successes[curr_gene][score_sum] = list(diag_clean.items())[dist_dum:dist]
                
            #And finally here we keep all of the scores that exceed 50% len(query).
            final_pass = [scsum for scsum in list(successes[curr_gene].keys()) if float(scsum) >= float(len(K)*final_score)]
            
            if final_pass:
                clean_succ[curr_gene]['Start'] = []
                
                for ele in final_pass:
                    clean_succ[curr_gene][ele] = successes[curr_gene][ele]
                    min_point = [i1 for i1,i2 in successes[curr_gene][ele]]
                    
                    if not clean_succ[curr_gene]['Start']:
                        min_point = [2,3,4]
                        for pt in min_point:        
                            clean_succ[curr_gene]['Start'].append(start_points_dummy[pt])
                
                for k,v in clean_succ.items():
                        counter = 0
                        gene_seqs[k] = []
        
                        for itm in list(v.items()):
                            
                            if type(itm[0]) == str:
                                continue
                            
                            for stp in v[itm[0]]:
                                        
                                total_l = v['Start'][counter] + stp[1]
                                while True:
                                    
                                    dummy_seq = K[v['Start'][counter]:total_l]
                                    
                                    if dummy_seq in i.rstrip('\n'):
                                        total_l += 1
                                        if K[v['Start'][counter]:total_l] in i.rstrip('\n'):
                                            continue
                                        else:
                                            gene_seqs[k].append(dummy_seq)
                                            break
                                    
                                    if dummy_seq not in i.rstrip('\n'):
                                        total_l -= 1
                                        if K[v['Start'][counter]:total_l] not in i.rstrip('\n'):
                                            continue
                                        else:
                                            gene_seqs[k].append(K[v['Start'][counter]:total_l])
                                            break
                                counter += 1       
                                

    
    #print({s1:successes[s1][s3] for s1,s2 in successes.items() for s3,s4 in s2.items() if float(s3) >= float(len(K)*final_score)})
    return {s1:s2 for s1,s2 in clean_succ.items() if s2}, successes, gene_seqs

            


pattern = "GTATATAACCTAAAAAGAAAATTTGATTAACGAATTCAGAACCCACACGCTGAGACAAGGGGCCCCTGATTATACACCTCTGGAAATGGAAATGGTTTGAATGGGCAAAGCCAGGAAGATATTGAGATTTTAACATATATATTATGTGACATTTGCCCTTCTTTTAGATGGTAGATAATGCTTATTTATGTAGCACGGGTCATTCATCATGGATTTGCGGAGTGGAAGACCAAATGATCATTAACCTGCTACCTGCTTTCAGAATCCTAATGATTTCAAAATGGATAACTTTACCTGTCTCTCAAATTCAGGTCTATTATCTCTCCACAAGATGCAAGCATCAATGTTGGCACCACTTTCGATATTGGGCTCACTCAACATGCTCATAACACTTAATAGAATTTTTTCTACACTTTGCACTGGCGACCATCTTTCTTCCGCTAATTCGTACATGTTAGGATCATCACCAGGGGAGTGTAGAATGGATATGCACACTTCCCCATTTGGATAAATATTTGGATGTAGTATGCTGGGTGTGAAAGTAAGTTTAGGTGGAGATAACGGATAGTCTTTAGGAAACTCTAGCTTAGCATTAAAAACACCATCAGCGTATGGCGTATCTGGAGGCCCTTGAATTAGGCAGTCCCAAATGAATATGTTATTCTCCGATTTGGGACCAGCCACTATACCAGGTGGAGAATCTTTAATTAACTGTTGAAGCTCCTTGAGGAGACGTTTCTGAGCGGTTTTCGACATGCTATGCCCTTCCAAATTACACTATTACTAGGGAAGTTCCTTTTAGGATAATCTCCTTCGTACGCTAAACGCCAGAGTTTACTTTGGCGCTTTTCGAGCTCTTGGTTTAGAGGCTGTAATCTCGTTTTCGGGTAATGGCGAAAAGGAGTTATAAGAAGGATCTCGAGACAATAAGCTGCTGCATCTTCGTGAGGTAGATGCGATGAGGCGCCTTGTTTTAAACTTGAAACAGTCTGGTAAGTTCTTCAAGCTCTATCAGAGATGATAATATTTAATGGGAACAAATATGCGTGTGCATCGTGCATCAGAGGGCATCGCTCTTCAACATGCAGGCATTCTCACCGAATGCTAATTAAGGTGAGAACTAGAGGAAGACCTTCACCCATGGCTATCAGAGACGCTATTTTAGTAGACTCTACATCGCAAAGTACAGAATACGAAAATGGTGCACAAATTGAAGGTGACTGCTGCAGCGCAATGAATCAACAGCCAATACTGTTTGTACGTGCCTCTGCTGTTAGAAAGGCAAGAATGATAAACGGAAAATTGCATATATTAATGGAGGAAGGTTTCACTGCTCATGAGCCTAAAGATATTAGCACATTTACCGATGATGGTAACAAATATATCACCGAAACGGAGTTTCTTAGGAAACACTCTCCCAAAGCTCCGGCAACAGGAACAATATCTCCGGACTCTACCAAGTCATCTTCTTCAAGTGAGAAGAAAGAGCGAAGCCGGCTTCAACAGGAGCCTATACGACATTTTTCAAACTGCTGTAAGAAAGATAAATCACAGAACCCAGCTTCTAATGGCAAGACGAACAAGGCACCGTCTGATGACATATTTACGCCATACGGCTCCTTGGAATCTACGTCCGCTTTTAACGATATTTTACAAGAAAACTACAATAGTTCTGTTCCTGGTGCGCATGACAGTTCAGAAACACTCACCCCACAAAGTACAACAACGATTGCTGCTCCTCATTCAAGCGACGTTGCTTCGAAAGTTGAAGTCCTGACTCATAAGGGCATTTTTTTAAGCACGCAGTGCTCTTGTGAAGATGAAAGCTGCCCATGTGTTAATTGTCTAATCCATAGAAGCGAAGAGGAACTGAATTCTTATATTCAACAAAGTGGTGTTCGGCCTTAATTAACTTAAGCAACATCCAATCTAAGATGAAGGATTATCTCCGGAGGATTGTAAATGCCCTGACAAGGACCGCATATGCCTCCTGGGGATAACCTTACTTGTGATGGAT"
seq = open('rep_cerevisiae.fa')

start_time = time.time()

clean_s, all_s, seq_dict = fasta(seq, pattern, kmer = 5, lthr = 20, gthr = 3, final_score = 0.5) #50%
#print(clean_s, all_s, seq_dict)
print(seq_dict)


print("--- %s seconds ---" % (time.time() - start_time))



# %%
print(len(seq_dict['>YMR021C'][0]))


# %%

#'Start': [12, 22, 82]

pattern = "GTATATAACCTAAAAAGAAAATTTGATTAACGAATTCAGAACCCACACGCTGAGACAAGGGGCCCCTGATTATACACCTCTGGAAATGGAAATGGTTTGAATGGGCAAAGCCAGGAAGATATTGAGATTTTAACATATATATTATGTGACATTTGCCCTTCTTTTAGATGGTAGATAATGCTTATTTATGTAGCACGGGTCATTCATCATGGATTTGCGGAGTGGAAGACCAAATGATCATTAACCTGCTACCTGCTTTCAGAATCCTAATGATTTCAAAATGGATAACTTTACCTGTCTCTCAAATTCAGGTCTATTATCTCTCCACAAGATGCAAGCATCAATGTTGGCACCACTTTCGATATTGGGCTCACTCAACATGCTCATAACACTTAATAGAATTTTTTCTACACTTTGCACTGGCGACCATCTTTCTTCCGCTAATTCGTACATGTTAGGATCATCACCAGGGGAGTGTAGAATGGATATGCACACTTCCCCATTTGGATAAATATTTGGATGTAGTATGCTGGGTGTGAAAGTAAGTTTAGGTGGAGATAACGGATAGTCTTTAGGAAACTCTAGCTTAGCATTAAAAACACCATCAGCGTATGGCGTATCTGGAGGCCCTTGAATTAGGCAGTCCCAAATGAATATGTTATTCTCCGATTTGGGACCAGCCACTATACCAGGTGGAGAATCTTTAATTAACTGTTGAAGCTCCTTGAGGAGACGTTTCTGAGCGGTTTTCGACATGCTATGCCCTTCCAAATTACACTATTACTAGGGAAGTTCCTTTTAGGATAATCTCCTTCGTACGCTAAACGCCAGAGTTTACTTTGGCGCTTTTCGAGCTCTTGGTTTAGAGGCTGTAATCTCGTTTTCGGGTAATGGCGAAAAGGAGTTATAAGAAGGATCTCGAGACAATAAGCTGCTGCATCTTCGTGAGGTAGATGCGATGAGGCGCCTTGTTTTAAACTTGAAACAGTCTGGTAAGTTCTTCAAGCTCTATCAGAGATGATAATATTTAATGGGAACAAATATGCGTGTGCATCGTGCATCAGAGGGCATCGCTCTTCAACATGCAGGCATTCTCACCGAATGCTAATTAAGGTGAGAACTAGAGGAAGACCTTCACCCATGGCTATCAGAGACGCTATTTTAGTAGACTCTACATCGCAAAGTACAGAATACGAAAATGGTGCACAAATTGAAGGTGACTGCTGCAGCGCAATGAATCAACAGCCAATACTGTTTGTACGTGCCTCTGCTGTTAGAAAGGCAAGAATGATAAACGGAAAATTGCATATATTAATGGAGGAAGGTTTCACTGCTCATGAGCCTAAAGATATTAGCACATTTACCGATGATGGTAACAAATATATCACCGAAACGGAGTTTCTTAGGAAACACTCTCCCAAAGCTCCGGCAACAGGAACAATATCTCCGGACTCTACCAAGTCATCTTCTTCAAGTGAGAAGAAAGAGCGAAGCCGGCTTCAACAGGAGCCTATACGACATTTTTCAAACTGCTGTAAGAAAGATAAATCACAGAACCCAGCTTCTAATGGCAAGACGAACAAGGCACCGTCTGATGACATATTTACGCCATACGGCTCCTTGGAATCTACGTCCGCTTTTAACGATATTTTACAAGAAAACTACAATAGTTCTGTTCCTGGTGCGCATGACAGTTCAGAAACACTCACCCCACAAAGTACAACAACGATTGCTGCTCCTCATTCAAGCGACGTTGCTTCGAAAGTTGAAGTCCTGACTCATAAGGGCATTTTTTTAAGCACGCAGTGCTCTTGTGAAGATGAAAGCTGCCCATGTGTTAATTGTCTAATCCATAGAAGCGAAGAGGAACTGAATTCTTATATTCAACAAAGTGGTGTTCGGCCTTAATTAACTTAAGCAACATCCAATCTAAGATGAAGGATTATCTCCGGAGGATTGTAAATGCCCTGACAAGGACCGCATATGCCTCCTGGGGATAACCTTACTTGTGATGGAT"

check = "GTATATAAATAAAAGAAAATTGATTAACGCTTCGGTACAGAGAGGTGAGAGTAGAGGGGCCCCTGATTATACACCAACGAAATGGAAATGGTTTGAATTCCTGAAAAGCCAGGAAGACATTGAGATTTTAACATATTTATTATGTGACATTTGCTTTCTTTTAGATGGTAGATATGCTTATTTATGTATATAGAGAACAGTTAAAAGGAAGACCAAATGATCATTAACCTGCTACCTGCTTTCAGAATCCTAATGATTTCAAAATGGATAACTTTACCTGTCTCTCAAATTCAGGTCTATTATCTCTCCACAAGATGCAAGCATCAATGTTGGCACCACTTTCGATATTGGGCTCACTCAACATGCTCATAACACTTAATAGAATTTTTTCTACACTTTGCACTGGCGACCATCTTTCTTCCGCTAATTCGTACATGTTAGGATCATCACCAGGGGAGTGTAGAATGGATATGCACACTTCCCCATTTGGATAAATATTTGGATGTAGTATGCTGGGTGTGAAAGTAAGTTTAGGTGGAGATAACGGATAGTCTTTAGGAAACTCTAGCTTAGCATTAAAAACACCATCAGCGTATGGCGTATCTGGAGGCCCTTGAATTAGGCAGTCCCAAATGAATATGTTATTCTCCGATTTGGGACCAGCCACTATACCAGGTGGAGAATCTTTAATTAACTGTTGAAGCTCCTTGAGGAGACGTTTCTGAGCGGTTTTCGACATGCTATGCCCTTCCAAATTACACTATTACTAGGGAAGTTCCTTTTAGGATAATCTCCTTCGTACGCTAAACGCCAGAGTTTACTTTGGCGCTTTTCGAGCTCTTGGTTTAGAGGCTGTAATCTCGTTTTCGGGTAATGGCGAAAAGGAGTTATAAGAAGGATCTCGAGACAATAAGCTGCTGCATCTTCGTGAGGTAGATGCGATGAGGCGCCTTGTTTTAAACTTGAAACAGTCTGGTAAGTTCTTCAAGCTCTATCAGAGATGATAATATTTAATGGGAACAAATATGCGTGTGCATCGTGCATCAGAGGGCATCGCTCTTCAACATGCAGGCATTCTCACCGAATGCTAATTAAGGTGAGAACTAGAGGAAGACCTTCACCCATGGCTATCAGAGACGCTATTTTAGTAGACTCTACATCGCAAAGTACAGAATACGAAAATGGTGCACAAATTGAAGGTGACTGCTGCAGCGCAATGAATCAACAGCCAATACTGTTTGTACGTGCCTCTGCTGTTAGAAAGGCAAGAATGATAAACGGAAAATTGCATATATTAATGGAGGAAGGTTTCACTGCTCATGAGCCTAAAGATATTAGCACATTTACCGATGATGGTAACAAATATATCACCGAAACGGAGTTTCTTAGGAAACACTCTCCCAAAGCTCCGGCAACAGGAACAATATCTCCGGACTCTACCAAGTCATCTTCTTCAAGTGAGAAGAAAGAGCGAAGCCGGCTTCAACAGGAGCCTATACGACATTTTTCAAACTGCTGTAAGAAAGATAAATCACAGAACCCAGCTTCTAATGGCAAGACGAACAAGGCACCGTCTGATGACATATTTACGCCATACGGCTCCTTGGAATCTACGTCCGCTTTTAACGATATTTTACAAGAAAACTACAATAGTTCTGTTCCTGGTGCGCATGACAGTTCAGAAACACTCACCCCACAAAGTACAACAACGATTGCTGCTCCTCATTCAAGCGACGTTGCTTCGAAAGTTGAAGTCCTGACTCATAAGGGCATTTTTTTAAGCACGCAGTGCTCTTGTGAAGATGAAAGCTGCCCATGTGTTAATTGTCTAATCCATAGAAGCGAAGAGGAACTGAATTCTTATATTCAACAAAGTGGTGTTCCTTTAACCAATATTGGTGAAGCTCAAATTACTGATAAGATGATGGATTATTTGGATGATTGTAAATGCACTGACAAGGAATGCATATGTCCTCCGGATAACTGCACTTGTGATGGAT"

pattern1 = pattern[223:1900]


print(len(pattern1))
#pattern1 = pattern[12:32]

pattern1 in check

# %% SIMPLE DRAFT
def fasta(L, K, kmer, lthr, gthr, final_score):
    K_kmer = indexKmers(K,kmer) #We start by getting the index of each kmer in our smaller sequence
    successes = {}
    clean_succ = {}
    
    for i in L: #Read each line of the file
        
        if ">" in str(i): #Skip lines without sequences
            curr_gene = str(i.rstrip('\n'))
            continue
        
        diag = {}
        
        L_kmer = indexKmers(i.rstrip('\n'),kmer) #Get the index of each kmer in each line (with sequences) of our file
                
        for j in range(len(K) - kmer + 1): #A simple sliding loop with a window of 1
            
            motif = K[j:j+kmer] #Assign each k-mer
                    
            L_match = L_kmer[motif] #And use it to get the indexes of that k-mer in our bigger sequence
            if not L_match: #We move on to the next motif if it is not found in the bigger sequence
                continue
            
            K_match = K_kmer[motif] #If we find it then we get the indexes of that k-mer
    
            pos = int(K_match[0]) #We will use the first element of our list since we will displace it later on
            count = 0
            
            while count != len(L_match): #We will repeat our loop for every element in the bigger sequence index list
                
                d = pos - int(L_match[count]) #We calculate the diagonal index
                if d not in diag.keys():
                    diag[d] = 1 #initialize it if it does not exist
                    count += 1 
                    continue
                
                diag[d] +=1 #If it exists, we just add 1
                count += 1    
            
            K_kmer[motif] += [K_kmer[motif].pop(0)]#Since we are using the first element always as the index, we can just move the element at the end, so we can keep the order by the end.
        
        #Here we only need to check if the max value is above the threshold
        if max(diag.values()) >= lthr:
            #Then we keep only the values that are above the score threshold
            diag_clean = dict(sorted({k:v for k,v in diag.items() if v >= lthr}.items(), key = lambda x : x[0], reverse=True))
            #This line can be used if we are also interested in the relative position of score assignment inside the loop (enumerate)
            #diag_clean = dict({k:v for k,v in enumerate(diag.items()) if v[1] >= lthr}.items(), key = lambda x : x[0])
            
            
            #Here are some dummy dictionaries to check if my chaotic code below works for different cases. Weirdly enough it seems to work.
            #diag_clean = {-2: 43, 3: 40, 4: 29, 17: 1685, 20:15, 50:100}
            #diag_clean = {2: 43, 3: 40, 4: 29, 17: 1685, 20:15, 50:100}
            #diag_clean = {2: 43, 3: 40, 4: 29, 17: 1685, 20:15}
            
            #print(diag_clean, end = f' for {curr_gene} \n') #Remove comment to check "real-time" values that exceed the threshold and for which gene
            
            dist = 1 
            dist_dum = 0
            score_sum = 0
            successes[curr_gene] = {} #This dictionary will hold all of the hits
            clean_succ[curr_gene] = {} #While this will only keep those that pass our final score threshold
            
            #And now the madness begins. What i tried to do here, is to check consecutive elements for the gap threshold.
            #After alot of experimenting i got it to work and so far seems to account for most cases i checked.
            #Of course, this part is entirely optional and can be removed, since this is only done to try and keep all of the successes
            # and not only the scores that exceed 50% of the length of the query
            while dist != len(diag_clean):
                
                #First we check if the consecutive elements have an index gap < threshold and save the score
                #Every time we check the next element (since we start with dist 1) and compare to previous
                if abs(list(diag_clean.keys())[dist] - list(diag_clean.keys())[dist-1]) <= gthr:
                    if score_sum == 0:
                        score_sum = list(diag_clean.values())[dist] + list(diag_clean.values())[dist-1]
                    
                    else:
                        score_sum += list(diag_clean.values())[dist]
                    
                    dist += 1
                    continue
                
                #This accounts for cases where gap threshold is not met and takes the previous element (since we start with position and not python indexing)
                #The dist_dum helps keep count of the consecutive elements we will add to the dictionary so this way we will have an overview of indexes and scores
                if score_sum == 0:    
                    score_sum = list(diag_clean.values())[dist-1]
                    successes[curr_gene][score_sum] = list(diag_clean.items())[dist_dum:dist]
                    dist_dum = dist
                    score_sum = 0
                    dist += 1
                    continue  
            
                successes[curr_gene][score_sum] = list(diag_clean.items())[dist_dum:dist]
                
                dist_dum = dist
                score_sum = 0
                dist += 1
                
            #To account for last element(s) from experimenting i found that this was necessary.
            
            #Since the loop will stop when reaching the last element, if no score was assigned we add it here
            if not score_sum:
                score_sum = list(diag_clean.values())[dist-1]
                successes[curr_gene][score_sum] = list(diag_clean.items())[dist_dum:dist]
            
            #And this part is necessary in case the last 2 elements had a gap <= gap_threshold
            if score_sum:
                successes[curr_gene][score_sum] = list(diag_clean.items())[dist_dum:dist]
                
            #And finally here we keep all of the scores that exceed 50% len(query).
            final_pass = [scsum for scsum in list(successes[curr_gene].keys()) if float(scsum) >= float(len(K)*final_score)]
            if final_pass:
                for ele in final_pass:
                    clean_succ[curr_gene][ele] = successes[curr_gene][ele]

    
    #print({s1:successes[s1][s3] for s1,s2 in successes.items() for s3,s4 in s2.items() if float(s3) >= float(len(K)*final_score)})

    return {s1:s2 for s1,s2 in clean_succ.items() if s2}, successes


pattern = "GTATATAACCTAAAAAGAAAATTTGATTAACGAATTCAGAACCCACACGCTGAGACAAGGGGCCCCTGATTATACACCTCTGGAAATGGAAATGGTTTGAATGGGCAAAGCCAGGAAGATATTGAGATTTTAACATATATATTATGTGACATTTGCCCTTCTTTTAGATGGTAGATAATGCTTATTTATGTAGCACGGGTCATTCATCATGGATTTGCGGAGTGGAAGACCAAATGATCATTAACCTGCTACCTGCTTTCAGAATCCTAATGATTTCAAAATGGATAACTTTACCTGTCTCTCAAATTCAGGTCTATTATCTCTCCACAAGATGCAAGCATCAATGTTGGCACCACTTTCGATATTGGGCTCACTCAACATGCTCATAACACTTAATAGAATTTTTTCTACACTTTGCACTGGCGACCATCTTTCTTCCGCTAATTCGTACATGTTAGGATCATCACCAGGGGAGTGTAGAATGGATATGCACACTTCCCCATTTGGATAAATATTTGGATGTAGTATGCTGGGTGTGAAAGTAAGTTTAGGTGGAGATAACGGATAGTCTTTAGGAAACTCTAGCTTAGCATTAAAAACACCATCAGCGTATGGCGTATCTGGAGGCCCTTGAATTAGGCAGTCCCAAATGAATATGTTATTCTCCGATTTGGGACCAGCCACTATACCAGGTGGAGAATCTTTAATTAACTGTTGAAGCTCCTTGAGGAGACGTTTCTGAGCGGTTTTCGACATGCTATGCCCTTCCAAATTACACTATTACTAGGGAAGTTCCTTTTAGGATAATCTCCTTCGTACGCTAAACGCCAGAGTTTACTTTGGCGCTTTTCGAGCTCTTGGTTTAGAGGCTGTAATCTCGTTTTCGGGTAATGGCGAAAAGGAGTTATAAGAAGGATCTCGAGACAATAAGCTGCTGCATCTTCGTGAGGTAGATGCGATGAGGCGCCTTGTTTTAAACTTGAAACAGTCTGGTAAGTTCTTCAAGCTCTATCAGAGATGATAATATTTAATGGGAACAAATATGCGTGTGCATCGTGCATCAGAGGGCATCGCTCTTCAACATGCAGGCATTCTCACCGAATGCTAATTAAGGTGAGAACTAGAGGAAGACCTTCACCCATGGCTATCAGAGACGCTATTTTAGTAGACTCTACATCGCAAAGTACAGAATACGAAAATGGTGCACAAATTGAAGGTGACTGCTGCAGCGCAATGAATCAACAGCCAATACTGTTTGTACGTGCCTCTGCTGTTAGAAAGGCAAGAATGATAAACGGAAAATTGCATATATTAATGGAGGAAGGTTTCACTGCTCATGAGCCTAAAGATATTAGCACATTTACCGATGATGGTAACAAATATATCACCGAAACGGAGTTTCTTAGGAAACACTCTCCCAAAGCTCCGGCAACAGGAACAATATCTCCGGACTCTACCAAGTCATCTTCTTCAAGTGAGAAGAAAGAGCGAAGCCGGCTTCAACAGGAGCCTATACGACATTTTTCAAACTGCTGTAAGAAAGATAAATCACAGAACCCAGCTTCTAATGGCAAGACGAACAAGGCACCGTCTGATGACATATTTACGCCATACGGCTCCTTGGAATCTACGTCCGCTTTTAACGATATTTTACAAGAAAACTACAATAGTTCTGTTCCTGGTGCGCATGACAGTTCAGAAACACTCACCCCACAAAGTACAACAACGATTGCTGCTCCTCATTCAAGCGACGTTGCTTCGAAAGTTGAAGTCCTGACTCATAAGGGCATTTTTTTAAGCACGCAGTGCTCTTGTGAAGATGAAAGCTGCCCATGTGTTAATTGTCTAATCCATAGAAGCGAAGAGGAACTGAATTCTTATATTCAACAAAGTGGTGTTCGGCCTTAATTAACTTAAGCAACATCCAATCTAAGATGAAGGATTATCTCCGGAGGATTGTAAATGCCCTGACAAGGACCGCATATGCCTCCTGGGGATAACCTTACTTGTGATGGAT"
seq = open('cerevisiae.fa')

start_time = time.time()

clean_s, all_s = fasta(seq, pattern, kmer = 5, lthr = 20, gthr = 3, final_score = 0.5) #50%
print(clean_s, all_s)

print("--- %s seconds ---" % (time.time() - start_time))