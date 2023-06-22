# %%
#The whole script takes ~300seconds to run against all the S.cerevisiae genes
import time

# Index the kmers of a sequence
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

    
def fasta(L, K, kmer, lthr, gthr, final_score):
    K_kmer = indexKmers(K,kmer) #We start by getting the index of each kmer in our smaller sequence
    successes = {} #This will contain matches that pass gap + score threshold
    clean_succ = {} #This will contain matches that pass all the thresholds
    gene_seqs = {} #We will append query sequences here

    
    for i in L: #Read each line of the file
        
        if ">" in str(i): #Skip lines without sequences
            curr_gene = str(i.rstrip('\n'))
            continue
        
        diag = {} #Diagonals dictionary
        start_points_dummy = {} #Starting points dictionary
        
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
                    start_points_dummy[d] = pos #Here we start by getting the first seen position
                    count += 1 
                    continue
                                
                if pos < start_points_dummy[d]: #Now we want to check where the identical sequence starts to match from and so we want the minimum position
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
            
            #This part is necessary in case the last 2 elements had a gap <= gap_threshold
            if score_sum:
                successes[curr_gene][score_sum] = list(diag_clean.items())[dist_dum:dist]
                
            #Here we keep all of the scores that exceed 50% len(query).
            final_pass = [scsum for scsum in list(successes[curr_gene].keys()) if float(scsum) >= float(len(K)*final_score)]
            
            
            #And now some more madness. The idea here is to try and get only the starting points that also passed our third parameter,
            #which is beeing identical to more than half of the query sequences length.
            #We will use Start to identify each starting position and we can use it later for indexing.
            #Since more than 1 points is possible, we should have as values of Start a list with elements
            #equal to the number of diagonals we found
            if final_pass:
                clean_succ[curr_gene]['Start'] = []
                
                for ele in final_pass:
                    clean_succ[curr_gene][ele] = successes[curr_gene][ele]
                    min_point = [i1 for i1,i2 in successes[curr_gene][ele]]
                    
                    if not clean_succ[curr_gene]['Start']:
                        #min_point = [2,3,4] #This is a test to check that it works for more than 1 point
                        for pt in min_point:        
                            clean_succ[curr_gene]['Start'].append(start_points_dummy[pt])
                
                #Closer to the end of this madness. This is an extra addition to get the actual sequence from the query.
                #The idea here is to manipulate our diagonals that passed all of the thresholds with the help of all the previously used dictionaries.
                
                for k,v in clean_succ.items(): #We start by creating a dictionary with keys beeing the names of genes.
                        counter = 0
                        gene_seqs[k] = []
        
                        for itm in list(v.items()): #Since the values of the previous dictionary is a dictionary in itself
                            
                            if type(itm[0]) == str: #We use this to skip the "Start"
                                continue
                            
                            for stp in v[itm[0]]: #This will return a tuple with index 0 beeing diagonal index, and index 1 beeing the diagonal score
                                        
                                total_l = v['Start'][counter] + stp[1] #So now we will use counter to keep track of indexes and get the total length
                                
                                #I noticed that while using different k-mers, besides having an exponential increase in time,
                                #the score changed, so to account for that, i use the loop seen here.
                                #The idea is that we try to find the exact position of the sub-sequence that matches,
                                #using the starting point and the total length.
                                
                                while True: 
                                    
                                    dummy_seq = K[v['Start'][counter]:total_l]
                                    
                                    #I should probably have used the position of the gene sequence to be completely sure,
                                    #but i already was lost in the script so i kept it simple, moving based on the consecutive sequence found (1685)
                                    
                                    #Aside from that, we check if the sub-sequence is contained within the gene,
                                    #and if so, we add 1 to the total length, until it is no longer contained, where we then keep the last sub-sequence used.
                                     
                                    if dummy_seq in i.rstrip('\n'):
                                        total_l += 1
                                        if K[v['Start'][counter]:total_l] in i.rstrip('\n'):
                                            continue
                                        else:
                                            gene_seqs[k].append(dummy_seq)
                                            break
                                    
                                    #The idea is same here, but instead it take account for sub-sequences longer.
                                    if dummy_seq not in i.rstrip('\n'):
                                        total_l -= 1
                                        if K[v['Start'][counter]:total_l] not in i.rstrip('\n'):
                                            continue
                                        else:
                                            gene_seqs[k].append(K[v['Start'][counter]:total_l])
                                            break
                                counter += 1  #used in case we have more than 1 consecutive diagonals.
                                

    
    #print({s1:successes[s1][s3] for s1,s2 in successes.items() for s3,s4 in s2.items() if float(s3) >= float(len(K)*final_score)})
    return {s1:s2 for s1,s2 in clean_succ.items() if s2}, successes, {s3:s4 for s3,s4 in gene_seqs.items() if s4}


pattern = "GTATATAACCTAAAAAGAAAATTTGATTAACGAATTCAGAACCCACACGCTGAGACAAGGGGCCCCTGATTATACACCTCTGGAAATGGAAATGGTTTGAATGGGCAAAGCCAGGAAGATATTGAGATTTTAACATATATATTATGTGACATTTGCCCTTCTTTTAGATGGTAGATAATGCTTATTTATGTAGCACGGGTCATTCATCATGGATTTGCGGAGTGGAAGACCAAATGATCATTAACCTGCTACCTGCTTTCAGAATCCTAATGATTTCAAAATGGATAACTTTACCTGTCTCTCAAATTCAGGTCTATTATCTCTCCACAAGATGCAAGCATCAATGTTGGCACCACTTTCGATATTGGGCTCACTCAACATGCTCATAACACTTAATAGAATTTTTTCTACACTTTGCACTGGCGACCATCTTTCTTCCGCTAATTCGTACATGTTAGGATCATCACCAGGGGAGTGTAGAATGGATATGCACACTTCCCCATTTGGATAAATATTTGGATGTAGTATGCTGGGTGTGAAAGTAAGTTTAGGTGGAGATAACGGATAGTCTTTAGGAAACTCTAGCTTAGCATTAAAAACACCATCAGCGTATGGCGTATCTGGAGGCCCTTGAATTAGGCAGTCCCAAATGAATATGTTATTCTCCGATTTGGGACCAGCCACTATACCAGGTGGAGAATCTTTAATTAACTGTTGAAGCTCCTTGAGGAGACGTTTCTGAGCGGTTTTCGACATGCTATGCCCTTCCAAATTACACTATTACTAGGGAAGTTCCTTTTAGGATAATCTCCTTCGTACGCTAAACGCCAGAGTTTACTTTGGCGCTTTTCGAGCTCTTGGTTTAGAGGCTGTAATCTCGTTTTCGGGTAATGGCGAAAAGGAGTTATAAGAAGGATCTCGAGACAATAAGCTGCTGCATCTTCGTGAGGTAGATGCGATGAGGCGCCTTGTTTTAAACTTGAAACAGTCTGGTAAGTTCTTCAAGCTCTATCAGAGATGATAATATTTAATGGGAACAAATATGCGTGTGCATCGTGCATCAGAGGGCATCGCTCTTCAACATGCAGGCATTCTCACCGAATGCTAATTAAGGTGAGAACTAGAGGAAGACCTTCACCCATGGCTATCAGAGACGCTATTTTAGTAGACTCTACATCGCAAAGTACAGAATACGAAAATGGTGCACAAATTGAAGGTGACTGCTGCAGCGCAATGAATCAACAGCCAATACTGTTTGTACGTGCCTCTGCTGTTAGAAAGGCAAGAATGATAAACGGAAAATTGCATATATTAATGGAGGAAGGTTTCACTGCTCATGAGCCTAAAGATATTAGCACATTTACCGATGATGGTAACAAATATATCACCGAAACGGAGTTTCTTAGGAAACACTCTCCCAAAGCTCCGGCAACAGGAACAATATCTCCGGACTCTACCAAGTCATCTTCTTCAAGTGAGAAGAAAGAGCGAAGCCGGCTTCAACAGGAGCCTATACGACATTTTTCAAACTGCTGTAAGAAAGATAAATCACAGAACCCAGCTTCTAATGGCAAGACGAACAAGGCACCGTCTGATGACATATTTACGCCATACGGCTCCTTGGAATCTACGTCCGCTTTTAACGATATTTTACAAGAAAACTACAATAGTTCTGTTCCTGGTGCGCATGACAGTTCAGAAACACTCACCCCACAAAGTACAACAACGATTGCTGCTCCTCATTCAAGCGACGTTGCTTCGAAAGTTGAAGTCCTGACTCATAAGGGCATTTTTTTAAGCACGCAGTGCTCTTGTGAAGATGAAAGCTGCCCATGTGTTAATTGTCTAATCCATAGAAGCGAAGAGGAACTGAATTCTTATATTCAACAAAGTGGTGTTCGGCCTTAATTAACTTAAGCAACATCCAATCTAAGATGAAGGATTATCTCCGGAGGATTGTAAATGCCCTGACAAGGACCGCATATGCCTCCTGGGGATAACCTTACTTGTGATGGAT"
seq = open('cerevisiae.fa')

start_time = time.time()

clean_s, all_s, seq_dict = fasta(seq, pattern, kmer = 5, lthr = 20, gthr = 3, final_score = 0.5) #50%

#Clean s is a dictionary with genes as keys. Each gene has a dictionary with different keys.
#'Start' Shows the starting position of each diagonal within the query sequence. The rest of the keys are 
#the sums of the consecutive (or not) identical matches, and their values are the (diagonal_index, diagonal_score)
print(clean_s, end = ' genes that pass all threshold and are identical to more than half of query sequence length \n')

#Same logic as above only here we do not include the starting positions.
print(all_s, end = ' genes that passed the length and gap threshold but not query sequence length threshold \n')

#This contains the query sequences that passed all the thresholds and were identical.
#Again it is a dictionary with genes as keys and the values is a list of the sequences. Since there was
#only 1 diagonal that passed all the thresholds, i cannot confirm that it works properly, but the idea is
#that we will have a list of sequences in the order of the starting points and diagonal indexes as seen above
print(seq_dict, end = ' \nAbove we can see the sequences found to be identical in the query sequence \n')

print("--- %s seconds ---" % (time.time() - start_time))

print(f'The highest matching sequence is for the gene YMR021C and can be seen here with a length of {len(seq_dict[">YMR021C"][0])}: \n{seq_dict[">YMR021C"][0]} .')
