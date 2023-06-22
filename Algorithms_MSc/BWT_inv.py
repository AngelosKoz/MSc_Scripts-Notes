# %% Part 1 of the exercise
def BWT(seq):
    counter = 0
    bwt_ext = {}
    
    
    
    #This will account for sequences that implement their own terminals and/or suffixes
    if seq.isalpha():
        seq += '$'
    
    #A simple for loop to get all the cyclic rotations. 
    #The idea is we just displace the last index at the beginning.
    for rotation in range(len(seq)):
        seq = seq[-1] + seq[:-1] #And this is the displacement
        bwt_ext[counter] = seq #We use a dummy counter to add all the possible strings
        counter += 1

    bwt_ext = sorted(bwt_ext.items(), key = lambda x: x[1]) #Here we sort based on the string
    
    return [s[-1] for c,s in bwt_ext] #This will return only the last letter of the string

#print(BWT('^BANANA|'))
#S = 'GAGTAAGTCA'
#print(BWT(S))


# %% Part 2 of the exercise

#The function will give use the original sequence from a CRT

def BWTinverse(bwt):
    
    #We use these transformations to keep track of the original positions of the CRT
    #We will use dictionaries with keys beeing positions and values the string
    bwt_enu = dict(sorted(dict(enumerate(bwt)).items(), key = lambda x: x[1]))
    #Here we sort the CRT while keeping tracking of the sorted positions to get the first column
    inv_sort = dict(enumerate(sorted(bwt)))
    
    got_seq = '' #Initialize the original sequence
    counter = 0
    increment = '1' #This will be used to add endings to the letters in the order appeared
    
    #Now that we have both the BWT and the first column enumerated, we need to assign values at the end of the letters
    #In the loop i chose to go through the BWT to use later as indexing
    for enu,let in bwt_enu.items():

        if not let.isalpha(): #Since it is sorted the first element will always be '$' (or other indicators although we not implement them here)
            counter += 1
            continue
        
        if counter == len(bwt)-1: #We take into account the last letter since our loop uses a consecutive letter check (index == index +1)
            
            if inv_sort[counter] == inv_sort[counter-1][0]: #If letter is the same we use the increment to increase it
                
                inv_sort[counter] += increment
                bwt_enu[enu] = inv_sort[counter] #At the same time we use the counter as index of the 'first column'.
                #Since the bwt is first enumerated then sorted we assure the letter assignment is proper
                break
            
            #This accounts for the last digit beeing different (maybe for words rather than sequences)
            inv_sort[counter] += '1'
            bwt_enu[enu] += '1'
            
            break
        
        if inv_sort[counter][0] == inv_sort[counter+1][0]: #We check if consecutive letters are the same
            #And we apply the same logic as for the last letter    
            inv_sort[counter] += increment
            bwt_enu[enu] = inv_sort[counter]
            #Since it is not the last letter we need to keep track of the occurences.
            increment = str(int(increment)+1) #So we increase them by 1
            counter += 1
            continue 
        
        #If the current letter does not match the next, we assign the increment
        inv_sort[counter] += increment
        bwt_enu[enu] = inv_sort[counter]
        #And then we set the increment back to 1 to start counting the next letters
        increment = '1'
        counter += 1
    
    #Since the first element of the first column is '$' we set it to a dummy variable to begin the sequence building
    dum = list(bwt_enu.keys())[0] 

    while len(got_seq) < len(bwt_enu)-1: #Since '$' is included in the sequence, we want 1 less
        import re
        #The logic follows the steps at increment assignment.
        #Since the properties of BWT - First column is opposite letter correspondance, we use the BWT key (index)
        got_seq += inv_sort[dum] #as key to the first column and assign the value to the string, 
        dum = inv_sort[dum] #then set the dummy variable to a temporary index value
        got_seq = got_seq[:-len(dum[1:])] #Remove the last elements of the string which will always be an integer    
        dum = list(bwt_enu.keys())[list(bwt_enu.values()).index(dum)] #And now we use the dummy index to find the corresponding position in the BWT
        
    return got_seq

#I did not implement a solution to terminals besides '$'.


# %%
# Accession name: A2VEY9
import time
start_time = time.time()


file = open('protein.fa', 'r')
protein = ''

for line in file:
    if '>' in line:
        continue
        
    protein += line.replace("\n", "")        

bwt = BWT(protein)
seq = BWTinverse(bwt)

print(protein == seq) #Checks out, it is the same sequence


print("--- %s seconds ---" % (time.time() - start_time))



# %% A test with a large protein to check time. ~1 minute for titin. (38,126 length)
#Propably the inverse can be optimised.
#Accession name: A0A6P7YNV3
import time
start_time = time.time()
file = open('titin.fa', 'r')
titin = ''



for line in file:
    if '>' in line:
        continue
    
    titin += line.replace("\n", "")

bwt_titin = BWT(titin)

print(len(bwt_titin))
print("--- %s seconds ---" % (time.time() - start_time))

seq_titin = BWTinverse(bwt_titin)

print(titin in seq_titin) #Checks out it is the same sequence.

print("--- %s seconds ---" % (time.time() - start_time))