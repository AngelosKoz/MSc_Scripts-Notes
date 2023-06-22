sequence = 'ACACAGTACACGTATACCCAGTTTGCACAGTTTT'
pattern = 'AGTT'
matches = 0
for i in range(len(sequence)-3): # consider that any string has n-k+1 substrings of length k
    string = sequence[i:i+4]
    if (string == pattern):
        matches += 1
        print("match found at position", i)
print("Pattern was found ", matches, " times.") 


#----------------------------#


# Reading ecoli genome

file = open('files/ecoli.fa', 'r')
ecoli = ''
count = 0

for line in file:
    count += 1
    if (count > 1): # the first line contains the non-sequence header so we discard it 
        ecoli += line.replace("\n", "") # we string the newline character from the end of each line

# Using time to measure time of execution
import time
start_time = time.time()

# Pattern search

pattern = 'AGTTAGGCCT'
#pattern = 'TCGGCATCAG'
matches = 0
k = 10

for i in range(len(ecoli)- k + 1): # consider that any string has n-k+1 substrings of length k
    string = ecoli[i:i+k]
    if (string == pattern):
        matches += 1
print("Pattern was found ", matches, " times.")

print("--- %s seconds ---" % (time.time() - start_time))
