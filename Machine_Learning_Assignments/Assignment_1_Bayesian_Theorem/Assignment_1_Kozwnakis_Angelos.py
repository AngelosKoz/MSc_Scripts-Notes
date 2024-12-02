import sys
print(f"Python version: {sys.version}\n")

#Created using pyenv virtualenv 3.6.15

# Exercise 2

##### a) #####
# P(photon passes Earth atmosphere) = P(Photon) = 1e-7
# Since P(Photon) + P(Not Photon) = 1 <--> P(Not Photon) = 1 - 1e-7
# Detector FPR : 10% --> P(Detection | Not Photon)
# Detector TPR : 85% --> P(Detection | Photon )

# P(Photon | Detection) = ? (1)
# We can use Bayes theorem on (1) so it becomes:
# (1) = P(Detection | Photon)*P(Photon) / P(Detection)
# We know P(Detection | Photon) = 0.85, P(Photon) = 1*(10**-7)

# We need to find P(Detection). 
# Let Photon + Not_Photon be the only two possible outcomes (We assume mutually exclusive and exhaustive events)
# Then from the Law of Total Probability we can say that:
# P(Detection) = P(Detection | Photon)*P(Photon) + P(Detection | Not Photon)*P(Not Photon)
# We know all the above terms so we just need to substitute now.

P_Photon = 1e-7
P_Not_Photon = 1 - P_Photon
P_Detection_Photon = 0.85
P_Detection_NotPhoton = 0.1

P_Detection = (P_Detection_Photon*P_Photon) + (P_Detection_NotPhoton*P_Not_Photon)
print(f'Probability of detection = {P_Detection}')
P_Photon_Detection = (P_Detection_Photon*P_Photon)/P_Detection
print(f'The probability that a photon package is actually received when a detector reports detection is: \n{P_Photon_Detection}\nor simplified\n{P_Photon_Detection:.10f}\n')


##### b) #####
# Photon package: 100 photons
# Photon energy (PE):  10, 20, 30 or 40. 
# P(PE_10) = P(PE_20) = P(PE_30) = P(PE_40) = 25% = 0.25

import random
import matplotlib.pyplot as plt
from collections import Counter

PE = [10,20,30,40]
P_energy = [0.25, 0.25, 0.25, 0.25]# P( PE_10 = P(PE_20) = P(PE_30) = P(PE_40) )
photons = 100

#pack_energy = random.choices(PE, k=photons) #Here we assume they all have the same probability to be chosen.
pack_energy = random.choices(PE, P_energy, k=photons) #By adding the P_energy list here, we could give different probabilities for each energy.
#Since all energy levels are equiprobable, the two snippets above produce same results in our case
print(f'Amount of photons for each energy level {dict(Counter(pack_energy))}\n\t\t\t\t\t\t(Energy: Photons)')


#For 1 photon package
plt.hist(pack_energy, bins = 60, color='blue', edgecolor = 'black')
plt.xlabel('Electronvolts')
plt.ylabel('Number of Photons (Total: 100)')
plt.title('Single Photon Package Energy Distribution')
plt.show()

#For different number of total packages
package_list = [2_000, 20_000, 200_000, 1_000_000]

#Each loop will add the sum of the package (after energies have been distributed) to the list
#We will run it for different number of total packages to see how it differentiates.

for packages in package_list:
    total_energy = [] #Initialize an empty list
    for _ in range(packages):
        pack_energy = random.choices(PE, P_energy, k=photons)
        total_energy.append(sum(pack_energy))
        
    plt.hist(total_energy, bins = 60, density = True, color='blue', edgecolor = 'black')
    plt.xlabel('Package Total Energy (Electronvolts)')
    plt.ylabel(f'Relative Frequency of Packages (Total Packages: {packages})')
    plt.title('Probability Density of Photon Package Total Energy')
    plt.show()
    

##### c) #####
#Mean = 1e-7 and standard deviation = 9e-8
import numpy as np
import matplotlib.pyplot as plt

mu = 1e-7  # Mean of the normal distribution
sigma = 9e-8  # Standard deviation of the normal distribution
P_Detection_Photon = 0.85  # True Positive Rate (TPR)
P_Detection_NotPhoton = 0.1 # False Positive Rate (FPR) 

# Randomly sampling from N(mu=1e-7, sigma=9e-8) for different sample sizes
sample_size_list = [10_000, 100_000, 1_000_000]


# Two aproaches follow. One (c_a) is first sampling our sample size and then filtering out the negatives.
# This aproach also plots with the negative values included to check how the distribution would look like if probabilities could be negative.
# The second aproach (c_b) will sample as many values as the sample size, while dropping the negative ones.

##### c_a) #####

for sample_n in sample_size_list:
    norm_samples_unfilter = np.random.normal(mu, sigma, sample_n)#First we sample for our parameters and sample size
    norm_samples = [i for i in norm_samples_unfilter if i >= 0] #Then we filter out the negative probabilities
    print(f'For n={sample_n}: Removed {sample_n - len(norm_samples)} negative values\n')
    
    
    P_Photon_Detection = np.zeros(shape = len(norm_samples), dtype = float) #Initialize empty array to fill in with P(Photon | Detection) values
    index = 0
    for P_Photon in norm_samples:
        P_Not_Photon = 1 - P_Photon
        P_Detection = (P_Detection_Photon*P_Photon) + (P_Detection_NotPhoton*P_Not_Photon)
        P_Photon_Detection[index] = (P_Detection_Photon*P_Photon)/P_Detection
        index += 1

    plt.hist(P_Photon_Detection, bins=60, color='blue', edgecolor = 'black')
    plt.title('Distribution of P(Photon|Detection) -- Filtered')
    plt.xlabel(f'P(Photon|Detection) (Sample size: {len(norm_samples)})')
    plt.ylabel('Probabilitiy Counts')
    plt.show()
    
    #This snippet will do the same, this time for the unfiltered (we keep negative probabilities for the sake of visualization)
    P_Photon_Detection_unfilt = np.zeros(shape = len(norm_samples_unfilter), dtype = float)
    index = 0
    for P_Photon in norm_samples_unfilter:
        P_Not_Photon = 1 - P_Photon
        P_Detection = (P_Detection_Photon*P_Photon) + (P_Detection_NotPhoton*P_Not_Photon)
        P_Photon_Detection_unfilt[index] = (P_Detection_Photon*P_Photon)/P_Detection
        index += 1

    plt.hist(P_Photon_Detection_unfilt, bins=80, color='blue', edgecolor = 'black')
    plt.title('Distribution of P(Photon|Detection) (including negatives)')
    plt.xlabel(f'P(Photon|Detection) (Sample size: {sample_n})')
    plt.ylabel('Probability Counts')
    plt.show()

##### c_b #####

for sample_n in sample_size_list:
    pos_norm_samples = []
    #The difference here is we use a while loop to loop through the random normal distribution until we get the specified sample size while dropping negative probabilities.
    while len(pos_norm_samples) != sample_n:    
        norm_samples_no_filtering = np.random.normal(mu, sigma, sample_n-len(pos_norm_samples))
        norm_samples_no_filtering = [i for i in norm_samples_no_filtering if i >= 0]
        pos_norm_samples.extend(norm_samples_no_filtering)
    
    P_Photon_Detection = np.zeros(shape = sample_n, dtype = float)
    index = 0    
    for P_Photon in pos_norm_samples:
        P_Not_Photon = 1 - P_Photon
        P_Detection = (P_Detection_Photon*P_Photon) + (P_Detection_NotPhoton*P_Not_Photon)
        P_Photon_Detection[index] = (P_Detection_Photon*P_Photon)/P_Detection
        index += 1
    
    plt.hist(P_Photon_Detection, bins=60, color='blue', edgecolor = 'black')
    plt.title('Distribution of P(Photon|Detection) ')
    plt.xlabel(f'P(Photon|Detection) (Sample size: {sample_n})')
    plt.ylabel('Probabilitiy Counts')
    plt.show()

