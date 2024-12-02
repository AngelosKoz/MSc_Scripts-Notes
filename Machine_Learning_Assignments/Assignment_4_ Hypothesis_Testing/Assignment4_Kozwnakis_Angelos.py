import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2

# In the tables we have the joint probability distribution of two random variables X and Y.
# What this shows us, is the probability of having specific values taken for the variables. 
# Since X and Y take only 0 and 1 as values, the table contains the probability P(X = x, Y = y).
# The tables can be interpreted as follows. We will use table 1 as a reference.
# There is a 5% probability that both X=0 and Y=0 at the same time --> P( X = 0, Y = 0) = 0.05
# P(X = 0, Y = 1) = 0.1 (10% that both X=0 and Y=1)
# P(X = 1, Y = 0) = 0.25 (25% that both X=1 and Y=0)
# P(X = 1, Y = 1) = 0.6 (60% that both X=1, Y=1)
# Since the above represent all possible outcomes, adding those up will give us 1.

# To check for independence of a random variable, we will use a statistical test (in our case chi-squared test).
# For that, we need a null hypothesis (H0) and an alternative hypothesis (H1).
# The H0 in our case is that X is independent of Y, for alpha = 0.05.
# X and Y are independent if the joint probability equals the product of the individual probabilities for all values of x and y.
# The above can be written as H0: P(X = x, Y = y) = P(X = x)*P(Y = y) (for all x,y values)
# Our alternative hypothesis H1 is that there is some dependence between the variables. This means that the outcome of one variable
# influences the outcome of the other and happens when the product of the individual probabilities does not equal the joint probability.
# The above can be written as H1: P(X = x, Y = y) not equal to P(X = x)*P(Y = y) (for some x,y values). 
# Note that for dependence we dont need the condition to be true for all x,y values.

# First we pass the tables as arrays.
jp_distr_T1 = np.array([[0.05, 0.1], [0.25, 0.6]])
jp_distr_T2 = np.array([[0.08, 0.23], [0.17, 0.52]])
jp_distr_T3 = np.array([[0.1, 0.2], [0.2, 0.5]])
jp_distr_T4 = np.array([[0.1, 0.1], [0.1, 0.7]])
jp_distr_T5 = np.array([[0.15, 0.1], [0.1, 0.65]])
jp_distr_T6 = np.array([[0.45, 0.05], [0.05, 0.45]])
all_jp_T = [jp_distr_T1, jp_distr_T2, jp_distr_T3, jp_distr_T4, jp_distr_T5, jp_distr_T6]
sample_list = [25, 100, 500]

##### A) #####
# To sample from a joint probability distribution. We use the distribution of probabilities and the sample size:
def sample_jd(distr, sample):
    # We use flatten to get all the probabilities in a single array (N,)
    flat_distr = distr.flatten()
    
    # Since we only have two random variables that can take the values 0 and 1, the possible outcomes are:
    # (0,0), (0,1), (1,0), (1,1). 
    # We can create the above possible values like so:
    xy_values = [(i, j) for i in range(distr.shape[0]) for j in range(distr.shape[1])]
    
    #Now we randomly sample from probability array. This will create a list that correspond to the 4 possible x and y values
    samples = np.random.choice(range(len(flat_distr)), size=sample, p=flat_distr) # We change p to alter the probabilities
    # In the list above, 0 corresponds to (0,0), 1 : (0,1), 2: (1,0), 3: (1,1), with probabilities as seen in the distribution 
    
    # Now we find based on the index the corresponding outcomes (X,Y values)
    return [xy_values[i] for i in samples]


##### B.i) #####
# Function to perform the chi-squared test for independence
def chi_square_ind(distr, n):
    # Sample from joint distribution
    samples = sample_jd(distr, n)
    
    # So the formula for chi-square is :
    # T = Σ_{i,j} [ (Ο_{i,j} - Ε_{i,j})^2 / Ε_{i,j} ] (i,j here is the values of X and Y)
    # O_{i,j}: Observed frequencies for x and y values (as generated from our sampling)
    # E_{i,j}: Expected frequencies for x and y values (Considering H0 is true, so X and Y independent)
    
    # Calculate observed (O_{i,j}) frequencies
    obs = np.zeros((2, 2), dtype=int)
    for x, y in samples: # Access the tuple values (X , Y) {0,1}
        obs[x, y] += 1 # Create the observe frequency table (total sum will be equal to n)

    # Calculate expected (E_{i,j}) frequencies (assuming H0 is true -- independence)
    row_sum = obs.sum(axis=1) # Get row sums
    col_sum = obs.sum(axis=0) # Get column sums
    # Multiply the row sum and column sum, then divide by all sum (n)
    exp = np.outer(row_sum, col_sum) / n
    # To avoid division by zero errors, we can add a very small value to not alter the outcomes
    # Compute chi-squared statistic
    T = ((obs - (exp + 1e-10)) ** 2 / (exp + 1e-10)).sum()
    
    # DF for 2x2 contigency table is (num_rows - 1) * (num_cols - 1)
    df = (len(distr) - 1) * (len(distr[0]) - 1)
    
    # For the p-value, we have : p_val = 1 − P(T ≤ t_obs )
    # Calculate the p-value using the cumulative distribution function (CDF)
    p_val = 1 - chi2.cdf(T, df)
    
    return T, p_val, obs, exp, samples


# We can use a dictionary to append all of the values generated in our function and access them in an easy manner,
sample_dict = {}
alpha = 0.05
for n in sample_list:
    sample_dict[f'n={n}'] = {} # Initiate the first layer containing sample size
    
    for j in range(len(all_jp_T)):
        distr = all_jp_T[j]
        T, p_val, obs, exp, samples = chi_square_ind(distr, n)
        
        # The second layer is the table, containing within all the outcomes of our independence test function
        sample_dict[f'n={n}'][f'Table_{j+1}'] = {
            'T_stat': T, # Chi-square T statistic
            'p_val': p_val, # Independence test p-value
            'obs': obs.tolist(), # Observed frequencies (.tolist() makes it more readable by remove the "array()"" in the output)
            'exp': exp.tolist(), # Expected frequencies
            #'samples': samples  # The samples (as asked in part A of exercise). Skipped for memory usage/readability 
        }
        
        # We (probably) reject the null hypothesis when our p-value is lower than the threshold alpha.
        # In our case alpha = 0.05, so we probably reject the null hypothesis when p-val < 0.05 meaning that we probably accept
        # the alternative hypothesis, that X and Y are not independent.
        if p_val < alpha:
            print(f'We probably reject null hypothesis for sample size n={n} and Table {j+1} since p-value = {p_val}')

        # From the results, we notice a consistency when sample size increases for Tables : 4,5 and 6. This is true in some cases
        # for smaller sample sizes, but not always.


##### B.ii) #####
# Steps for permutation: 
# Randomly shuffle the observed frequencies (In case of potential dependency between variables)
# Recalculate chi-squared statistic (T) for each permutation
# Check if the new statistic is higher than the previously calculated ("real")
# Repeat as many times is the permutation loops (in our case 200)
# Calculate the pvalue as : (Number of times perm statistic is higher than "real" + 1) / (Total permutations +1)
# The reason we used Permutations - 1 as initial loop counter is to make the division more efficient (200 instead of 201)

def chi2_perm_test(n, obs, exp= None, T=None, perms=199):
    
    # Same logic as the chi square independent test function. We can skip it if we provide the expected frequencies and statistic
    if exp is None and T is None:
        row_sum = obs.sum(axis=1)
        col_sum = obs.sum(axis=0)
        exp = np.outer(row_sum, col_sum) / n
        T = ((obs - (exp + 1e-10)) ** 2 / (exp + 1e-10)).sum()

    counter_stat = 0 # Initial counter of how many times the permutation statistic is higher than the previously calculated

    for _ in range(perms):
        # Shuffle the observed frequencies within each row.
        shuf_obs = np.array([np.random.permutation(row) for row in obs])

        # Calculate chi-squared statistic using shuffled observed frequencies
        perm_T = ((shuf_obs - (exp+1e-10)) ** 2 / (exp+1e-10)).sum()

        # If permuted statistic is higher than the previously calculated, add 1
        if perm_T >= T:
            counter_stat += 1
            
    # Calculate the p-value
    p_val = (counter_stat + 1) / (perms + 1) # We use permutations - 1 as initial loop, to make the division better (thus for 200, we use 199)
    return p_val

# Again we use a dictionary for the results, this time we can only keep the p-values to compare non-perm p-value vs perm p-value
sample_dict_perm = {}
alpha = 0.05
sampling = 100

for n in sample_list:
    sample_dict_perm[f'n={n}'] = {}  # Initiate the first layer containing sample size
    
    for j in range(len(all_jp_T)):
        distr = all_jp_T[j]
        p_val_list = [] # List to include all the pvalues before permutation
        perm_p_val_list = [] # List to include corresponding pvalues after permutation
        
        # Loop to re-sample for "sampling" times (in our case 100)
        for _ in range(sampling):
            T, p_val, obs, exp, samples = chi_square_ind(distr, n)
            if p_val == 0:
                p_val = p_val + 1e-20
            p_val_list.append(p_val)
            
            perm_p_val = chi2_perm_test(n=n, obs=obs, exp=exp, T=T, perms=199)
            perm_p_val_list.append(perm_p_val)

        sample_dict_perm[f'n={n}'][f'Table_{j+1}'] = {
            'p_val': p_val_list,  # Independence test p-value
            'perm_p_val': perm_p_val_list  # Permutation p-value
        }

        # And we can plot at the same time both the non-permutation p-values and permutation p-values
        # 1) For chi-square 
        plt.subplot(1, 2, 1)
        plt.hist(p_val_list, bins=30, color='blue')
        plt.title(f'x² - Distribution of p-values\n(n={n}, Table_{j+1})')
        plt.xlabel('p-value')
        plt.ylabel('Frequency')
        #plt.xlim(-0.1, 1.1) # Did not include since it disrupts the bins
        #plt.ylim(0, 100) # Did not include since the results are not as readable
                
        # 2) For permutation 
        plt.subplot(1, 2, 2)
        plt.hist(perm_p_val_list, bins=30, color='green')
        plt.title(f'x² - Permutation p-values\n(n={n}, Table_{j+1})')
        plt.xlabel('p-value')
        plt.ylabel('Frequency')
        #plt.xlim(-0.1, 1.1) 
        #plt.ylim(0, 100)
        
        plt.tight_layout()
        plt.show()
        
        # I did not include xlim or ylim because it distorted the plots, since most, if not all of the values, were 0.
