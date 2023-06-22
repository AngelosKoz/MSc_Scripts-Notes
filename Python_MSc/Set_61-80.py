# %% Askhsh 61
from collections import Counter

def f61():
    cnt = Counter()

    with open('words_greek_normalized.txt') as greek_words:
        for words in greek_words:
                cnt[words[0]] +=1
        return cnt
            
            
f61()
# %% Askhsh 62
from collections import defaultdict
def  f62():
    dict_62 = defaultdict(int)
    with open('words_greek_normalized.txt') as greek_words:
        for words in greek_words:
            dict_62[words[0]] +=1
        return dict_62
f62()
# %% Askhsh 63
import itertools
def f63():
    dict_63 = {}
    with open('words_greek_normalized.txt') as greek_words:
        for words in greek_words:
            iter_63 = itertools.groupby(words)
            for key, group in iter_63:
                if not key in dict_63:
                    dict_63[key] =  1
                    break
                else:
                    dict_63[key] +=1
                    break
        return dict_63
f63()

# %% Askhsh 64
def f64():
    dict_64 = {}
    with open('words_greek_normalized.txt') as greek_words:
        for words in greek_words:
            if not words[0] in dict_64:
                dict_64[words[0]] = 1
                continue
            else:
                dict_64[words[0]] += 1
    return dict_64
f64()
# %%  Askhsh 65
import random, math
def prime_num(): #Tha mporoysa na valw san orisma to 1000 kai meta ena 1000+9 gia na parw tous prwtous alla afoy exoyme fixed arithmo thewrw pws den xreiasthke. An eixame orisma diaforetiko gia to range tou n tote tha mporoyse na ginei auth h allagh
    primes_65 = []
    prime = True
    for i in range(2,1010):
        for j in range(2, int(math.sqrt(i) + 1)):
            if (i % j) == 0:
                prime = False
                break
            else:
                prime = True
        if prime:
            primes_65.append(i)
            
    return primes_65

def f65():
    primes = set(prime_num())
    how_many_not_primes = 0
    for i in range(10_000):
        check_condition = set()
        rand_65 = random.randint(1,1_001)
        range_rand_65 = set(list(range(rand_65,(rand_65+10))))
        check_condition = primes & range_rand_65
        if not check_condition:
            how_many_not_primes +=1
    return how_many_not_primes/10_000
f65()


# %% Askhsh 66 Proud for this short code! 

import itertools

def f66(a,b,c):
    pre_ind = min(itertools.product(a,b,c), key = lambda x: (0.5)*abs( ((x[0][0] - x[2][0]) * (x[1][1] - x[0][1])) - ((x[0][0]-x[1][0]) * (x[2][1] - x[0][1])) ))    
    return a.index(pre_ind[0]), b.index(pre_ind[1]), c.index(pre_ind[2])

   
S_1 = [
 (-44.01, -37.17),
 (-40.15, 22.65),
 (35.69, -25.12),
 (-2.2, 36.35),
 (-37.96, 12.24),
 (10.65, 29.59),
 (-22.78, 48.35),
 (-42.13, 33.88),
 (-49.78, -42.15),
 (-11.16, -4.25)
]

S_2 = [
 (-15.1, -44.87),
 (5.62, -44.22),
 (-8.6, -38.89),
 (2.22, 41.26),
 (-43.73, 0.64),
 (24.69, -4.72),
 (6.94, 39.15),
 (-0.62, 24.1),
 (45.6, -11.61),
 (-27.31, -11.51),
]

S_3 = [
 (25.67, -47.55),
 (-7.55, 44.2),
 (-32.86, 0.58),
 (-39.38, 11.36),
 (-20.08, 25.73),
 (-5.55, 43.33),
 (37.67, 41.97),
 (-21.57, 26.77),
 (8.84, 12.54),
 (24.12, 3.97),
]    
print(f66(S_1,S_2,S_3))
# %% Askhsh 67
import itertools
def f67(x):
    pre_ind = min(itertools.combinations(x,3), key = lambda x: (0.5)*abs( ((x[0][0] - x[2][0]) * (x[1][1] - x[0][1])) - ((x[0][0]-x[1][0]) * (x[2][1] - x[0][1])) ))    
    return x.index(pre_ind[0]), x.index(pre_ind[1]), x.index(pre_ind[2])
    


S = [
 (-15.01, 6.31),
 (3.48, 18.51),
 (-1.93, -4.41),
 (-2.54, -36.13),
 (-10.62, 35.12),
 (-15.29, 31.85),
 (15.02, -1.13),
 (47.69, -41.41),
 (48.33, 39.15),
 (8.99, 10.17),
 (19.13, 29.09),
 (22.71, 49.55),
 (1.81, -32.67),
 (5.05, -37.53),
 (-44.16, -34.83),
 (-16.93, 8.36),
 (1.39, -15.56),
 (-43.82, 11.6),
 (-6.97, 7.53),
 (30.39, -21.45),
 (45.29, 8.12),
 (-9.02, -48.07),
 (-7.44, 43.51),
 (3.74, 25.33),
 (26.97, 12.57),
 (32.59, -43.21),
 (37.79, -2.49),
 (-8.71, 32.03),
 (-42.4, 17.28),
 (-36.22, 36.38),
 (-24.13, -0.94),
 (30.56, 33.76),
 (-4.91, 41.0),
 (5.73, -42.27),
 (42.46, -38.52),
 (-28.19, 25.26),
 (42.07, -5.42),
 (43.16, -33.38),
 (47.76, -41.65),
 (-25.64, -13.43),
 (-33.28, -49.91),
 (-21.35, -24.52),
 (24.87, -44.45),
 (-21.44, -40.46),
 (29.08, 6.97),
 (-19.77, 19.68),
 (14.17, -30.92),
 (-36.93, -28.72),
 (-41.14, -14.3),
 (21.14, -35.95),
 (-17.48, 22.8),
 (35.71, 49.83),
 (-15.76, -6.71),
 (20.43, -46.28),
 (-46.63, -34.48),
 (40.85, 16.15),
 (-17.89, 38.35),
 (-44.59, 29.8),
 (-5.86, -29.39),
 (-40.21, 43.27),
 (-30.18, 11.52),
 (33.91, 2.83),
 (39.39, -22.01),
 (-15.58, 12.68),
 (11.95, -49.36),
 (-17.09, 41.52),
 (-41.8, -14.44),
 (22.95, -10.94),
 (-16.21, -40.27),
 (43.91, 7.18),
 (-1.12, -25.06),
 (-17.88, -35.94),
 (18.94, -45.69),
 (-47.16, 20.15),
 (-26.55, 0.88),
 (-18.03, 18.6),
 (38.33, 10.81),
 (-23.4, -14.08),
 (-44.97, 6.51),
 (10.98, -21.74),
 (9.54, 20.15),
 (-37.17, 18.6),
 (9.23, 31.8),
 (-47.1, 2.04),
 (-16.45, -1.91),
 (25.31, -12.45),
 (-45.12, -5.21),
 (-47.99, -12.63),
 (36.58, -19.4),
 (34.72, 9.76),
 (-44.51, -11.83),
 (45.78, 13.84),
 (15.19, -2.97),
 (39.44, 49.41),
 (-9.23, 26.34),
 (13.55, 18.85),
 (-11.43, 14.66),
 (-32.99, 33.92),
 (-22.91, -16.26),
 (28.35, -3.41)]

print(f67(S))
# %%
# άσκηση 68 #Works with the file ask_68.py

from ask_68 import closer
print(closer(25,[20,50,60,60]))

# %% Askhsh 69 Mallon ginetai kai xwris th lista sto itertools, katalabainw pws einai bary
#Side note: Me tous generators poy eidame th deytera tha ginotan mallon alla kathws auth htan h arxikh moy lush to afhsa etsi

import itertools
def f69():
    pool_69 = list(range(1,50))
    to_remove = (3,13,23,33,43)
    counter = 0
    for i in to_remove:
        pool_69.remove(i)

    lotto = list(itertools.combinations(pool_69 , 6))
    for i in lotto:
        if any(j%2 == 0 for j in i):
            counter += 1

    return counter

print(f69())
# %% Askhsh 70
def f70():
    letters = set(['β','γ','δ','ζ','θ','κ','λ','μ','ν','ξ','π','ρ','σ','τ','φ','χ','ψ','α','ε','η','ι','ο','υ','ω'])
    dict_70 = {}
    with open('words_greek_normalized.txt') as greek_words:
        for words in greek_words:
            cnt_let = list(set(words))
            cnt_let.remove('\n')
            cnt_let = tuple(cnt_let)
            if len(cnt_let) > 4:
                continue
            else:
                if not cnt_let in dict_70:
                    dict_70[cnt_let] = 1
                    continue
                else:
                    dict_70[cnt_let] += 1
        return max(dict_70, key = lambda x : dict_70[x])
        
print(f70())

# %% Askhsh 71 #Sigoura exei ginei kapoio lathos me ta clusters kai isws sthn askhsh 74 olotela giati pairnw mono 1 gia ta clusters. 
# Tha ektimoysa an stelnate kai lyseis (opws panta) gia ayta giati fainetai poly endiaferon!
from math import sqrt
import numpy as np

def f71_noloop(ar):
    ele = 0
    ind = 0
    temp = []
    result = []
    while ele != len(ar):
        dist = sqrt((ar[ele][0] - ar[ind][0])**2 + (ar[ele][1] - ar[ind][1])**2)
        temp.append(dist)
        ind += 1
        if ind == len(ar):
            ele += 1
            ind = 0
            result.append(temp)
            temp = []
    return np.array(result)

a = [(1,2),(-1,4),(0,5)]

print(f71_noloop(a))


# %% Askhsh 72 
import numpy as np
def f72_noloop(c, minPTS, eps):
    all_distances = f71_noloop(c)
    index = 0
    results = []
    while index != len(all_distances):
        if sum(all_distances[index]<=eps)-1 > minPTS:
            results.append(index)
            index += 1
        else:
            index +=1
    return np.array(results)

C = np.array([
    [0.3, 0.5],
    [0.9, 0.2],
    [0.8, 0.2],
    [0.5, 0.3],
    [0.6, 0.5],
    [0.3, 0.4],
    [0.3, 0. ],
    [0.3, 0.2],
    [0.2, 0.1],
    [0. , 0.5]])
print(f72_noloop(C, minPTS=3, eps = 0.3)) # Τυπώνει 0, 3, 7 

# %% Askhsh 73
import numpy as np
def f73_noloop(c,minPTS,eps):
    CP = f72_noloop(c,minPTS,eps)
    non_CP = np.array(list(set(list(range(len(c)))) - set(CP)))
    return CP, non_CP


C = np.array([
    [0.3, 0.5],
    [0.9, 0.2],
    [0.8, 0.2],
    [0.5, 0.3],
    [0.6, 0.5],
    [0.3, 0.4],
    [0.3, 0. ],
    [0.3, 0.2],
    [0.2, 0.1],
    [0. , 0.5]])

CP, non_CP = f73_noloop(C,3,0.3)

print(CP)
print(non_CP)

# %% Askhsh 74 I did not know if i should append values when comparing to oneself
import numpy as np
import random


def f74_noloop(nn, cp, minPTS, eps, r):
    
    S = np.array([], dtype=int)
    r_all_val = list(range(len(nn)))
    r_val = r
    r_all_val.remove(r_val)
    ind = 0
    check = 1

    while check:

        while ind != len(cp):
            
            nn1 = nn[cp[ind]][r_val]
            
            if nn1<=eps:
                
                if cp[ind] != r_val:#we use this if we dont want the distance from themselves
                        
                    S = np.append(S,[cp[ind]])
                    ind +=1
                else:
                    
                    ind+=1
            else:
                ind+=1
        ind = 0
        condition = 1
        
        while condition:
            
            r_val = random.randint(0,len(nn)-1)
            
            if not len(r_all_val):
                
                condition = 0
                check = 0
            
            else:
            
                if r_val in r_all_val:
                    
                    r_all_val.remove(r_val)
                    condition = 0
    
    return np.array(list(set(S)))

C = np.array([
    [0.3, 0.5], #0
    [0.9, 0.2],
    [0.8, 0.2],
    [0.5, 0.3], #3
    [0.6, 0.5],
    [0.3, 0.4],
    [0.3, 0. ],
    [0.3, 0.2], #7
    [0.2, 0.1],
    [0. , 0.5]])


cp, non_cp = f73_noloop(C,3,0.3)
nn = f71_noloop(C)
r = random.randint(0,len(C)-1)

print(f74_noloop(nn,cp,3,0.3,r)) #i chose 6 because of the 7th element referenced in the exercise

print(C[f74_noloop(nn,cp,3,0.3,r)]) #And these are the points


# %% Askhsh 75
import numpy as np

def f75_noloop(s,non_CP,eps,D):
    S_non = np.array([], dtype=int)
    ind = 0
    while ind != len(non_CP):
        s_elements = 0
        
        while s_elements != len(s):
            point = D[non_CP[ind]][s[s_elements]]
            if point<=eps:
    
                S_non = np.append(S_non,[non_CP[ind]])
                s_elements +=1
        
            else:
                s_elements+=1
        
        ind += 1
    
    return np.array(list(set(S_non)))
C = np.array([
    [0.3, 0.5],
    [0.9, 0.2],#1
    [0.8, 0.2],#2
    [0.5, 0.3],
    [0.6, 0.5],#4
    [0.3, 0.4],#5
    [0.3, 0. ],#6
    [0.3, 0.2],
    [0.2, 0.1],#8
    [0. , 0.5]])#9

nn = f71_noloop(C)
cp, non_cp = f73_noloop(C,3,0.3)

set_75 =f74_noloop(nn,cp,3,0.3,6)

print(f75_noloop(set_75,non_cp,0.3,nn))

print(C[f75_noloop(set_75,non_cp,0.3,nn)]) #and these are the points


# %% Askhsh 76 
import numpy as np
import random
def f76_noloop(C,minPTS,eps):
    dists = f71_noloop(C)
    cp, non_cp = f73_noloop(C,minPTS,eps)
    cluster_name = 1
    clusters = {}
    clusters_ind = {}
    r = random.randint(0,len(C)-1)
    check = False
    while True:
        condition = 0
        S = f74_noloop(dists,cp,minPTS,eps,r)
        S_points = C[f74_noloop(dists,cp,minPTS,eps,r)]
        T = f75_noloop(S,non_cp,eps,dists)
        T_points = C[f75_noloop(S,non_cp,eps,dists)]
        
        clusters_ind[cluster_name] = {"S": S, "T": T} #replace with _points
        clusters[cluster_name] = {"S": S_points, "T": T_points} #replace with _points
        while condition != len(clusters):
            condition += 1
            point_diff = list(set(cp) - set(clusters_ind[condition]["S"]) - set(clusters_ind[condition]["T"]))#This way we make sure the point is not contained in either the S set and the T set
            if len(point_diff):
                print(point_diff)
                r = random.randint(0,len(point_diff)-1)
                cluster_name += 1
                check = False
                break
            else:
                check = True
        if check:
            return clusters
  
C = np.array([
    [0.3, 0.5],
    [0.9, 0.2],
    [0.8, 0.2],
    [0.5, 0.3],
    [0.6, 0.5],
    [0.3, 0.4],
    [0.3, 0. ],
    [0.3, 0.2],
    [0.2, 0.1],
    [0. , 0.5]])


(f76_noloop(C,3,0.3))
# %% Askhsh 77

import numpy as np
import matplotlib.pyplot as plt 

def create_blob():
    N = 1000
    mean_1 = np.array([2,3])
    cov_1 = np.array([[2.0, 0.3], [0.3, 0.5]])
    A = np.random.multivariate_normal(mean_1, cov_1, N)
    return A

def f77():
    C = create_blob()
    
    clusters = f76_noloop(C,minPTS=3,eps=0.3)
    max_clust = max(clusters.values(), key = lambda x: (int(len(x["S"])) + int(len(x["T"]))))
    concat_max = np.concatenate((max_clust["S"], max_clust["T"]), axis = 0)

    remaining_ind = (C[:,None]!= concat_max).any(-1).all(1) #https://stackoverflow.com/questions/66674537/python-numpy-get-difference-between-2-two-dimensional-array
    remaining = C[remaining_ind]

    plt.scatter(max_clust["S"][:, 0], max_clust["S"][:, 1], color = 'green', s = 14)
    plt.scatter(max_clust["T"][:, 0], max_clust["T"][:, 1], color = 'red', s = 14)
    plt.scatter(remaining[:, 0], remaining[:, 1], color = 'black', s = 14)



print(f77())


# %% Askhsh 78 

import numpy as np
import matplotlib.pyplot as plt 

def f78(K, to_plot):
    
    A = create_blob() #create initial array
    
    move_B = np.array([K,0])#create a 1x2 array with K as the value of X
    move_B = np.tile(move_B, (len(A),1))#And now replicate it by N (where N is the length of A)
    
    move_C = np.full((len(A),len(A[0])),K) #create an identical matrix with NxM dimensions ( in our case 1000x2 ) and K values
    
    move_D = np.array([0,K])#same as move_B
    move_D = np.tile(move_D, (len(A),1))
    
    B = A - move_B #move the matrix A by -K on the x axis 
    C = A - move_C #move the matrix A by -K on the x and y axis
    D = A - move_D #move the matrix A by -K on the y axis
    
    if to_plot:
        plt.scatter(A[:, 0], A[:, 1], color = 'black', s = 14)
        plt.scatter(B[:, 0], B[:, 1], color = 'black', s = 14)
        plt.scatter(C[:, 0], C[:, 1], color = 'black', s = 14)
        plt.scatter(D[:, 0], D[:, 1], color = 'black', s = 14)
        
    return np.vstack((A,B,C,D))
    
    
print(f78(10,True))


# %% Askhsh 79
import numpy as np
import matplotlib.pyplot as plt

def f79(K, minPTS):
    C = f78(K,False)
    clusters = f76_noloop(C,minPTS, eps = 0.3)
    for i in range(len(clusters)):#We add +1 since we initialized cluster name with 1 and then adding +1 each time
        conc =  np.concatenate((clusters[i+1]["S"], clusters[i+1]["T"]), axis = 0)
        
        remaining_ind = (C[:,None]!= conc).any(-1).all(1) 
        remaining = C[remaining_ind]
        
        
        plt.scatter(clusters[i+1]["S"][:, 0], clusters[i+1]["S"][:, 1], color = 'green', s = 14)
        plt.scatter(clusters[i+1]["T"][:, 0], clusters[i+1]["T"][:, 1], color = 'red', s = 14)
        plt.scatter(remaining[:, 0], remaining[:, 1], color = 'black', s = 14)
        plt.title(f"K={K}")
        plt.ylabel(f"P={minPTS}")
    
    return plt.show

print(f79(10,3))




# %% Askhsh 80
import matplotlib.pyplot as plt #https://jakevdp.github.io/PythonDataScienceHandbook/04.08-multiple-subplots.html
plt.style.use('seaborn-white')
import numpy as np

# I have messed up somewhere and this runs for way too long so i cant even debug it. This is just an implementation of exercise 79 with 2 loops and since exercise 79 takes some times also does this.
#My plan was to use the ax to add a plot for each but i cant seem to make it work so i just removed the implementation.
#Now i am really curious as for how this is solved! Looks very interesting
def f80():
    fig, ax = plt.subplots(10, 10)
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    for i in range(1,11):
        for j in range(1,11):
            
            C = f78(i,False)
            clusters = f76_noloop(C,j, eps = 0.3)
            for k in range(len(clusters)):#We add +1 since we initialized cluster name with 1 and then adding +1 each time
                conc =  np.concatenate((clusters[k+1]["S"], clusters[k+1]["T"]), axis = 0)

                remaining_ind = (C[:,None]!= conc).any(-1).all(1) 
                remaining = C[remaining_ind]
                plt.scatter(clusters[k+1]["S"][:, 0], clusters[k+1]["S"][:, 1], color = 'green', s = 14)
                plt.scatter(clusters[k+1]["T"][:, 0], clusters[k+1]["T"][:, 1], color = 'red', s = 14)
                plt.scatter(remaining[:, 0], remaining[:, 1], color = 'black', s = 14)
                plt.title(f"K={i}")
                plt.ylabel(f"P={j}")


#print(f80()) #There is no need to run this it takes too long, i know its wrong, the results are not pretty either.
