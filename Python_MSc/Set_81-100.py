# %% _____ Askhsh 81 _____
import matplotlib.pyplot as plt
import numpy as np


array_81_x = np.array([[0,0,0,0],[0,0.7,0.7,0.7]])
array_81_y = np.array([[0,0,1,2],[2,0,1,2]])

plt.plot(array_81_x,array_81_y,color='black')
plt.xlim([-0.05, 1.4])




# %% _____ Askhsh 82 _____
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt

def f1(a,b):
    res1 = (a+b)/2
    return res1

def f2(A,B):
    res2 = f1(A[0], B[0]), f1(A[1],B[1])
    return res2

def f3(A,B):
    if A==B:
        res3 = 'Δεν υπάρχει ευθεία'
    elif A[0]==B[0]:
        res3 = 'Κάθετη'
    else:
        res3 = (B[1]-A[1])/(B[0]-A[0])
    return res3

def f4(A,B):
    res3 = f3(A,B)
    if res3 == 'Δεν υπάρχει ευθεία':
        res4 = res3
    elif A[1]==B[1]:
        res4 = 'Οριζόντια'
    else:
        res4 = -1/res3
    return res4


def f5(a,b,r,c):
    K = c[0]
    L = c[1]
    k1 = (-sqrt(-(a**2)*(K**2) + (a**2)*(r**2) - 2*a*b*K + 2*a*K*L - b**2 + 2*b*L - L**2 + r**2 ) - a*b + a*L + K ) / ( a**2 + 1)
    k2 = (sqrt(-(a**2)*(K**2) + (a**2)*(r**2) - 2*a*b*K + 2*a*K*L - b**2 + 2*b*L - L**2 + r**2 ) - a*b + a*L + K ) / ( a**2 + 1)
    l1 = a*k1 + b
    l2 = a*k2 + b
    if L != a*K + b:
        return 'Το κέντρο δεν ανήκει στην ευθεία'
    else:
        return (k1,l1) , (k2,l2)

def f82(A,B,dist):
    if A == B:
        return f'Λάθος'
    tup_mean = f2(A,B)
    slope = f4(A,B)
    constant = tup_mean[1] - (tup_mean[0]*slope)
    perpendicular = f5(slope,constant,dist,tup_mean)    
    #Another approach would be the below
    #array_82_x = np.array([[A[0],points[0][0]],[B[0],points[1][0]]]) 
    #array_82_y = np.array([[A[1],points[0][1]],[B[1],points[1][1]]]) 
    points_x = np.array([A[0],B[0]])
    points_y = np.array([A[1],B[1]])
    perpendicular_x = np.array([perpendicular[0][0],perpendicular[1][0]])
    perpendicular_y = np.array([perpendicular[0][1],perpendicular[1][1]])    
    plt.plot(points_x,points_y,color='black',marker = 'o', markersize = 3)
    plt.plot(perpendicular_x,perpendicular_y,color='blue',marker = 'o', markersize = 3)
    plt.plot(tup_mean[0],tup_mean[1],color ='red',marker = 'o', markersize = 3)
    #plt.xlim([-3, 8.5]) #Ta exw peiraksei gia na moiazoyn me askhsh
    #plt.ylim([0,10.5]) #Ta exw peiraksei gia na moiazoyn me askhsh

A = (1,2)
B = (3,9)
dist = 4
f82(A,B,dist)





# %% _____ Askhsh 83 _____
import matplotlib.pyplot as plt
import numpy as np
import random


def random_stocks():
    days = 200
    price_1 = 100
    price_2 = 100
    
    ret = [(price_1, price_2)]
    for x in range(days-1):
        price_1 += random.random() - 0.5
        price_2 += random.random() - 0.5
        ret.append((price_1, price_2))

    return ret


    
def plot_stocks():
    Y = random_stocks()
    X = list(range(len(Y)))
    stocks = max(Y, key = lambda x: abs(x[0]-x[1]))
    stock_index = Y.index(stocks)
    points_83_x = np.array([stock_index,stock_index])
    points_83_y = np.array([min(stocks), max(stocks)])
    fig, ax = plt.subplots()
    
    ax.plot(X, [y[0] for y in Y], c='blue')
    ax.plot(X, [y[1] for y in Y], c='red')
    ax.plot(points_83_x,points_83_y, c='black')
plot_stocks()


# %% _____ Askhsh 84 _____
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

def f84(x):
    frequency = dict(Counter(x))
    frequency_sum = {}
    for k,v in frequency.items():
        if not v in frequency_sum:
            frequency_sum[v] = [k]
        else:
            frequency_sum[v].append(k)
            
    max_frequency = max(frequency.values())
    axis_lim_x = max_frequency*2 -1
    axis_lim_y = max(x)+2
    for k1,v1 in frequency_sum.items():
        if k1%2 == 0:
            freq_index = []
            sth = 1
            for i in range(int(k1/2)): 
                freq_index.append(max_frequency - sth)
                freq_index.append(max_frequency + sth)
                sth += 2
            for age in v1:
                age_arr = np.full((k1,1),age)
                freq_arr = np.array(freq_index)
                plt.scatter(freq_arr,age_arr,color='black')
                
        else:    
            freq_index = [max_frequency]
            sth = 2
            if k1 == 1:
                freq_index = [max_frequency]
                
            for j in range(int(k1/2)):
                freq_index.append(max_frequency - sth)
                freq_index.append(max_frequency + sth)
                sth+=2

            for age in v1:
                age_arr = np.full((k1,1),age)
                freq_arr = np.array(freq_index)
                plt.scatter(freq_arr,age_arr,color='black')

    plt.xlim([0.7, axis_lim_x+0.3]) #Ta exw peiraksei gia na moiazoyn me askhsh
    plt.ylim([min(x)-1,axis_lim_y]) #Ta exw peiraksei gia na moiazoyn me askhsh

l = [25, 27, 23, 28, 22, 25, 35, 30, 20, 27, 34, 34, 32, 24, 23, 25, 22, 31, 26, 27, 24, 25, 28, 32, 34, 20]
f84(l) 


# %% _____ Askhsh 85 _____
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

def get_data_85():
    dfs = pd.read_html('https://en.wikipedia.org/wiki/List_of_minimum_annual_leave_by_country')
    df = dfs[0]
    
    ret = {
        x['Country']: {
            'adeies': int(x['Paid vacation days by year (five-day workweek)[1][2]']), 
            'argies': int(x['Total paid leave (five-day workweek)'])
        }  for x in df.to_dict('records') 
           if not pd.isna(x['Total paid leave (five-day workweek)'])
           if re.fullmatch(r'\d+', x['Paid vacation days by year (five-day workweek)[1][2]'])
           if re.fullmatch(r'\d+', x['Total paid leave (five-day workweek)'])
           }  
    return ret

countries = get_data_85()
countries_days = {}
countries_list = []
for i in countries.keys():
    countries_days[i] = (countries[i]['adeies'] + countries[i]['argies'])
countries_days = dict(sorted(countries_days.items(), key = lambda x: x[1]))

for j in countries_days.keys():
    countries_list.append(j)


fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
moires = 360/len(countries_days) #140 len
for m in range(len(countries_days)):
    country = countries_list[m]
    c_adeies = countries[country]['adeies'] # vd
    c_argies = countries[country]['argies'] #ph
    c_moires = moires*m
    
    ax.plot([c_moires/180 * np.pi, c_moires/180 * np.pi], [30, 30+c_adeies], c='red') #30+vd
    ax.plot([c_moires/180 * np.pi, c_moires/180 * np.pi], [30+c_adeies,30+c_adeies+c_argies], c='blue') #
    ax.axis('off')



# %% _____ Askhsh 86 _____
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import det
import random


def bezier_curve(p0, p1, p2, p3, t):
    x = (1 - t)**3 * p0[0] + 3 * (1 - t)**2 * t * p1[0] + 3 * (1 - t) * t**2 * p2[0] + t**3 * p3[0]
    y = (1 - t)**3 * p0[1] + 3 * (1 - t)**2 * t * p1[1] + 3 * (1 - t) * t**2 * p2[1] + t**3 * p3[1]
    return (x, y)

def create_random_curve():
    p0_1 = 10*random.random()
    p0_2 = 10*random.random()
    
    p1_1 = 10*random.random()
    p1_2 = 10*random.random()
    
    p2_1 = 10*random.random()
    p2_2 = 10*random.random()

    p3_1 = 10*random.random()
    p3_2 = 10*random.random()

    p0 = (p0_1, p0_2)
    p1 = (p1_1, p1_2)
    p2 = (p2_1, p2_2)
    p3 = (p3_1, p3_2)

    d = np.array([bezier_curve(p0, p1, p2, p3, x) for x in np.linspace(0,1,100)])

    return d

def left(l, b, x1, y1):

    x2 = 0
    y2 = l*x2 + b
    x3 = 10
    y3 = l*x3 + b
    K = np.array([[x2-x1, x3-x1], [y2-y1, y3-y1]])
    return np.sign(det(K))


def take_random_curve():
    data_arr = create_random_curve()
    arr_left = random.randint(0,40)
    arr_right = random.randint(41,99)
    if arr_right - arr_left < 10: #so we will have atleast a distance of 10
        arr_right += 10
    small_curve = data_arr[arr_left:arr_right]

    line_l = random.uniform(-2,2)
    line_b = random.uniform(-2,2)
    
    x_2 = 0
    x_3 = 10
    rectilinear_y_2 = line_l*x_2 + line_b
    rectilinear_y_3 = line_l*x_3 + line_b
    rectilinear_arr_x = np.array([x_2,x_3]) 
    rectilinear_arr_y= np.array([rectilinear_y_2,rectilinear_y_3])
    
    fig, ax = plt.subplots()
    ax.plot(data_arr[:, 0], data_arr[:, 1], '-', color='black', alpha = 0.05 ) # []
    ax.plot(data_arr[arr_left:arr_right, 0], data_arr[arr_left:arr_right, 1], '-', color='black' ) # []
    ax.plot(rectilinear_arr_x,rectilinear_arr_y, c = 'blue')
    point_sum = []
    point_mean = False
    
    for point in small_curve:
        if not point_mean: # compensate first point
            point_mean = True
            point_arr = np.array(point)
            continue
        
        intersect = left(line_l,line_b,point[0],point[1])
        point_sum.append(intersect)
        point_arr = np.vstack((point_arr, point))
        
        if not intersect:
            ax.plot(point[0],point[1],marker = 'D', color = 'red', markersize = 3)
        if len(point_sum) == 2:
            if not sum(point_sum):
                int_point = np.sum(point_arr, axis = 0, dtype = float)/2
                ax.plot(int_point[0],int_point[1],marker = 'D', color = 'red', markersize = 6)    
            point_sum = point_sum[::-1]
            point_sum.pop()
            point_arr =  np.flipud(point_arr)  
            point_arr = point_arr[0]

take_random_curve()


# %% _____ Askhsh 87 _____
import matplotlib.pyplot as plt
import numpy as np

def bezier_curve(p0, p1, p2, p3, t):
    x = (1 - t)**3 * p0[0] + 3 * (1 - t)**2 * t * p1[0] + 3 * (1 - t) * t**2 * p2[0] + t**3 * p3[0]
    y = (1 - t)**3 * p0[1] + 3 * (1 - t)**2 * t * p1[1] + 3 * (1 - t) * t**2 * p2[1] + t**3 * p3[1]
    return (x, y)

def create_random_curve():
    p0_1 = 10*random.random()
    p0_2 = 10*random.random()
    
    p1_1 = 10*random.random()
    p1_2 = 10*random.random()
    
    p2_1 = 10*random.random()
    p2_2 = 10*random.random()

    p3_1 = 10*random.random()
    p3_2 = 10*random.random()

    p0 = (p0_1, p0_2)
    p1 = (p1_1, p1_2)
    p2 = (p2_1, p2_2)
    p3 = (p3_1, p3_2)

    d = np.array([bezier_curve(p0, p1, p2, p3, x) for x in np.linspace(0,1,100)])

    return d



def create_art():
    art_data = create_random_curve()
    fig, ax = plt.subplots()
    for i in range(4):
            
        most_left = min(art_data, key = lambda x: x[0])[0]
        most_down = min(art_data, key = lambda x : x[1])[1]

        y_axi = []
        x_axi = []
        sym_x_axi = []
        sym_y_axi = []
        for i in art_data:
            x_axi.append(i[0])
            y_axi.append(i[1])    
            val_x = i[0]
            val_y = i[1]
            sym_x_axi.append(2*most_left-val_x)
            sym_y_axi.append(2*most_down-val_y)

        resh_symy_x = np.reshape(sym_x_axi,(-1,1))
        resh_symy_y = np.reshape(y_axi,(-1,1))
        sym_arr_B  = np.concatenate([resh_symy_x,resh_symy_y], axis = 1) # B
        resh_symx_x = np.reshape(x_axi,(-1,1))
        resh_symx_y = np.reshape(sym_y_axi,(-1,1))
        sym_arr_C  = np.concatenate([resh_symx_x,resh_symx_y], axis = 1) # C
        most_left_C = min(sym_arr_C, key = lambda x: x[0])[0]

        sym_x_axi_C = []
        for j in art_data:
            val_y_C = j[0]
            sym_x_axi_C.append(2*most_left_C-val_y_C)

        resh_symc_x = np.reshape(sym_x_axi_C,(-1,1))
        resh_symc_y = np.reshape(resh_symx_y,(-1,1))
        sym_arr_D  = np.concatenate([resh_symc_x,resh_symc_y], axis = 1) # C
        ax.plot(art_data[:, 0], art_data[:, 1], '-', color='black') # A
        ax.plot(sym_arr_B[:, 0], sym_arr_B[:, 1], '-', color='black') # B
        ax.plot(sym_arr_C[:, 0], sym_arr_C[:, 1], '-', color='black') # C
        ax.plot(sym_arr_D[:, 0], sym_arr_D[:, 1], '-', color='black') # D
        art_data = np.concatenate([art_data,sym_arr_B,sym_arr_C,sym_arr_D], axis=0)

create_art()


# %% _____ Askhsh 88 _____
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binom 

def B_wins(N, p): 
    if N%2 == 0:
        return binom.cdf(N//2-1, N, p)
    
    return binom.cdf(N//2, N, p)


def min_games(p, a): 
    for N in range(3,10_000):
        if B_wins(N, p)<a:
            return N 

def f88():
    range_a = (np.array(list(range(1,10,1))))/100 
    range_p = (np.array(list(range(530,650,5))))/1000 
    fig, ax = plt.subplots()
    colors =plt.cm.viridis(np.linspace(0,1,9))
    color_count = 0
    for a in range_a:
        
        results = []
        for p in range_p:
            games = min_games(p,a)
            results.append(games)
        results_arr = np.array(results)
        ax.plot(range_p, results_arr, '-', c=colors[color_count])
        color_count += 1

f88()
# %% _____ Askhsh 89 _____
import numpy as np
import pandas as pd
from collections import Counter
import re


def get_df_89():
    return pd.read_csv('https://www.dropbox.com/s/5e4btcptal5pp3w/ask_89.csv?dl=1')

def f89():
    data_89 = get_df_89()
    tr_list = (data_89['GENE.NAME']).tolist()
    tr_count = sorted((dict(Counter(tr_list))).items(), key = lambda x: x[1], reverse = True)
    return tr_count[0][0]
f89()

# %% _____ Askhsh 90 _____
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import Counter


def f90():
    data_90 = get_df_89()
    tr_count = Counter((data_90['GENE.NAME']).tolist())
    tr_arr = []

    for i in data_90['GENE.NAME'].tolist():
        tr_arr.append(tr_count[i])

    tr_col = np.column_stack([tr_arr])
    tr_col_pd = pd.DataFrame(tr_col, columns = ['TR.COUNT'])
    data_90['TR.COUNT'] = tr_col
    data_90['AV.EXON.L'] = data_90['EXONS.L']/data_90['N.EXONS']
    data_90['AV.INTRONS.L'] = data_90['INTRONS.L']/data_90['N.INTRONS']
    data_90['AV.EXON.DIV.INTRON'] = data_90['AV.EXON.L']/data_90['AV.INTRONS.L']
    data_90['MEAN.AV.EX.IN.DIV'] = data_90['AV.EXON.DIV.INTRON']/data_90['TR.COUNT']

    data_90_drop = data_90.dropna()
    df_sorted = data_90_drop.sort_values('MEAN.AV.EX.IN.DIV', ascending=True)
    plot_df = df_sorted[['GENE.NAME', 'MEAN.AV.EX.IN.DIV']].copy()

    plot_df.plot(kind='line', logy = True)

f90()
# %% _____ Askhsh 91 _____

import math
def is_prime(n):
    for x in range(2, int(math.sqrt(n))+1):
        if n% x == 0:
            return False
    return True


def prime_generator():
    yield 2
    current = 3
    while True:
        if is_prime(current):
            yield current
        current += 2

def gen_91(n):
    count = 1
    a = n
    b = n*2   
    while b-a == n:
        rng = list(range(a,b+1))
        for j in prime_generator():
            if j in rng:
                a += 1
                b += 1
                break
            if j > b:
                yield(a,b)
                a += 1
                b += 1
                break
            continue

for i,pair in enumerate(gen_91(3)):
    if i>10:
        break
    print (pair) # to evala edw gia na bgazei to idio apotelesma me esas alliws vgazei +1 tuple


# %% _____ Askhsh 92 _____
from collections import Counter

def greek_word_gen():
    fn = 'words_greek_normalized.txt'
    with open(fn) as f:
        for l in f:
            if len(l) == 6:
                yield l.strip()


def leksli(gen,tried,pattern):
    patrn_pos = list(enumerate(zip(list(tried),list(pattern))))
    for i in gen:
        ind = 0
        word_list = list(i)
        word_cnt = Counter(i)
        for j in patrn_pos:
            existence = int(j[1][1])
            
            if not existence:
                if str(j[1][0]) in i:
                    break 
                
                else:
                    ind += 1
                    if ind == 5:
                        yield i
                
                    
            if existence == 1: 
                if str(j[1][0]) not in i:
                    break
                
                try: 
                    excpt = word_list.index(j[1][0],ind)
                except:
                    excpt = word_list.index(j[1][0])
                
                 
                if word_cnt[str(j[1][0])] == 2: 
                    if int(j[0]) == excpt:
                        break
                    else:
                        ind += 1
                        if ind == 5:
                            yield i
                if int(j[0]) == excpt:
                    break
                    
                else:
                    ind += 1
                    if ind == 5:
                        yield i

            
            if existence == 2:

                    
                if str(j[1][0]) not in i:
                    break
                
                
                try: 
                    excpt = word_list.index(j[1][0],ind)
                except:
                    excpt = word_list.index(j[1][0])
                
                
                if word_cnt[str(j[1][0])] == 2: 
                    if int(j[0]) != excpt:
                        break
                    else:
                        ind += 1
                        if ind == 5:
                            yield i
                
                if int(j[0]) != excpt:
                    break
                
                else:
                    ind += 1
                    if ind == 5:
                        yield i
g1 = greek_word_gen()
g2 = leksli(g1, 'στανη', '01000')
#g3 = leksli(g2, 'χωροσ', '10002')
#g4 = leksli(g3, 'ιχθυσ', '02202')              

for wrds in g4:
    print(wrds)
# %% ----- This is Extra ------
# this can be used to take a random word to start your own leksli
import random

def random_word():
    gen = greek_word_gen()  
    len_gen = len(list(gen)) # gia kapoio logo prepei na ksanaorisw to gen an kanw ayto
    gen = greek_word_gen()  
    random_word = random.randint(1,len_gen+1) 
    for i in range(1,len_gen+2):
        chosen_word = next(gen)
        if i == random_word:
            yield chosen_word
            break
    
the_word = (list(random_word()))[0]

# %%  ----- This is Extra ------
# this can be used to play your own leksli, by inputting a word and it gives back a string
def leksli_start_game(gen,tried): 
    attempt = list(tried)
    result = ''
    chosen_word = gen
    ind = 0
    for letter in attempt:
        if letter not in chosen_word:
            result = result+'0'
            ind+=1
            continue
        
        if letter in chosen_word:    
            if ind == int(list(enumerate(chosen_word))[ind][0]) and letter == chosen_word[ind]:
                result = result+'2'
                ind+=1
                continue
            else:
                result = result+'1'
                ind+=1
    yield result

first_try = leksli_start_game(the_word,'πλωτε')

for wrds in first_try:
    print(wrds)

# %% _____ Askhsh 93 _____
import re

def gen93(ex_set):
    regex = re.compile('^# exercise [\d]+')
    exc_dict = {}
    num = 0
    quot_sign = 0
    quotation_start = '```\n'
    quotation_end = '```'
    with open(ex_set) as file:
        for i in file:
            whatline = str(regex.findall(i))
            whatline = re.sub('[^\d]','',whatline)
            if whatline:
                num = int(whatline)
                exc_dict[num] = [quotation_start]
                if not quot_sign:
                    continue
                if quot_sign != int(whatline):
                    exc_dict[quot_sign].append(quotation_end)
                continue
            if not num:
                continue    
            else:
                if i == '\n':
                    continue
                quot_sign = num
                exc_dict[num].append(i)
    exc_dict[quot_sign].append(quotation_end)
    exc_dict = dict(sorted(exc_dict.items(), key = lambda x : x[0]))
    for k,v in exc_dict.items():
        solution = ''.join(exc_dict[k])
        yield (k, solution)


corrector = gen93('file_93.txt')


# %% _____ Askhsh 94 _____
#I have implemented 2 possible solutions, as specified on the right, but i have kept 
import re
def line_generator(fn):
    with open(fn) as f:
        for l in f:
            yield l.strip()

def gen94(fasta):
    regex_q = re.compile('^@.*') # for fast q
    rejectex = re.compile('^time.*') # for fast q
    #regex = re.compile('^>.*') #for fasta
    #exc_dict = {} #for fasta
    #current = '' #for fasta
    fastq_dict = {} #for fast q
    current_q = '' # for fast q
    mediator = '\n+\n' #for fast q
    fasta_gen = line_generator(fasta)
    for i in fasta_gen:
        
        reject = rejectex.findall(i)
        if reject: #for fast q
            continue
        
        #whatline = regex.findall(i) #for fasta
        fastqline = regex_q.findall(i) #for fast q
        
        if fastqline: #whatline: #for fasta
            #current = whatline[0] #for fasta
            #exc_dict[current] = [] #for fasta
            current_q = fastqline[0] #for fast q
            fastq_dict[current_q] = [] #for fastq
            continue
        
        else:
            #exc_dict[current].append(i) #for fasta
            if i == '+':
                fastq_dict[current_q].append(mediator) #for fastq
                continue
            else:
                fastq_dict[current_q].append(i)
    
    for k,v in fastq_dict.items(): #exc_dict: #for fasta
        #solution = ''.join(exc_dict[k]) #for fasta
        solution_q = ''.join(fastq_dict[k]) #for fastq
        yield (k, solution_q)

file_g = gen94('small.file_94.fastq')


# %% _____ Askhsh 95 _____
class Line:
    def __init__(self,a,b):
        self.a = a
        self.b = b
        
    def __str__(self):
        return f'y = {self.a}*x + {self.b}'

l = Line(3,5)
print(l)


# %% _____ Askhsh 96 _____
class Line:
    def __init__(self,a,b):
        self.a = a
        self.b = b
        
    def __str__(self):
        return f'y = {self.a}*x + {self.b}'

    def is_parallel(self,other):
        return self.a == other.a
        
l_1 = Line(3,5)
l_2 = Line(3,7)
l_3 = Line(4,5)

print(l_1.is_parallel(l_2))

# %% _____ Askhsh 97 _____

class Line:
    def __init__(self,a,b):
        self.a = a
        self.b = b
        
    def __str__(self):
        return f'y={self.a}*x + {self.b}'

    def is_parallel(self,other):
        return self.a == other.a

    def __or__(self,other):
        return self.a == other.a
    
l_1 = Line(3,5)
l_2 = Line(3,7)
l_3 = Line(4,5)

print(l_1 | l_3)
# %% _____ Askhsh 98 _____
from itertools import combinations

class Line:
    def __init__(self,a,b):
        self.a = a
        self.b = b
        
    def __str__(self):
        return f'y={self.a}*x + {self.b}'

    def is_parallel(self,other):
        return self.a == other.a

    def __or__(self,other):
        return self.a == other.a
    
class LineCollection:
    
    def __init__(self,lines):
        self.lines = lines
        self.linelist = [str(lein) for lein in lines]
        self.combs = [(a,b) for a,b in combinations(self.linelist, 2) if a[2]==b[2]]
    
    def __iter__(self):
        self.current_index = 0
        return self
    
    def __next__(self):
        if self.current_index < len(self.combs):
            line_1,line_2 = self.combs[self.current_index]
            self.current_index += 1
            return line_1,line_2
        
        raise StopIteration


l = [
 (1, 96),
 (2, 92),
 (1, 11),
 (1, 25),
 (1, 90),
 (3, 41),
 (3, 75),
 (2, 37),
 (1, 43),
 (3, 68),
 (1, 32),
 (2, 4),
 (3, 25),
 (2, 2),
 (2, 70),
 (1, 98),
 (3, 51),
 (3, 65),
 (2, 52),
 (2, 85),
]

lines = [Line(x[0], x[1]) for x in l]

lc = LineCollection(lines)

for a,b in lc:
    print(f'{a} is parallel to {b}')


# %% _____ Askhsh 99 _____
def f99(A,B):

    if not ( type(A) and type(B) ) is tuple:
        raise TypeError('Λάθος τύπος')

    if ( len(list(A))  or len(list(B)) ) != 2:
            raise TypeError('Λάθος τύπος')
    
    if not ( type(A[0]) or type(A[1]) or type(B[0]) or type(B[1]) ) is ( int or float or tuple ):
        raise TypeError('Λάθος τύπος')
        
    if not (A[0] - B[0]) + (A[1] - B[1]):
        raise ValueError('Δεν υπάρχει ευθεία')

    if not (A[1]-B[1]):
        raise ValueError('Οριζόντια ευθεία')
    
    else:
        return f'Η απάντηση είναι {-1/(B[1]-A[1])/(B[0]-A[0])}'


# %% _____ Askhsh 100 _____
#Epistrefei to error alla meta den to epistrefei sth synarthsh den katalabainw


def f_100(a,b):
    try:
        f99(a,b)
    except Exception as exc:
        return exc
    return f99(a,b)

f_100( 5, 7.8)


# %% 
# Epistrefei oles tis times alla oxi me to error
def f_100(a,b):
    
    try:
        print(f99(a,b))

    except Exception as e:
        print(e)
    
