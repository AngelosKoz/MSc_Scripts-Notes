# %%
#Askhsh 21

def f21(x):
    res21 = 0
    vowels = 'AaEeIiOoUu'

    for i in range(len(x)):
        for j in vowels:           
            if x[i][0][0] == j:
                for k in vowels:
                    if x[i][1][0] == k:
                        res21 = 1
                        break
    
    if res21:
        return 2%2==0
    if not res21:
        return 2%2 == 1

l  =[('Alexandros', 'Kanterakis'), ('Spiros', 'Mpimpilas'), ('Maria', 'Aliferh')]
print (f21(l)) # Τυπώνει False (κανένα ζευγάρι δεν έχει και τα 2 πρώτα γράμματα φωνήεντα )

l  =[('Alexandros', 'Kanterakis'), ('Spiros', 'Mpimpilas'), ('Eleftheria', 'Arvanitaki'), ('Maria', 'Aliferh')]
print(f21(l)) # Τυπώνει True 
    

# %%
#Αskhsh 21 allos tropos
def f21(x):
    vowels = 'AaEeIiOoUu'
    for i in range(len(x)):
        for j in vowels:           
            if x[i][0][0] == j:
                for k in vowels:
                    if x[i][1][0] == k:
                        return 2%2==0
                        break
    
    else:
        return 2%2==1

l  =[('Alexandros', 'Kanterakis'), ('Spiros', 'Mpimpilas'), ('Maria', 'Aliferh')]
print (f21(l)) # Τυπώνει False (κανένα ζευγάρι δεν έχει και τα 2 πρώτα γράμματα φωνήεντα )

l  =[('Alexandros', 'Kanterakis'), ('Spiros', 'Mpimpilas'), ('Eleftheria', 'Arvanitaki'), ('Maria', 'Aliferh')]
print(f21(l)) # Τυπώνει True 

# %%
#Askhsh 22

def f22(x):
    res22 = 0
    vowels = 'AaEeIiOoUu'

    for i in range(len(x)):
        for j in vowels:           
            if x[i][0][0] == j:
                for k in vowels:
                    if x[i][1][0] == k:
                        res22 += 1
                        continue
    return res22

l = [('Alexandros', 'Kanterakis'), ('Spiros', 'Mpimpilas'), ('Maria', 'Aliferh')]
print (f22(l))

l = [('Alexandros', 'Kanterakis'), ('Spiros', 'Mpimpilas'), ('Eleftheria', 'Arvanitaki'), ('Maria', 'Aliferh')]
print (f22(l))

l = [('Alexandros', 'Kanterakis'), ('Spiros', 'Mpimpilas'), ('Eleftheria', 'Arvanitaki'), ('Maria', 'Aliferh'),  ('Aristotelhs', 'Onassis')]
print (f22(l))

# %%
#Askhsh 23

def f23(a,b):
    res23 = 0
    vowels = 'AaEeIiOoUu'
    x = list(zip(a,b))
    for i in range(len(x)):
        for j in vowels:           
            if x[i][0][0] == j:
                for k in vowels:
                    if x[i][1][0] == k:
                        res23 += 1
                        continue
    return res23

a = ['Alexandros', 'Spiros', 'Maria']
b = ['Kanterakis', 'Mpimpilas', 'Aliferh']
print (f23(a,b)) # Τυπώνει 0 (κανένα ζευγάρι δεν έχει και τα 2 πρώτα γράμματα φωνήεντα )


a = ['Alexandros', 'Spiros', 'Eleftheria', 'Maria']
b = ['Kanterakis', 'Mpimpilas', 'Arvanitaki', 'Aliferh']
print (f23(a,b)) # Τυπώνει 1 

a = ['Alexandros', 'Spiros', 'Eleftheria', 'Maria', 'Aristotelhs']
b = ['Kanterakis', 'Mpimpilas', 'Arvanitaki', 'Aliferh', 'Onassis']
print (f23(a,b)) # Τυπώνει 2

# %%
#Askhsh 24

def f24(x):
    vowels = 'AaEeIiOoUu'
    new_list=[]
    for i in range(len(x)):
        for j in vowels:           
            if x[i][0][0] == j:
                for k in vowels:
                    if x[i][1][0] == k:
                        new_list.append(x[i][1])
    if new_list:
        return min(new_list)                 
                        
l = [('Alexandros', 'Kanterakis'), ('Spiros', 'Mpimpilas'), ('Maria', 'Aliferh')]
print (f24(l)) # Τυπώνει `None` (κανένα ζευγάρι δεν έχει και τα 2 πρώτα γράμματα φωνήεντα )

l = [('Alexandros', 'Kanterakis'), ('Spiros', 'Mpimpilas'), ('Aristotelhs', 'Onassis'), ('Maria', 'Aliferh')]
print (f24(l)) # Τυπώνει `Onassis`

l = [('Alexandros', 'Kanterakis'), ('Spiros', 'Mpimpilas'), ('Eleftheria', 'Arvanitaki'), ('Maria', 'Aliferh'),  ('Aristotelhs', 'Onassis')]
print (f24(l)) # Τυπώνει 'Arvanitaki'  

# %%
#Askhsh 25
import math
def f25(x):
    return math.exp(x)-3*x

def sol(f25,a,b):
    lpv = a
    while lpv<b:
        res25 = f25(lpv)
        prev = f25(lpv-0.01)
        if res25 < 0 < prev or prev < 0 < res25:
            break
        lpv+=0.01
    return lpv
        
print (sol(f25,0,1))
print (sol(f25,1,2))  

# %%
#Askhsh 26

def f26():
    capital=10_000
    years = 0
    while capital < 20_000:
        years += 1
        capital -= 1_000
        capital = capital + capital*(15/100) 
        if capital >= 20_000:
            break
    return years
print(f26())        

# %%
#Askhsh 27

def f27():
    n = 1
    thresh = set(range(0,10))
    while n >= 1:
        res27 = 3**n
        cap = {int(i) for i in str(res27)}
        if cap == thresh:
            break
        n+=1
    return n
print(f27())

# %%
#Askhsh 27: Allos tropos

def f27():
    n = 1
    cap = 0
    thresh = set(range(0,10))
    while cap != thresh:
        res27 = 3**n
        cap = {int(i) for i in str(res27)}
        n+=1
    return n-1
print(f27())

# %%
#Askhsh 28

import math

def f28():
    n = 0
    k = 2
    while True:
        res28 = True    
        for i in range(2, int(math.sqrt(k) + 1)):
            if k % i == 0: 
                res28 = False
                break
        if res28:
            p = k
            n += 1
            if n/p <= 0.1:
                break
        k += 1
    return (n, p, n/p,)  
print(f28())

# %%
#Askhsh 29
import math

def f29():
    k = 2
    consec = 0
    list_29 = []
    got_it = False
    while True:
        res29 = True    
        for i in range(2, int(math.sqrt(k) + 1)):
            if k % i == 0: #den einai prwtos
                res29 = False
                break
        if res29: #einai prwtos
            p = k
            consec2 = consec + 1
            list_29.append(p)
            len_29 = len(list_29)
            if consec2 - p == 0: #den eiai synexomenos
                list_29 = []
            if len_29 == 7: #einai lista me 7 stoixeia
                prime = sum(list_29)
                if prime > 1:
                        for j in range(2, int(prime/2)+1):
                            if (prime % j) == 0: #den einai prwtos
                                list_29.remove(list_29[0])
                                got_it = False
                                break
                            elif (prime % j ) != 0: #einai prwtos
                                got_it = True
            consec = p 
        if got_it: #einai prwtos
            break
        k += 1
    return [sum(list_29), list_29]
print(f29())

# %%
#Askhsh 30

def f30(x):
    topstar = 0
    bot = 2
    botstar = 1
    pattern = ""
    blank = ""
    for i in range(1, x + 1): 
        for j in range (1, (x - i) + 1): 
            pattern = " ".join([pattern, blank])         
        while topstar != (( 2 * i ) - 1):
            pattern = "*".join([pattern, blank]) 
            topstar += 1
        topstar = 0   
        pattern = "\n".join([pattern, blank]) 
    for k in range(1, x): 
        for l in range (1, bot):
            pattern = " ".join([pattern, blank])  
        bot += 1	  
        while botstar <= (2 * (x - k) - 1): 
            pattern = "*".join([pattern, blank]) 
            botstar += 1
        botstar = 1	
        pattern = "\n".join([pattern, blank]) 
        if botstar > (2 * (x - k) - 1):
            break
    return pattern
print(f30((11)))



# %%
#Askhsh 31

import string

def letters():
    return set(string.ascii_lowercase)

def f31(ls):
    list_31 = []
    pool = letters()
    for i in range(len(ls)):
        for j in ls[i]:
            list_31.append(j)
            res31 = set(list_31)
    return pool-res31
            
ls = [
        'agapi',
        'olympia',
        'andreas',
        'nikos',
        'elena',
        'aggelos',
        'anna',
        'katerina',
    ]
print (f31(ls)) # Τυπώνει: {'b', 'c', 'f', 'h', 'j', 'q', 'u', 'v', 'w', 'x', 'z'}

# %%
#Askhsh 32

def f32(x):
    val_32 = list(x.values())
    key_32 = list(x.keys())
    list_32 = []
    for i in val_32:
        frst = i[0]
        scnd = i[1]
        list_32.append((scnd,frst))
    res32 = dict(zip(key_32, list_32))
    return res32

d = {
 '1': (5, 4),
 '2': (6, 3),
 '3': (2, 9),
 '4': (4, 9),
 '5': (9, 2),
 '6': (10, 10),
 '7': (8, 6),
 '8': (10, 4),
 '9': (4, 3),
 'X': (1, 6),
 }

print(f32(d))

# %%
#Askhsh 33
def f33(n):
    key_33 = list(range(1,n+1))
    val_33 = []
    div_33 = 0
    
    for j in range(1,n+1):
        for i in range(2,j):
            if j%i == 0:
                div_33 += 1
        val_33.append(div_33)
        div_33 = 0
    res33 = dict(zip(key_33, val_33))
    return res33

print(f33(50))          
# %%
#Askhsh 34

def gene_names(n):
    names_34 = []
    for i in n:
        for k,v in i.items():
            names_34.append(v)
            break
    return names_34

def f34(x):
    start = []
    end = []
    sum_start = []
    sum_end = []
    key_34 = gene_names(x)
    for i in range(len(x)):
        start.append(list(x[i].values())[1])
        end.append(list(x[i].values())[2])
    for s in range(len(start)):
        sum_start.append(sum(start[s]))
    for e in range(len(end)):
        sum_end.append(sum(end[e]))
    val_34 = [k - l for k, l in zip(sum_end, sum_start)]
    res34 = dict(zip(key_34, val_34))
    return res34

genes = [
 {'name': 'VLDZ',
  'start': [1335, 2287, 3395, 4344],
  'end': [1568, 2727, 3976, 4864]
 },
 {'name': 'SUVM',
  'start': [1014, 2064, 3131, 4335],
  'end': [1608, 2670, 3583, 4775]
 },
 {'name': 'AMPM',
  'start': [1386, 2305, 3306, 4010, 5394, 6260],
  'end': [1972, 2786, 3710, 4601, 5871, 6702]
 },
 {'name': 'ZNYP',
  'start': [1025, 2272, 3185, 4225, 5159, 6362],
  'end': [1942, 2679, 3552, 4584, 5961, 6725]
 },
 {'name': 'DVLY',
  'start': [1221, 2330, 3013, 4240, 5386, 6230],
  'end': [1533, 2571, 3737, 4747, 5895, 6922]
 },
 {'name': 'SSMU', 
  'start': [1096, 2253, 3170],
  'end': [1923, 2801, 3611]},
 {'name': 'KAEX',
  'start': [1013, 2010, 3398, 4292],
  'end': [1973, 2771, 3864, 4982]
 },
 {'name': 'ZTDU',
  'start': [1177, 2250, 3225, 4132, 5040, 6099],
  'end': [1676, 2682, 3915, 4688, 5723, 6956]
 },
 {'name': 'NQSP',
  'start': [1399, 2377, 3163, 4229],
  'end': [1513, 2589, 3695, 4550]},
 {'name': 'MXWY', 
  'start': [1082, 2052, 3292], 
  'end': [1726, 2581, 3865]
 }
]

print(f34(genes))


# %%
#Askhsh 35

def f35(d):
    list_35 = []
    for i in d:
        for j in d[i]:
            res35 = j
            if d[i][j] % 2 == 1:
                list_35.append(res35)
                
    
    return list_35

d = {
 'fwdy': {'yldh': 5, 'ujaw': 16},
 'zeiw': {'oyof': 14, 'wnbt': 4, 'sbzh': 16, 'yoke': 13},
 'mjnf': {'iatr': 12, 'owuc': 13, 'izvb': 12, 'axdk': 8},
 'vpvp': {'pqwe': 2, 'gtil': 19, 'qzjs': 15, 'lrgk': 19},
 'flkz': {'ykfj': 20},
 'vawu': {'anbv': 3, 'jfmq': 11, 'gpet': 11, 'bnru': 17},
 'oyzw': {'mclq': 19},
 'raip': {'jrzs': 4},
 'wxwv': {'vzoa': 16, 'inuw': 9},
 'rfcj': {'bmdi': 8, 'zurf': 2, 'sbcj': 17}
 }

print(f35(d))

# %%
#Askhsh 36

def f36(d):
    list_36 = []
    sum_36 = []
    for i in d:
        res36 = i
        for j in d[i]:
            val_36 = d[i][j]
            sum_36.append(val_36)
        if sum(sum_36) % 2 == 1:
            list_36.append(res36)
        sum_36 = []
    return list_36

d = {
 'fwdy': {'yldh': 5, 'ujaw': 16},
 'zeiw': {'oyof': 14, 'wnbt': 4, 'sbzh': 16, 'yoke': 13},
 'mjnf': {'iatr': 12, 'owuc': 13, 'izvb': 12, 'axdk': 8},
 'vpvp': {'pqwe': 2, 'gtil': 19, 'qzjs': 15, 'lrgk': 19},
 'flkz': {'ykfj': 20},
 'vawu': {'anbv': 3, 'jfmq': 11, 'gpet': 11, 'bnru': 17},
 'oyzw': {'mclq': 19},
 'raip': {'jrzs': 4},
 'wxwv': {'vzoa': 16, 'inuw': 9},
 'rfcj': {'bmdi': 8, 'zurf': 2, 'sbcj': 17}
 }

print(f36(d))
# %%
#Askhsh 37: Me all()

def f37(d):
    list_37 = []
    res_37 = []
    for i in d:
        for k,v in d[i].items():
            val_37 = v
            if v % 2 == 0:
                list_37.append(0)
            if v % 2 == 1:
                list_37.append(1)
        if all(list_37):
            res_37.append(i)
        list_37 = []
    return res_37
            


d = {
 'fwdy': {'yldh': 5, 'ujaw': 16},
 'zeiw': {'oyof': 14, 'wnbt': 4, 'sbzh': 16, 'yoke': 13},
 'mjnf': {'iatr': 12, 'owuc': 13, 'izvb': 12, 'axdk': 8},
 'vpvp': {'pqwe': 2, 'gtil': 19, 'qzjs': 15, 'lrgk': 19},
 'flkz': {'ykfj': 20},
 'vawu': {'anbv': 3, 'jfmq': 11, 'gpet': 11, 'bnru': 17},
 'oyzw': {'mclq': 19},
 'raip': {'jrzs': 4},
 'wxwv': {'vzoa': 16, 'inuw': 9},
 'rfcj': {'bmdi': 8, 'zurf': 2, 'sbcj': 17}
 }

print(f37(d))

# %%
#Askhsh 37: Allos Tropos

def f37(d):
    list_37 = []
    allsingle = True
    for i in d:
        for j in d[i]:
            val_37 = d[i][j]
            if val_37 % 2 == 0:
                allsingle = False
                break
            if val_37 % 2 == 1:
                allsingle = True
        if allsingle:
            list_37.append(i)
    return list_37


d = {
 'fwdy': {'yldh': 5, 'ujaw': 16},
 'zeiw': {'oyof': 14, 'wnbt': 4, 'sbzh': 16, 'yoke': 13},
 'mjnf': {'iatr': 12, 'owuc': 13, 'izvb': 12, 'axdk': 8},
 'vpvp': {'pqwe': 2, 'gtil': 19, 'qzjs': 15, 'lrgk': 19},
 'flkz': {'ykfj': 20},
 'vawu': {'anbv': 3, 'jfmq': 11, 'gpet': 11, 'bnru': 17},
 'oyzw': {'mclq': 19},
 'raip': {'jrzs': 4},
 'wxwv': {'vzoa': 16, 'inuw': 9},
 'rfcj': {'bmdi': 8, 'zurf': 2, 'sbcj': 17}
 }

print(f37(d))


# %%
#Askhsh 38: Me any()

def f38(d):
    list_38 = []
    res_38 = []
    for i in d:
        for k,v in d[i].items():
            val_38 = v
            if v % 2 == 0:
                list_38.append(0)
            if v % 2 == 1:
                list_38.append(1)
        if any(list_38):
            list_38 = []
        else:            
            res_38.append(i)
    return res_38


d = {
 'fwdy': {'yldh': 5, 'ujaw': 16},
 'zeiw': {'oyof': 14, 'wnbt': 4, 'sbzh': 16, 'yoke': 13},
 'mjnf': {'iatr': 12, 'owuc': 13, 'izvb': 12, 'axdk': 8},
 'vpvp': {'pqwe': 2, 'gtil': 19, 'qzjs': 15, 'lrgk': 19},
 'flkz': {'ykfj': 20},
 'vawu': {'anbv': 3, 'jfmq': 11, 'gpet': 11, 'bnru': 17},
 'oyzw': {'mclq': 19},
 'raip': {'jrzs': 4},
 'wxwv': {'vzoa': 16, 'inuw': 9},
 'rfcj': {'bmdi': 8, 'zurf': 2, 'sbcj': 17}
 }

print (f38(d))

# %%
#Askhsh 38: Allos tropos

def f38(d):
    list_38 = []
    for i in d:
        for j in d[i]:
            val_38 = d[i][j]
            if val_38 % 2 == 1:
                nosingle = False
                break
            if val_38 % 2 == 0:
                nosingle = True
        if nosingle:
            list_38.append(i)
    return list_38


d = {
 'fwdy': {'yldh': 5, 'ujaw': 16},
 'zeiw': {'oyof': 14, 'wnbt': 4, 'sbzh': 16, 'yoke': 13},
 'mjnf': {'iatr': 12, 'owuc': 13, 'izvb': 12, 'axdk': 8},
 'vpvp': {'pqwe': 2, 'gtil': 19, 'qzjs': 15, 'lrgk': 19},
 'flkz': {'ykfj': 20},
 'vawu': {'anbv': 3, 'jfmq': 11, 'gpet': 11, 'bnru': 17},
 'oyzw': {'mclq': 19},
 'raip': {'jrzs': 4},
 'wxwv': {'vzoa': 16, 'inuw': 9},
 'rfcj': {'bmdi': 8, 'zurf': 2, 'sbcj': 17}
 }

print (f38(d))

# %%
#Askhsh 39

def f39(d):
    return { k:v for i in d for k,v in list(d[i].items()) if v % 2 == 1}

d = {
 'fwdy': {'yldh': 5, 'ujaw': 16},
 'zeiw': {'oyof': 14, 'wnbt': 4, 'sbzh': 16, 'yoke': 13},
 'mjnf': {'iatr': 12, 'owuc': 13, 'izvb': 12, 'axdk': 8},
 'vpvp': {'pqwe': 2, 'gtil': 19, 'qzjs': 15, 'lrgk': 19},
 'flkz': {'ykfj': 20},
 'vawu': {'anbv': 3, 'jfmq': 11, 'gpet': 11, 'bnru': 17},
 'oyzw': {'mclq': 19},
 'raip': {'jrzs': 4},
 'wxwv': {'vzoa': 16, 'inuw': 9},
 'rfcj': {'bmdi': 8, 'zurf': 2, 'sbcj': 17}
 }

print(f39(d))

# %%
#Askhsh 40

def f40(opn,clsd):
    list_opn = []
    list_clsd = []
    pool = set(list(range(1,101)))
    for j in clsd:
        list_clsd.append( list(range(j[0], j[1]+1)))
    for i in opn:
        list_opn.append( list(range(i[0], i[1]+1)))
        
    set_clsd = set(sum(list_clsd, []))
    set_opn = set(sum(list_opn, []))
    unwant = set_clsd & set_opn
    res_40 = pool & set_opn
    return res_40 - unwant 

a = [(5, 15), (20, 25), (30, 35)]
b = [(7,22), (25, 34)]

print (f40(a,b)) # Τυπώνει {5,6,23,24,35}