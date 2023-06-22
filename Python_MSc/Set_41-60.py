
# %% Askhsh 41


def f41(x):                    # Get Initial List
  k = ["".join(sorted(i)) for i in x] # Sort Individual Strings in List l
  l_k = [(a,b) for (a,b) in zip(x,k)]   # Create a map for l and k as it preserves position
  sort_with_k = sorted(l_k, key = lambda x : x[1]) # Sort the map with values at index [1] as they are from k
  #As we already sorted the map with k, the sorted maps first values will contain the elements of rearranged_list
  n = [m[0] for m in sort_with_k]
  return n #Result returned.


l = ['cbc', 'adb', 'dab', 'acb', 'bbc', 'aca', 'bbb', 'aab', 'cad', 'bba']
# Τυπώνει: ['aab', 'aca', 'bba', 'acb', 'adb', 'dab', 'cad', 'bbb', 'bbc', 'cbc'] 

print(f41(l))


# %% #Askhsh 41_Tropos_2
def f41_2(x):
    k = ["".join(sorted(el)) for el in x]
    o = sorted(k)
    index_map = {el:i for i,el in enumerate(o)}
    n = [0]*len(x)
    for i,el in enumerate(k):
        n[index_map[el]] = x[i]
    return n

l = ['cbc', 'adb', 'dab', 'acb', 'bbc', 'aca', 'bbb', 'aab', 'cad', 'bba']

print(f41_2(l))


# %%
#Askhsh 42

def f42(x):
    n = 0
    list_42 = []
    for i in range(len(x)-2):
        diff_42 = abs(x[i] - x[i+1]) + abs(x[i+1] - x[i+2])
        list_42.append(diff_42)
    lowest_diff = min(list_42)
    index_42 = list_42.index(lowest_diff) + 1
    return index_42
    
    
l = [78, 14, 11, 45, 87, 43, 58, 52, 49, 63, 50, 31, 72, 34, 32, 10, 17, 67, 7, 53]
print(f42(l))


# %% 
# Askhsh 42 Allos tropos
def f42(x):
    list_42 = []
    for i in range(len(x)-2):
        diff_42 = abs(x[i] - x[i+1]) + abs(x[i+1] - x[i+2])
        tup_42 = (diff_42, i+1)
        list_42.append(tup_42)
    return min(list_42, key = lambda x : x[0])[1]
    
l = [78, 14, 11, 45, 87, 43, 58, 52, 49, 63, 50, 31, 72, 34, 32, 10, 17, 67, 7, 53]
print(f42(l))
# %% 
# Askhsh 42 Allos tropos tou allou tropou -- Pairnei upopshn thn periptwsh pou bgalw parapanw apo ena min
def f42(x):
    list_42 = []
    check_dups_42 = []
    tup_42 = ()
    for i in range(len(x)-2):
        diff_42 = abs(x[i] - x[i+1]) + abs(x[i+1] - x[i+2])
        check_dups_42.append(diff_42)
        tup_42 = (diff_42, i+1)
        list_42.append(tup_42)
    count_42 = check_dups_42.count(min(check_dups_42))
    if count_42 > 1:
        output_42 = [j for k, j in list_42 if k == min(check_dups_42)]
    else:
        output_42 = min(list_42, key = lambda x : x[0])[1]
        
    
    return output_42


l = [78, 14, 11, 45, 87, 43, 58, 52, 49, 63, 50, 31, 72, 34, 32, 10, 17, 67, 7, 53]
l_dup = [78, 14, 11, 78, 14, 11,] #Dinw auto san paradeiga ti ennow gia thn periptwsh auth
print(f42(l))
print(f42(l_dup))
# %%
# Askhsh 43

def int_to_roman(num):
    m = ["", "M", "MM", "MMM"]
    c = ["", "C", "CC", "CCC", "CD", "D",
         "DC", "DCC", "DCCC", "CM"]
    x = ["", "X", "XX", "XXX", "XL", "L",
         "LX", "LXX", "LXXX", "XC"]
    i = ["", "I", "II", "III", "IV", "V",
         "VI", "VII", "VIII", "IX"]
  
    # Converting to roman
    thousands = m[num // 1000]
    hundreds = c[(num % 1000) // 100]
    tens = x[(num % 100) // 10]
    ones = i[num % 10]
  
    ans = (thousands + hundreds +
           tens + ones)
  
    return ans


def f43():
    len_str = len(int_to_roman(1))
    for i in range(2, 4000):
        roman = int_to_roman(i)
        if len(roman) > len_str:
            list_43 = []
            len_str = len(roman)
            tup_43 = (roman,i)
            list_43.append(tup_43)
            continue
        if len(roman) == len_str: #koitaw th periptwsh opou exw duo arithmous me to idio len()
            tup_43 = (roman,i)
            list_43.append(tup_43)
            
    return list_43


print (f43())


# %% 
# Askhsh 44
import itertools, operator
from math import sqrt


def f44(x):
    list_44 = []
    distance_of_par = []
    paral = sorted(x, key = lambda x : x[0])
    for key,group in itertools.groupby(paral,operator.itemgetter(0)): #This uses index
        list_44.append(list(group))
    for i in list_44: # Each i here is a list with only the parallels
        a = i[0][0] #Since i is a list of the parallels, a is consistent
        for j in (range(len(i) -1)): #Take one element at a time
            for k in range(len(i) - 1): #And compare it to each consecutive element
                if j != k+1: #Compensates the case where it finds distance with itself a.k.a 0
                    dist = ( abs(i[j][1] - i[k+1][1]) ) / (sqrt( (a**2) + 1) )
                    distance_of_par.append(dist)
    return min(distance_of_par)

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

print (f44(l))
# %% Askhsh 45 

def f45(x):
    points = [1,2,3,4,5,6,7,8,10,12]
    points_rev = list(reversed(points))
    lets_see_who_wins = []
    score_sorted = []
    thresh_check = 1
    one_more_list_bcs_i_cant_do_comperhensions = []
    for k,v in x.items():
        score_sorted.append(sorted(v))
    for i in points_rev:
        for j in score_sorted:
            if i in j:
                lets_see_who_wins.append(j)
        
        if lets_see_who_wins:
            score_sorted = lets_see_who_wins
            lets_see_who_wins = []
        if not lets_see_who_wins:
            pass
        if len(score_sorted) == 1:
            who_won = [ k for k,v in x.items() if sorted(v) == score_sorted[0] ]
            break
        for c in score_sorted:
            thresh = (c.count(i))
            if thresh > thresh_check:
                thresh_check = thresh
                one_more_list_bcs_i_cant_do_comperhensions = []
                one_more_list_bcs_i_cant_do_comperhensions.append(c)
                
            elif thresh == thresh_check:
                one_more_list_bcs_i_cant_do_comperhensions.append(c)
            else:
                pass
        thresh_check = 1   
        if len(one_more_list_bcs_i_cant_do_comperhensions) == 1:
            who_won = [ k for k,v in x.items() if sorted(v) == one_more_list_bcs_i_cant_do_comperhensions[0] ]
            break
        score_sorted = one_more_list_bcs_i_cant_do_comperhensions #this is extra in case of bug
        one_more_list_bcs_i_cant_do_comperhensions = []
    print(one_more_list_bcs_i_cant_do_comperhensions) #This is the countrys sorted result
    return who_won


score = {
    'Country_1': [2, 4, 1, 7, 10, 7, 12, 7, 2, 10, 12, 10, 6, 7, 6, 5, 5, 5, 1, 7, 12, 10, 10, 12, 12, 3, 3, 3, 3, 2, 12, 7, 2, 8, 8, 6, 12, 7, 4, 2, 8, 7, 1, 6, 4], 
     'Country_2': [8, 1, 8, 4, 7, 8, 3, 5, 7, 4, 1, 8, 10, 12, 2, 5, 1, 8, 8, 2, 4, 12, 2, 3, 5, 5, 10, 4, 12, 2, 4, 8, 8, 10, 1, 4, 7, 4, 3, 12, 1, 4, 10, 4, 10, 6, 4, 5, 8, 6], 
     'Country_3': [10, 6, 7, 12, 2, 5, 2, 5, 8, 3, 8, 10, 7, 2, 1, 4, 1, 4, 1, 2, 4, 4, 12, 4, 2, 12, 1, 10, 3, 5, 12, 1, 2, 6, 8, 8, 8, 1, 3, 10, 12, 12, 5, 2, 3, 2, 5, 3, 4, 7, 8, 5, 2, 4], 
     'Country_4': [2, 5, 1, 7, 2, 12, 2, 6, 3, 10, 4, 3, 8, 6, 1, 6, 2, 7, 10, 7, 7, 5, 3, 7, 8, 12, 12, 5, 2, 8, 7, 10, 6, 7, 4, 12, 1, 8, 8, 10, 4, 6, 6, 12, 3, 1, 4, 8], 
     'Country_5': [2, 2, 3, 3, 6, 5, 10, 10, 8, 10, 6, 6, 8, 4, 7, 10, 3, 5, 2, 7, 8, 5, 12, 4, 3, 10, 12, 2, 6, 3, 6, 4, 6, 12, 1, 2, 10, 6, 2, 7, 6, 2, 2, 2, 12, 7, 12, 5, 3, 1], 
     'Country_6': [5, 8, 4, 1, 10, 5, 5, 7, 5, 4, 10, 4, 10, 7, 3, 3, 5, 5, 4, 6, 3, 6, 3, 4, 6, 3, 3, 4, 10, 1, 1, 10, 6, 3, 1, 1, 6, 3, 6, 2, 8, 6, 5, 12, 10, 8, 1, 3, 4, 8, 7, 5, 8, 1, 5, 6], 
     'Country_7': [12, 5, 10, 12, 3, 4, 8, 3, 2, 1, 1, 3, 3, 12, 7, 4, 7, 1, 7, 12, 5, 10, 1, 5, 5, 10, 6, 7, 2, 12, 6, 10, 2, 3, 5, 1, 12, 8, 3, 6, 1, 7, 1, 10, 8, 12, 6, 7, 2], 
     'Country_8': [5, 8, 4, 3, 7, 7, 4, 4, 1, 4, 10, 8, 6, 12, 12, 3, 1, 7, 8, 5, 4, 12, 10, 7, 4, 3, 1, 2, 2, 7, 7, 1, 5, 1, 7, 3, 3, 12, 2, 1, 3, 8, 7, 10, 6, 7, 1, 10, 10, 4, 4, 6, 1], 
     'Country_9': [4, 10, 6, 8, 8, 2, 8, 12, 5, 7, 1, 1, 7, 12, 2, 12, 5, 3, 6, 5, 12, 3, 10, 12, 3, 3, 7, 10, 1, 6, 6, 6, 3, 6, 4, 5, 1, 7, 12, 1, 5, 10, 6, 10, 5, 12], 
     'Country_10': [6, 4, 10, 1, 1, 2, 10, 12, 5, 8, 6, 10, 2, 8, 5, 8, 5, 6, 8, 6, 2, 4, 5, 3, 2, 8, 12, 8, 2, 12, 12, 6, 8, 5, 3, 2, 8, 2, 10, 7, 6, 4, 4, 2, 10, 5, 7, 1, 7],
}
print (f45(score))

score_old = {
 'country_1': [4, 8, 10, 2, 10, 10, 4, 1, 1,],
 'country_2': [7, 2, 1, 7, 7, 5, 3, 8, 4, 2, 4],
 'country_3': [7, 12, 8, 10, 7, 4, 1, 1,10],
 'country_4': [3, 2, 6, 6, 2, 1, 12, 2, 1, 7, 2, 1, 1, 4],
 'country_5': [3, 3, 1, 1, 3, 7, 10, 6, 2, 12, 2],
 'country_6': [6, 1, 8, 8, 6, 2, 3, 6, 1, 2, 1, 2, 3, 1],
 'country_7': [10, 5, 7, 12, 8, 1, 3, 1, 1, 2],
 'country_8': [8, 7, 3, 6, 8, 12, 5, 1,],
 'country_9': [2, 3, 6, 12, 6, 4, 6, 4, 4, 3],
 'country_10': [3, 5, 5, 3, 4, 5, 4, 5, 12, 4],
}
print (f45(score_old))

# %%
# Askhsh 45.5 ena arxiko draft...
import random as r

def random_point():
    one_country = []
    all_countries = []
    step = 0
    all_points = [y for x in range(50) for y in [1,2,3,4,5,6,7,8,10,12]]
    print(sum(all_points))
    random.shuffle(all_points)
    
    for i in all_points:
        if sum(one_country) == 290:
            all_countries.append(one_country)
            one_country = []
        elif sum(one_country) < 290:
            one_country.append(i)
        else:
            all_points.append(i)
            one_country = one_country[0:-1]
    for j in all_countries:
        print(sum(j)) # Fainetai pws einai 290 to athroisma omws ligoteres listew
    return len(all_countries) # Kathe fora diaferei, apo oti katalabainw apla ta noumera pou vazw sto telos den mporoun na ahtrisoun 50

print(random_point())


# %%
# askhsh 46
def fmedian(g):
    a = (int((sum(g[0:3]))/3), g[3])
    return list(a)

def f46(x):

    periferies = sorted(list(set([  v[3] for k,v in x.items() ] )))
    search_cities = { k:fmedian(v) for k,v in x.items() }
    time_to_sort = sorted(search_cities.items(), key = lambda x: (x[1][1], x[1][0]))
    bottom_two = []
    city_two = []
    step = 0
    for i in periferies:
        for j in time_to_sort:
            if i == j[1][1]:
                bottom_two.append(j[0])
                step +=1
            if step == 2:
                bot_tup = (i, bottom_two)
                bottom_two = []
                city_two.append(bot_tup)
                step = 0
                break
    return dict(city_two)


cities = {
 'Athens': [772072, 745514, 664046, 'Attica'],
 'Thessaloniki': [383967, 363987, 315196, 'Central Macedonia'],
 'Patras': [152570, 160400, 167446, 'Western Greece'],
 'Piraeus': [182671, 175697, 163688, 'Attica'],
 'Larissa': [112777, 124394, 144651, 'Thessaly'],
 'Heraklion': [115270, 130914, 140730, 'Crete'],
 'Peristeri': [137288, 137918, 139981, 'Attica'],
 'Kallithea': [194233, 109609, 100641, 'Attica'],
 'Acharnes': [61052, 75329, 99346, 'Attica'],
 'Kalamaria': [80698, 87255, 91279, 'Central Macedonia'],
 'Nikaia': [87597, 93086, 89380, 'Attica'],
 'Glyfada': [63306, 80409, 87305, 'Attica'],
 'Volos': [77192, 82439, 86046, 'Thessaly'],
 'Ilio': [78326, 80859, 84793, 'Attica'],
 'Ilioupoli': [75037, 75904, 78153, 'Attica'],
 'Keratsini': [71982, 76102, 77077, 'Attica'],
 'Evosmos': [28821, 52624, 74686, 'Central Macedonia'],
 'Chalandri': [66285, 71684, 74192, 'Attica'],
 'Nea Smyrni': [69749, 73986, 73076, 'Attica'],
 'Marousi': [64092, 69470, 72333, 'Attica'],
 'Agios Dimitrios': [57574, 65173, 71294, 'Attica'],
 'Zografou': [80492, 76115, 71026, 'Attica'],
 'Egaleo': [78563, 74046, 69946, 'Attica'],
 'Nea Ionia': [27904, 30804, 32661, 'Thessaly'],
 'Ioannina': [56699, 61629, 65574, 'Epirus'],
 'Palaio Faliro': [61371, 64759, 64021, 'Attica'],
 'Korydallos': [63184, 67456, 63445, 'Attica'],
 'Trikala': [45835, 48686, 61653, 'Thessaly'],
 'Vyronas': [58523, 61102, 61308, 'Attica'],
 'Agia Paraskevi': [47463, 56836, 59704, 'Attica'],
 'Galatsi': [57230, 58042, 59345, 'Attica'],
 'Agrinio': [52081, 54523, 59329, 'Western Greece'],
 'Chalcis': [51646, 53584, 59125, 'Central Greece'],
 'Petroupoli': [38278, 48327, 58979, 'Attica'],
 'Serres': [50017, 54266, 58287, 'Central Macedonia'],
 'Alexandroupoli': [37904, 48885, 57812, 'Eastern Macedonia and Thrace'],
 'Xanthi': [37430, 45111, 56122, 'Eastern Macedonia and Thrace'],
 'Katerini': [43613, 50510, 55997, 'Central Macedonia'],
 'Kalamata': [43625, 49154, 54100, 'Peloponnese'],
 'Kavala': [56571, 58663, 54027, 'Eastern Macedonia and Thrace'],
 'Chania': [50077, 53373, 53910, 'Crete'],
 'Lamia': [44084, 46406, 52006, 'Central Greece'],
 'Komotini': [37036, 43326, 50990, 'Eastern Macedonia and Thrace'],
 'Irakleio': [42905, 45926, 49642, 'Attica'],
 'Rhodes': [42400, 52318, 49541, 'South Aegean'],
 'Kifissia': [39166, 43929, 47332, 'Attica'],
 'Stavroupoli': [37596, 41653, 46008, 'Central Macedonia'],
 'Chaidari': [44831, 45227, 45642, 'Attica'],
 'Drama': [37604, 42501, 44823, 'Eastern Macedonia and Thrace'],
 'Veria': [37858, 42794, 43158, 'Central Macedonia'],
 'Alimos': [32024, 38047, 41720, 'Attica'],
 'Kozani': [31553, 35242, 41066, 'Western Macedonia'],
 'Polichni': [27894, 36146, 39332, 'Central Macedonia'],
 'Karditsa': [30067, 32031, 38554, 'Thessaly'],
 'Sykies': [34059, 41726, 37753, 'Central Macedonia'],
 'Ampelokipoi': [40093, 40959, 37381, 'Central Macedonia'],
 'Pylaia': [20785, 22744, 34625, 'Central Macedonia'],
 'Agioi Anargyroi': [30739, 32957, 34168, 'Attica'],
 'Argyroupoli': [31530, 33158, 34097, 'Attica'],
 'Ano Liosia': [21397, 26423, 33565, 'Attica'],
 'Rethymno': [23420, 27868, 32468, 'Crete'],
 'Ptolemaida': [25125, 28679, 32127, 'Western Macedonia'],
 'Tripoli': [22429, 25520, 30866, 'Peloponnese'],
 'Cholargos': [33691, 32166, 30840, 'Attica'],
 'Vrilissia': [16571, 25582, 30741, 'Attica'],
 'Aspropyrgos': [15715, 27741, 30251, 'Attica'],
 'Corinth': [27412, 29787, 30176, 'Peloponnese'],
 'Gerakas': [8512, 13921, 29939, 'Attica'],
 'Metamorfosi': [21052, 26448, 29891, 'Attica'],
 'Giannitsa': [22504, 26296, 29789, 'Central Macedonia'],
 'Voula': [17998, 25532, 28364, 'Attica'],
 'Kamatero': [17410, 22234, 28361, 'Attica'],
 'Mytilene': [23971, 27247, 27871, 'North Aegean'],
 'Neapoli': [30568, 29995, 27084, 'Central Macedonia'],
 'Eleftherio-Kordelio': [16549, 21630, 27067, 'Central Macedonia'],
 'Chios': [22894, 23779, 26850, 'North Aegean'],
 'Agia Varvara': [28706, 30562, 26550, 'Attica'],
 'Kaisariani': [26701, 26323, 26370, 'Attica'],
 'Nea Filadelfeia': [25261, 24112, 25734, 'Attica'],
 'Moschato': [22039, 23153, 25441, 'Attica'],
 'Perama': [24119, 25720, 25389, 'Attica'],
 'Salamina': [22567, 25730, 25370, 'Attica'],
 'Eleusis': [22793, 25863, 24910, 'Attica'],
 'Corfu': [31359, 28185, 24838, 'Ionian Islands'],
 'Pyrgos': [28465, 23274, 24359, 'Western Greece'],
 'Megara': [20403, 23032, 23456, 'Attica'],
 'Kilkis': [12139, 17430, 22914, 'Central Macedonia'],
 'Dafni': [24152, 23674, 22913, 'Attica'],
 'Thebes': [19505, 21211, 22883, 'Central Greece'],
 'Melissia': [13469, 19526, 22741, 'Attica'],
 'Argos': [21901, 24239, 22209, 'Peloponnese'],
 'Arta': [19087, 19435, 21895, 'Epirus'],
 'Artemida': [9485, 17391, 21488, 'Attica'],
 'Livadeia': [18437, 20061, 21379, 'Central Greece'],
 'Pefki': [17987, 19887, 21352, 'Attica'],
 'Oraiokastro': [5458, 11896, 20852, 'Central Macedonia'],
 'Aigio': [22178, 21061, 20422, 'Western Greece'],
 'Kos': [14714, 17890, 19432, 'South Aegean'],
 'Koropi': [12790, 15860, 19164, 'Attica'],
 'Preveza': [13695, 16321, 19042, 'Epirus'],
 'Naousa': [19794, 19870, 18882, 'Central Macedonia'],
 'Orestiada': [12691, 15246, 18426, 'Eastern Macedonia and Thrace'],
 'Peraia': [2949, 13306, 18326, 'Central Macedonia'],
 'Edessa': [17128, 18253, 18229, 'Central Macedonia'],
 'Florina': [12355, 14279, 17686, 'Western Macedonia'],
 'Panorama': [10275, 14552, 17444, 'Central Macedonia'],
 'Nea Erythraia': [12993, 15439, 17379, 'Attica'],
 'Elliniko': [13517, 16740, 17259, 'Attica'],
 'Amaliada': [15232, 18261, 16763, 'Western Greece'],
 'Pallini': [8021, 12552, 16415, 'Attica'],
 'Sparta': [13011, 14817, 16239, 'Peloponnese'],
 'Agios Ioannis Rentis': [14218, 15060, 16050, 'Attica'],
 'Thermi': [5156, 11360, 16004, 'Central Macedonia'],
 'Vari': [8488, 10998, 15855, 'Attica'],
 'Nea Makri': [12120, 13986, 15554, 'Attica'],
 'Tavros': [15456, 14963, 14972, 'Attica'],
 'Alexandreia': [12109, 13229, 14821, 'Central Macedonia'],
 'Menemeni': [12932, 14910, 14746, 'Central Macedonia'],
 'Paiania': [9710, 12855, 14595, 'Attica'],
 'Kalyvia Thorikou': [8488, 12202, 14424, 'Attica'],
 'Nafplio': [11897, 13822, 14203, 'Peloponnese'],
 'Drapetsona': [13094, 12944, 13968, 'Attica'],
 'Efkarpia': [3480, 6598, 13905, 'Central Macedonia'],
 'Papagou': [13974, 13207, 13699, 'Attica'],
 'Nafpaktos': [10854, 12924, 13415, 'Western Greece'],
 'Kastoria': [14775, 14813, 13387, 'Western Macedonia'],
 'Grevena': [9345, 10177, 13137, 'Western Macedonia'],
 'Pefka': [3561, 6434, 13052, 'Central Macedonia'],
 'Nea Alikarnassos': [10683, 11551, 12925, 'Crete'],
 'Missolonghi': [10916, 12225, 12785, 'Western Greece'],
 'Gazi': [1395, 8018, 12606, 'Crete'],
 'Ierapetra': [9541, 11678, 12355, 'Crete'],
 'Kalymnos': [10543, 10149, 12324, 'South Aegean'],
 'Rafina': [7752, 11352, 12168, 'Attica'],
 'Loutraki': [9388, 11383, 11564, 'Peloponnese'],
 'Agios Nikolaos': [8093, 10080, 11421, 'Crete'],
 'Ermoupoli': [13030, 11799, 11407, 'South Aegean'],
 'Ialysos': [7193, 10107, 11331, 'South Aegean'],
 'Mandra': [10012, 10947, 11327, 'Attica'],
 'Tyrnavos': [12028, 11116, 11069, 'Thessaly'],
 'Glyka Nera': [5813, 6623, 11049, 'Attica'],
 'Ymittos': [11671, 11139, 10715, 'Attica'],
 'Neo Psychiko': [12023, 10848, 10137, 'Attica'],
}

print(f46(cities))

# %%
# Askhsh 47

def f47(x):
    
    contained = 0
    tup_con = ()
    lets_count = []
    test = list(enumerate(x))
    for i in test:
        for j in test:
            if i[0] == j[0]:
                print('Please dont compare me with myself\n')
            else:
                if i[1][0] < j[1][0] and i[1][1] > j[1][1]:
                    contained += 1
                    print(f'{j[1]} is contained in {i[1]}')
        if contained:
            tup_con = (contained, i[0])
            lets_count.append(tup_con)
            contained = 0
    print(lets_count)
    return x[max(lets_count, key = lambda x: x[0])[1]] #Apo to max sto metrhma vash tou contained, dwse mou to index toy x poy periexetai sto tuple [1] kathe stoixeioy ths listas
    

l = [(1, 2), (3, 4), (5, 6), (0, 7), (-1, 5.5)]
print(f47(l))

# %% Askhsh 47... alla dinei kai thn periptwsh opoy exw parapanw apo 1 apotelesmata

def f47(x):
    
    contained = 0
    tup_con = ()
    test = list(enumerate(x))
    check_dup_47 = [(0,0)]
    for i in test:
        for j in test:
            if i[0] == j[0]:
                print(f'Please dont compare me a.k.a i={i[0]} with my other self j={j[0]}\n')
            else:
                if i[1][0] < j[1][0] and i[1][1] > j[1][1]:
                    contained += 1
                    print(f'{j[1]} is contained in {i[1]}')
        if contained > check_dup_47[0][0]:
                tup_con = (contained, i[0]) #Shows how many times it is contained , index
                check_dup_47[0] = tup_con
                del(check_dup_47[1:-1])
                contained = 0
        elif contained and contained == check_dup_47[0][0]:
            tup_con = (contained, i[0]) #Shows how many times it is contained , index
            check_dup_47.append(tup_con)
            contained = 0
        else:
            pass
    if len(check_dup_47) == 1:
        output_47 = x[check_dup_47[0][1]]
    else:
        output_47 = [x[check_dup_47[k][1]] for k in range(len(check_dup_47)) ]    
    return output_47

            
l = [(1, 2), (3, 4), (5, 6), (0, 7), (-2, 5.5),(-1,0)] #Ekana ena test wste na dw an paizei
print(f47(l))

# %% 
#Askhsh 48
a = lambda x: lambda y: lambda z: (x+y+z)/3
print (a(2)(4)(9))


# %%
#Askhsh 49

from math import sqrt

def f49():
    lets_sort = []
    sum_digits = 0
    for i in range(1,101):
        division = 1/sqrt(i)
        digits = str(division)[2:7]
        for j in digits:
            sum_digits += int(j)
        
        tup_49 = (i, sum_digits)
        lets_sort.append(tup_49)
        sum_digits = 0
    return sorted(lets_sort, key= lambda x: x[1])
print(f49())
    
# %%
# Askhsh 50
def f50(x,y):
    who_won = {}
    candi_score = []
    for ky, vy in y.items():
        failed = False
        for k,v in vy.items():
            if v >= sports[k]['threshold']: #prints {thresh.. quot..}
                candi_score.append((v*sports[k]['quotient']))
            else:
                print(f'{ky} has failed in {k} with a score of {v}')
                candi_score = []
                failed = True
                break
        if not failed:
            who_won[ky]= sum(candi_score)
            candi_score = []
    winner = max(who_won, key = lambda x: who_won[x])
    return f' {winner} won with {who_won[winner]} points across all sports'
sports = {
 'sport_1': {'threshold': 21, 'quotient': 9.858},
 'sport_2': {'threshold': 23, 'quotient': 3.558},
 'sport_3': {'threshold': 29, 'quotient': 0.143},
 'sport_4': {'threshold': 24, 'quotient': 2.852},
 'sport_5': {'threshold': 28, 'quotient': 3.461},
 'sport_6': {'threshold': 21, 'quotient': 4.071},
 'sport_7': {'threshold': 27, 'quotient': 3.279},
 'sport_8': {'threshold': 23, 'quotient': 8.091},
 'sport_9': {'threshold': 21, 'quotient': 0.199},
 'sport_10': {'threshold': 30, 'quotient': 5.088}
}

results = {
    "condidate_1": {
        "sport_1": 121,
        "sport_2": 98,
        "sport_3": 25,
        "sport_4": 20,
        "sport_5": 28,
        "sport_6": 29,
        "sport_7": 82,
        "sport_8": 80,
        "sport_9": 20,
        "sport_10": 122
    },
    "condidate_2": {
        "sport_1": 103,
        "sport_2": 95,
        "sport_3": 129,
        "sport_4": 68,
        "sport_5": 141,
        "sport_6": 96,
        "sport_7": 30,
        "sport_8": 90,
        "sport_9": 21,
        "sport_10": 100
    },
    "condidate_3": {
        "sport_1": 98,
        "sport_2": 144,
        "sport_3": 63,
        "sport_4": 139,
        "sport_5": 25,
        "sport_6": 81,
        "sport_7": 18,
        "sport_8": 118,
        "sport_9": 146,
        "sport_10": 148
    },
    "condidate_4": {
        "sport_1": 87,
        "sport_2": 81,
        "sport_3": 138,
        "sport_4": 114,
        "sport_5": 42,
        "sport_6": 95,
        "sport_7": 142,
        "sport_8": 140,
        "sport_9": 98,
        "sport_10": 52
    },
    "condidate_5": {
        "sport_1": 48,
        "sport_2": 100,
        "sport_3": 103,
        "sport_4": 66,
        "sport_5": 17,
        "sport_6": 14,
        "sport_7": 27,
        "sport_8": 148,
        "sport_9": 68,
        "sport_10": 19
    },
    "condidate_6": {
        "sport_1": 116,
        "sport_2": 54,
        "sport_3": 89,
        "sport_4": 12,
        "sport_5": 125,
        "sport_6": 148,
        "sport_7": 138,
        "sport_8": 66,
        "sport_9": 27,
        "sport_10": 35
    },
    "condidate_7": {
        "sport_1": 146,
        "sport_2": 49,
        "sport_3": 49,
        "sport_4": 148,
        "sport_5": 19,
        "sport_6": 119,
        "sport_7": 49,
        "sport_8": 14,
        "sport_9": 146,
        "sport_10": 79
    },
    "condidate_8": {
        "sport_1": 127,
        "sport_2": 28,
        "sport_3": 107,
        "sport_4": 55,
        "sport_5": 107,
        "sport_6": 135,
        "sport_7": 123,
        "sport_8": 75,
        "sport_9": 29,
        "sport_10": 137
    },
    "condidate_9": {
        "sport_1": 50,
        "sport_2": 40,
        "sport_3": 107,
        "sport_4": 131,
        "sport_5": 67,
        "sport_6": 74,
        "sport_7": 48,
        "sport_8": 80,
        "sport_9": 116,
        "sport_10": 41
    },
    "condidate_10": {
        "sport_1": 125,
        "sport_2": 98,
        "sport_3": 84,
        "sport_4": 129,
        "sport_5": 97,
        "sport_6": 39,
        "sport_7": 42,
        "sport_8": 68,
        "sport_9": 67,
        "sport_10": 52
    },
    "condidate_11": {
        "sport_1": 39,
        "sport_2": 37,
        "sport_3": 69,
        "sport_4": 81,
        "sport_5": 126,
        "sport_6": 39,
        "sport_7": 31,
        "sport_8": 67,
        "sport_9": 119,
        "sport_10": 45
    },
    "condidate_12": {
        "sport_1": 143,
        "sport_2": 12,
        "sport_3": 118,
        "sport_4": 114,
        "sport_5": 109,
        "sport_6": 108,
        "sport_7": 95,
        "sport_8": 82,
        "sport_9": 137,
        "sport_10": 64
    },
    "condidate_13": {
        "sport_1": 131,
        "sport_2": 85,
        "sport_3": 11,
        "sport_4": 146,
        "sport_5": 41,
        "sport_6": 97,
        "sport_7": 64,
        "sport_8": 145,
        "sport_9": 117,
        "sport_10": 110
    },
    "condidate_14": {
        "sport_1": 103,
        "sport_2": 67,
        "sport_3": 49,
        "sport_4": 100,
        "sport_5": 111,
        "sport_6": 73,
        "sport_7": 107,
        "sport_8": 94,
        "sport_9": 90,
        "sport_10": 73
    },
    "condidate_15": {
        "sport_1": 96,
        "sport_2": 17,
        "sport_3": 22,
        "sport_4": 83,
        "sport_5": 59,
        "sport_6": 99,
        "sport_7": 21,
        "sport_8": 18,
        "sport_9": 52,
        "sport_10": 87
    },
    "condidate_16": {
        "sport_1": 137,
        "sport_2": 105,
        "sport_3": 36,
        "sport_4": 115,
        "sport_5": 38,
        "sport_6": 102,
        "sport_7": 17,
        "sport_8": 56,
        "sport_9": 116,
        "sport_10": 64
    },
    "condidate_17": {
        "sport_1": 134,
        "sport_2": 14,
        "sport_3": 38,
        "sport_4": 50,
        "sport_5": 65,
        "sport_6": 51,
        "sport_7": 47,
        "sport_8": 128,
        "sport_9": 132,
        "sport_10": 85
    },
    "condidate_18": {
        "sport_1": 124,
        "sport_2": 117,
        "sport_3": 55,
        "sport_4": 37,
        "sport_5": 14,
        "sport_6": 79,
        "sport_7": 128,
        "sport_8": 66,
        "sport_9": 35,
        "sport_10": 19
    },
    "condidate_19": {
        "sport_1": 127,
        "sport_2": 54,
        "sport_3": 94,
        "sport_4": 26,
        "sport_5": 109,
        "sport_6": 49,
        "sport_7": 63,
        "sport_8": 123,
        "sport_9": 52,
        "sport_10": 127
    },
    "condidate_20": {
        "sport_1": 10,
        "sport_2": 57,
        "sport_3": 136,
        "sport_4": 37,
        "sport_5": 67,
        "sport_6": 113,
        "sport_7": 69,
        "sport_8": 87,
        "sport_9": 145,
        "sport_10": 52
    },}

print(f50(sports,results))

# %%
#Askhsh 51

import re
def f51(x):
    results_51 = []
    pattern_51 = re.compile('\(\d+\)') #thelei to escape
    with open(x) as alien_file:
        for lines in alien_file:
            result = pattern_51.findall(lines)
            if result:
                for i in range(len(result)):
                    corrected = re.sub(r'[()]', '', result[i])
                    
                    results_51.append(int(corrected))
    return(max(results_51))
print(f51('alien.txt'))

# %%
#Askhsh 51 Tropos 2 


import re
def f51(x):
    results_51 = []
    with open(x) as alien_file:
        for lines in alien_file:
            result = re.findall(r'\(\d+\)', lines)
            if result:
                for i in range(len(result)):
                    test = re.sub(r'[()]', '', result[i])
                    results_51.append(int(test))
    return(max(results_51))
print(f51('alien.txt'))

# %%
# Askhsh 52
def f52(x):
    results_52 = []
    pattern_52 = re.compile('\(\s*\+?\s*\d+[ \d]*\)')
    with open(x) as alien_file:
        for lines in alien_file:
            result = pattern_52.findall(lines)
            if result:
                for i in range(len(result)):
                    corrected = re.sub(r'[+()\s]', '', result[i])
                    results_52.append(int(corrected))
    return(max(results_52))
print(f52('alien.txt'))


# %%
# Askhsh 53
import re
def f53(x):
    results_53 = []
    pattern_53 = re.compile('^\^.*\$$') #thelei to escape
    with open(x) as alien_file:
        for lines in alien_file:
            result = pattern_53.findall(lines)
            for i in range(len(result)):
                results_53.append(i)

    return(len(results_53))
print(f53('alien.txt'))

#%% 
# Askhsh 54
import re
def f54(x):
    results_54 = []
    pattern_54 = re.compile('\s+[0-9]*[02468]+\s+[0-9]*[02468]+\s+') #thelei to escape
    sec_pat = 0
    with open(x) as alien_file:
        for lines in alien_file:
            result = pattern_54.findall(lines)
            if result:
                for i in range(len(result)):
                    sum_54 = re.findall('\d+', result[i])
                    for j in sum_54:
                        sec_pat += int(j)
                    results_54.append(sec_pat)
                    sec_pat = 0
                    
    return max(results_54)
print(f54('alien.txt'))

# %% 
# Askhsh 55 
#Me diafora regex opws 
#1) \S[a-zA-Z]{1}\d{1}[a-zA-Z]{1}\d{1}[a-zA-Z]{1}\d{1}[a-zA-Z]{1}\d{1}\S
#2) \W[a-zA-Z]{1}\d{1}[a-zA-Z]{1}\d{1}[a-zA-Z]{1}\d{1}[a-zA-Z]{1}\d{1}\W
#3) [a-zA-Z]{4}\d{4}
# Parathrhsa pws to L0E4L9v7 einai to mono apotelesma pou mprosta kai pisw tou periexei whitespace enw ta alla paroysiazan xarakthres opws @ { ! klp
import re
def f55(): 
    pattern_55 = re.compile('\ {1}[a-zA-Z]{1}\d{1}[a-zA-Z]{1}\d{1}[a-zA-Z]{1}\d{1}[a-zA-Z]{1}\d{1}\ {1}')
    with open('alien.txt') as alien_file:
        for lines in alien_file:
            result = pattern_55.findall(lines)
            if result:
                corrected = re.sub(r'[ ]', '', result[0]) # Kathws dieukrinizetai pws yparxei 1 fora den xreiazetai na psaksw an exei perissotero apo 1 luseis
    return corrected
                  
print(f55())  

# %% 
# Askhsh 56
def f56():
    pattern_56_1 = re.compile('69\d{8}') #Dinei 2 apotelesmata poy shmainei mallon ena apo auta exei arihmo prin h meta
    pattern_56_2 = re.compile('.69\d{8}.') #Epivevaiwnei pws sto deutero apotelesma exw arithmo meta
    pattern_56 = re.compile('\D69\d{8}\D') #Fainetai pws auto dinei th swsth akolouthia
    
    with open('alien.txt') as alien_file:
        for lines in alien_file:
            result = pattern_56.findall(lines)
            if result:
                corrected = re.sub(r'[\D]', '', result[0])
    return corrected
print(f56())

# %% Askhsh 57
import re
def f57():
    regex_list = ['\ +[a-zA-Z_]+[\w_]*\ *\=\ *\d+\ +',
                  '^[a-zA-Z_]+[\w_]*\ *\=\ *\d+\ +',
                  '\ +[a-zA-Z_]+[\w_]*\ *\=\ *\d+$',]
    result_57 = []
    with open('alien.txt') as alien_file:
        for lines in alien_file:
            for pat in regex_list:
                result = re.findall(pat,lines)
                if result:
                    for i in range(len(result)):
                        split_result = result[i].split("=")
                        corrected = re.sub(r' ', '', split_result[1])
                        corrected = int(corrected)
                        result_57.append(corrected)
    return(sum(result_57))

print(f57())

# %%
# Askhsh 58 Den kserw oute egw ti ekana alla paizei!

def ask_58_create_file():


    with open('ask_58_1.txt', 'w') as f:

        f.write(f'A=5,6,7,8,9|B=8,6,7|C=5,4,5,3,2\n')
        f.write(f'A=1,6,7,6,5|B=6,8,9|C=7,8,8,8,8\n')
        f.write(f'A=7,6,9,2,1|B=1,2,3|C=7,8,1,1,2\n')

ask_58_create_file()


import re, itertools, operator
def f58(x):
    letter_set = []
    result_58 = []
    another_list = []
    third_list = []
    with open(x) as ask_file:
        open('ask_58_2_for_kantale.txt', 'w').close() #PROSOXH AN EXETE TETOIO ARXEIO ADEIAZEI CONTENTS
        for let in ask_file:
            letters_58 = re.findall('[A-Za-z]+', let)
            for i in letters_58:
                letter_set.append(i)
        got_letters = (sorted(list(set(letter_set)))) #finding every letter used as variable
        ask_file.seek(0)
        for l in got_letters:
            for line in ask_file:
                numbers_58 = re.findall(f'{l}=[\d\,]*', line) #me to f' pernaw to gramma san variable
                numbers_58_split_1 = re.sub('[A-Za-z]*=', '', numbers_58[0])
                numbers_58_split_2 = numbers_58_split_1.split(',')
                num_enu = list(enumerate(numbers_58_split_2))
                for j in num_enu:
                    result_58.append(j)
            result_58 = sorted(result_58, key=operator.itemgetter(0))
            for key,group in itertools.groupby(result_58,operator.itemgetter(0)):
                temp = (list(group))
                another_list.append(temp)
            ask_file.seek(0)
            result_58 = []
            for dont_know_anymore in another_list:
                for really_i_dont_know in range(len(dont_know_anymore)):
                    int_58 = int(dont_know_anymore[really_i_dont_know][1])
                    result_58.append(int_58)
                temp_tup = tuple(result_58)
                third_list.append(temp_tup)
                result_58 = []
            another_list = [] 
            result_58 = []
            test = str(third_list).strip('[]')
            test = re.sub(' ', '', test)
            finally_its_over = re.sub('\)\,\(' , '), (', test)
            third_list = []
            with open('ask_58_2_for_kantale.txt', 'a') as f:
                f.write(f'{l}={finally_its_over}\n')
f58('ask_58_1.txt')

with open('ask_58_2_for_kantale.txt') as g:
    for lin in g:
        print (lin.strip())
        
        
# %% 
# Askhsh 59 - Thewroume pws sto arxeio kathe grammh exei apo ena gramma, alliws tha eftiaxna to arxeio sth morfh pou thelw wste na treksei to parakatw

def ask_59_create_file():


    with open('ask_59_1.txt', 'w') as f:

        f.write(f'A=(5,1,7), (6,6,6), (7,7,9), (8,6,2), (9,5,1)\n')
        f.write(f'B=(8,6,1), (6,8,2), (7,9,3)\n')
        f.write(f'C=(5,7,7), (4,8,8), (5,8,1), (3,8,1), (2,8,2)\n')
ask_59_create_file()


import re, itertools, operator

def f59(x):
    letter_set = []
    result_59 = []
    another_list = []
    third_list = []
    enu_list = []
    with open(x) as ask_file:
        open('ask_59_2_for_kantale.txt', 'w').close() #PROSOXH AN EXETE TETOIO ARXEIO ADEIAZEI CONTENTS. To exw valei logw tou append sto open with opote se periptwsh pou ksanatrseksei to script then tha diplasiasei ta contents
        for let in ask_file:
            letters_59 = re.findall('[A-Za-z]+', let)
            for i in letters_59:
                letter_set.append(i)
        got_letters = (sorted(list(set(letter_set))))
        ask_file.seek(0)
        for l in got_letters:
            for line in ask_file:
                numbers_59 = re.findall(f'{l}=.*', line)
                if not numbers_59:
                    continue
                numbers_59_split_1 = re.sub(r'[A-Za-z]*=', '', numbers_59[0])
                numbers_59_split_2 = numbers_59_split_1.split(', ')
                res = [eval(ele) for ele in numbers_59_split_2] #https://www.geeksforgeeks.org/python-convert-string-tuples-to-list-tuples/
                for enmt in res:
                    enu_59 = list(enumerate(enmt))
                    for enu_ap in enu_59:
                        result_59.append(enu_ap)
                result_59 = sorted(result_59, key=operator.itemgetter(0))
                for key,group in itertools.groupby(result_59,operator.itemgetter(0)):
                    temp = (list(group))
                    another_list.append(temp)
            ask_file.seek(0)
            result_59 = []
            for dont_know_anymore in another_list:
                for really_i_dont_know in range(len(dont_know_anymore)):
                    int_59 = int(dont_know_anymore[really_i_dont_know][1])
                    result_59.append(int_59)
                temp_tup = tuple(result_59)
                third_list.append(temp_tup)
                result_59 = []
            another_list = [] 
            result_59 = []
            test = list(enumerate(third_list))
            for whynot in test:
                enu_list.append(whynot)
            third_list = []
        enu_list = sorted(enu_list, key = lambda x: x[0])
        for enukey,enugroup in itertools.groupby(enu_list,operator.itemgetter(0)):
                temp = (list(enugroup))
                another_list.append(temp)
        enu_list = []
        for enu_l in another_list:
            for enu_ind in range(len(enu_l)):
                str_59 = str(enu_l[enu_ind][1])
                result_59.append(str_59)
            enu_list.append(result_59)
            result_59 = []
        for enu_2 in enu_list:
            test = list(map('='.join, zip(got_letters, enu_2)))
            test = str(test).strip('[]')
            test = re.sub(r'\'', '', test)
            test = re.sub(r' ', '', test)
            test = re.sub('\)\,' , ')|', test)
            finally_its_over_again = re.sub(r'[()]' , '', test)
            with open('ask_59_2_for_kantale.txt', 'a') as nightmare:
                nightmare.write(f'{finally_its_over_again}\n')
f59('ask_59_1.txt')


with open('ask_59_2_for_kantale.txt') as reading_nightmare:
    for lin in reading_nightmare:
        print (lin.strip())
        
        
# %% 
# Askhsh 60

import re
def f60(x):
    result_60 = 0
    regex_list = ['^[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}$', #to kserw enai pole mpempe h lush
                  '^[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}$',
                  '^[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}$',
                  '^[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}$',
                  '^[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}$',
                  '^[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}$',
                  '^[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}$',
                  '^[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}$',
                  '^[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}[βγδζθκλμνξπρστφχψ]{1}[αεηιουω]{1}$',
                  ]
    with open(x) as greek_words:
        open('easy_words_for_children.txt', 'w').close() #PROSOXH AN EXETE TETOIO ARXEIO ADEIAZEI CONTENTS
        for words in greek_words:
            for w_pat in regex_list:
                result = re.findall(w_pat,words)
                if result:
                    #result_60 += 1 #gia metro sugrishs
                    corrected = str(result).strip('[]\'')
                    with open('easy_words_for_children.txt', 'a') as g:
                        g.write(f'{corrected}\n')
                    result = False
        #return result_60 #extra apla gia na tsekarw to arxeio meta
f60('words_greek_normalized.txt')

with open('easy_words_for_children.txt') as f:
    c = 0
    for l in f:
        c += 1
    print (f'Σύνολο: {c} λέξεις με συλλαβές')
# %%
