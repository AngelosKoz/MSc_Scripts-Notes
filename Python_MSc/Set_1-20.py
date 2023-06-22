from math import sqrt


#Askhsh 1

def f1(a,b):
    res1 = (a+b)/2
    return res1

print(f1(3,5))

#Askhsh 2


def f2(A,B):
    res2 = f1(A[0], B[0]), f1(A[1],B[1])
    return res2


#Askhsh 3

def f3(A,B):
    if A==B:
        res3 = 'Δεν υπάρχει ευθεία'
    elif A[0]==B[0]:
        res3 = 'Κάθετη'
    else:
        res3 = (B[1]-A[1])/(B[0]-A[0])
    return res3

A = (3,7)
B = (2,5)
print (f3(A,B))

A = (3,7)
B = (3,10)
print (f3(A,B))

A = (3,7)
B = (3,7)
print (f3(A,B))


#Askhsh 4

def f4(A,B):
    res3 = f3(A,B)
    if res3 == 'Δεν υπάρχει ευθεία':
        res4 = res3
    elif A[1]==B[1]:
        res4 = 'Οριζόντια'
    else:
        res4 = -1/res3
    return res4
    


A = (3,7)
B = (2,5)
print (f4(A,B))

A = (3,7)
B = (0,7)
print (f4(A,B))

A = (3,7)
B = (3,7)
print (f4(A, B))


#Askhsh 5

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

c = (0, 0)
r = 1
a = 1
b = 0

point_1, point_2 = f5(a,b,r,c)
print(point_1)
print(point_2)


c = (0, 1)
r = 1
a = 1
b = 0

print(f5(a,b,r,c))


#Askhsh 6 

def f6(a,b,P,d):
    K = P[0]
    L = P[1]
    if L != K*a + b:
        res6 = 'Λάθος'
    else:
        res6 = f5(a,b,d,P)
    return res6


a=1
b=0
P=(0,0)
d=1

point_1, point_2 = f6(a,b,P,d)
print(point_1)
print(point_2)

a=1
b=0
P=(0,1)
d=1

print(f6(a,b,P,d))


#Askhsh 7

def f7(A,B):
    distance = sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
    a = (B[1] - A[1])/(B[0] - A[0])
    b = A[1] - a*A[0]
    if A==B:
        res7 = 'Λάθος'
    elif A[0] == B[0]:
        res7 = (A[0], -distance)
    else:

        res7 = f5(a,b,distance,A)
    return res7

A = (0,0)
B = (1,0)
posf7, negf7 = f7(A,B)
print(negf7)


#Askhsh 8

def f8(A,B,d):
    Γ = f2(A,B)
    λ = f4(A,B)
    b = Γ[1] - λ*Γ[0]
    if A==B:
        res8 = 'Λαθος'
    else:
        res8 = f5(λ,b)
    return res8

A = (-0.7071067811865476, 0.7071067811865476)
B = (0.7071067811865476, -0.7071067811865476)
d = 1

pointK, pointL = f8(A,B,d)
print(pointK)
print(pointL)


#Askhsh 9

def f9(a,b,tolerance = 0.0001):
    return abs(a-b) < tolerance

print(f9(1, 0))
print(f9(1.2345678, 1.234587))


#Askhsh 10

def f10(x):
    x1 = x.split('/')
    res10 = "/".join((x1[1],x1[0],x1[2]))
    return res10

print (f10('12/01/2022'))


from math import sqrt, cos, sin, pi


#Askhsh 11

def f11a(x1):
    
    def f11b(x2):
        if  x2 >= 10 and x2 <= 20:
            return x2
        
    def f11c(x3):
        return x3*2
    
    return sum(list(map(f11c, filter(f11b,x1))))

l = [3,15,7,8,12,20,3]
print (f11a(l)) 


#Askhsh 12

def f12(y):
    return sum([
        x*2 for x in y if x>= 10 and x<=20
    ])

l = [3,15,7,8,12,20,3]
print (f12(l))


#Askhsh 13

def f13(a,b):
    return max([
        sqrt((b1[0] - a1[0])**2 + (b1[1] - a1[1])**2) for a1 in a for b1 in b
    ])

a = [(9.68, 9.55), (4.99, 8.97), (8.67, 6.28), (5.98, 8.99), (7.54, 0.94), (0.04, 7.76), (5.47, 7.83), (8.78, 0.17), (1.4, 4.45), (6.74, 2.76)]
b = [(8.02, 6.82), (2.08, 4.8), (1.85, 5.42), (9.92, 8.86), (0.84, 7.62)]

print (f13(a,b))


#Askhsh 14

#https://www.quora.com/How-do-I-move-a-circle-on-a-graph , xrhsimopoihsa auto elpizw na einai swsto

def f14a(K,r,phi):
    x = cos(phi/360 * 2*pi) *r + K[0]  
    y = sin(phi/360 * 2*pi) *r + K[1]
    return x,y
    

def f14b(K,r,n):
    return [f14a(K,r, (i*90)) for i in range(n)]
    
#Apotelesmata f14a

K = (0,0)
r = 1
phi = 45

print (f14a(K,r,phi))


K = (0,0)
r = 1
phi = 90

print(f14a(K,r,phi))

#Apotelesmata f14b

K = (0,0)
r = 1

print (f14b(K, r, n=3))

print (f14b(K, r, n=4))

print (f14b(K, r, n=5))


#Askhsh 15
#S = Skip, V = Reverse, D = Draw Two, W = Wild, F = Four (eixa thema me ta symbola)

def f15():
    first = list('RYGB')
    second = list('123456789SVD0')
    return list(
        [i+j for i in first for j in second] + [i+j for i in first for j in second[0:12]] + 4*[('W'),('F')]
        )

print(f15())
print(len(f15()))


#Askhsh 16

def f16():
    orektika = ['Σαγανάκι', 'Τυροκαυτερή', 'Πατάτες_Τηγανητές']
    salates = ['Σίζαρ', 'Κόλσλο']
    kyria = ['Πανσέτα', 'Χοχλιοί', 'Αμελέτητα']
    return list([
        (o,s,k) for o in orektika for s in salates for k in kyria
        if (o,k) != ('Σαγανάκι', 'Αμελέτητα')
        if (s,k) != ('Σίζαρ', 'Χοχλιοί')
        if (s,k) != ('Σίζαρ','Αμελέτητα')
    ])

print(f16())


#Askhsh 17

def f17(y):
    return ''.join([
        x[0] for x in y.split()
    ])

print (f17('Alexandros Kanterakis'))
print (f17('Μεταπτυχιακό Βιοπληροφορικής Ηρακλείου Κρήτης'))
print (f17('Περνάω Όμορφα Λύνοντας Υστερικές Και Αρκετά Λετζίτ Ασκήσεις'))


#Askhsh 18

def f18():
    
    return [
        a+b+c+d+e for a in 'ACGT' for b in 'ACGT' for c in 'ACGT' for d in 'ACGT' for e in 'ACGT'
        if a == e and b == d
    ]

print(f18())


#Askhsh 19

def f19(n):
    return [
        [x for x in range(1,y+1)] for y in range(1,n+1)
    ]

print(f19(6))
print(f19(14))


#Askhsh 20

def f20(α,β):
    return [x for x in α for y in β if x in range(y[0], y[1])]

α = [3, 10, 17]
β = [(4,7), (0,1), (13,18), (21,25)]
print (f20(α,β))

α = [3, 10, 17]
β = [(4,7), (0,1), (13,18), (16,25)]

print (f20(α,β))

α = [17, 3, 10, 17]
β = [(4,7), (0,1), (13,18), (16,25)]

print (f20(α,β))
