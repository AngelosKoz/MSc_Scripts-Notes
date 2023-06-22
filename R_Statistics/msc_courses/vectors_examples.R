#Using the : operator
#This is useful if we want to generate vector of consecutive integers

a <- 1:10
b <- 50:70
d <- 70:50
a
##  [1]  1  2  3  4  5  6  7  8  9 10
b
##  [1] 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70
d
##  [1] 70 69 68 67 66 65 64 63 62 61 60 59 58 57 56 55 54 53 52 51 50
#Notice, that in the last case of vector d, the order of the elements is from the greatest to the smallest.

#Using the seq function

#The seq function allows to define the first and the last elements of a vector as well as the step OR the number of elelements. In the later case, the elements are distributed at equal distances. For example:

## a vector for which the first element is 1, the step is 2 and the last is NOT GREATER than 20. Thus here: 
a <- seq(from=1, to=20, by=2)
a
##  [1]  1  3  5  7  9 11 13 15 17 19
## a vector from 1 to 20 with 30 elements in total
b <- seq(from = 1, to=20, length.out = 30)
b
##  [1]  1.000000  1.655172  2.310345  2.965517  3.620690  4.275862  4.931034
##  [8]  5.586207  6.241379  6.896552  7.551724  8.206897  8.862069  9.517241
## [15] 10.172414 10.827586 11.482759 12.137931 12.793103 13.448276 14.103448
## [22] 14.758621 15.413793 16.068966 16.724138 17.379310 18.034483 18.689655
## [29] 19.344828 20.000000



#Using the rep function
#The rep
# function is very useful because it allows to generate vectors with repetitive elements. For example, if we want to generate a vector of 10 1s then:

## a vector of ten 1s
a <- rep(1, 10)
a
##  [1] 1 1 1 1 1 1 1 1 1 1
## a vector of ten "R"
a <- rep("R", 10)
a
##  [1] "R" "R" "R" "R" "R" "R" "R" "R" "R" "R"
## a vector of ten "WORDS"
a <- rep("WORDS", 10)
a
##  [1] "WORDS" "WORDS" "WORDS" "WORDS" "WORDS" "WORDS" "WORDS" "WORDS"
##  [9] "WORDS" "WORDS"


#Using the c() function
#The c
# function (from the word concatenate), allows to concatenate elements to a vector (or vectors to another vector). For example:

## a simple example
a <- c(1,2,3,4,5,6,6)

## a new vector from two other vectors
b <- c(c(1,2,3), c(34,5,6))
It is of course possible to combine these methods

a <- 1:10
b <- rep(1:5, 4)
d <- c(a,b)
d
##  [1]  1  2  3  4  5  6  7  8  9 10  1  2  3  4  5  1  2  3  4  5  1  2  3
## [24]  4  5  1  2  3  4  5


#Vectors of random numbers
#Since R is a statistical language, it is very useful to be able to draw random numbers from distributions

#a vector of random elements from a uniform
## 100 random values from a uniform distribution with min:0 and max:10
a <- runif(100, 1, 10)
a
##   [1] 7.464879 6.886708 6.203531 3.992212 7.263983 7.207444 7.445366
##   [8] 2.848153 4.877889 8.626188 7.426321 2.534370 6.272095 1.644620
##  [15] 2.837092 6.429066 2.983051 9.633578 1.146730 6.964132 9.067826
##  [22] 2.907615 8.708384 3.544351 9.764700 8.006667 7.439672 1.069004
##  [29] 2.768810 6.392695 9.598843 5.177910 3.929282 6.404419 6.756739
##  [36] 8.307353 1.026679 5.908097 3.315683 2.827687 5.155643 8.252820
##  [43] 5.564817 1.220608 2.775947 8.037676 8.977219 5.765221 3.561910
##  [50] 5.497651 8.491035 9.157577 2.752755 8.124427 6.520487 6.490044
##  [57] 6.429695 4.286225 7.691840 9.430905 4.594378 9.152012 2.886251
##  [64] 4.817335 6.987277 9.316862 1.773513 7.130113 3.018482 1.869190
##  [71] 6.702309 8.200955 4.304943 9.457422 2.634071 6.523368 1.174948
##  [78] 1.677302 9.630140 4.122888 2.611649 6.613924 5.990192 1.748773
##  [85] 6.205209 2.470271 4.265555 5.350173 6.888210 2.809349 9.967937
##  [92] 7.416900 6.046988 6.419863 1.084900 5.480109 7.613014 2.712113
##  [99] 5.286768 4.626499
a vector of random values from a Gaussian distribution with mean 0 and variance being 1
a <- rnorm(100, 0, 1)
a
##   [1]  2.0946337942  0.6554531425  0.8380659292 -0.9475945647 -2.3514440081
##   [6] -1.6122461388  1.3189475181 -2.6763645010 -0.8899014682  1.3202070338
##  [11]  1.3196816745  0.4404675251  0.0041998971  0.6937386854  1.0010217059
##  [16] -0.0641711023  0.6238180048  0.0004438392 -0.7026667789 -1.5037723829
##  [21]  0.9600142727  0.8024330473 -0.7802737095 -0.7698345228 -1.8436766229
##  [26]  1.0270551291  0.1962214601  0.0705065223 -0.7202309180 -1.2055164057
##  [31]  0.1006982879 -1.6292702915 -2.5957777140  0.4532173200 -0.1603272313
##  [36]  1.0079714417 -0.0128366742  0.3350258860 -0.4318295366 -1.8745025690
##  [41] -1.3460439102 -0.9439307976  0.2993822797 -0.2957475056  0.3980882626
##  [46]  0.9064310254 -0.0295384253 -0.4558752780  0.2831411252 -0.5597299169
##  [51]  0.9170813459 -0.8681504175  0.8002739925 -1.9125461109 -0.6113538372
##  [56] -0.1998753916  1.2556359075 -0.1515829362  0.9897277272  0.3655384962
##  [61] -0.1836337938  0.0634031522 -0.2755568594 -1.2442723708  0.2838318045
##  [66]  1.7505358804 -0.3229662018 -0.7440265731  0.4344032327 -1.1005365597
##  [71]  0.5075000019 -0.0377518960  0.8391601866 -1.5950507259 -0.0714672203
##  [76]  1.4117886950  0.2495724679  0.1977833122  1.0793356309  0.8106645950
##  [81]  0.8282520930  1.3652674128  0.0517836466  1.6313902450  0.2317626571
##  [86]  0.8415510019  0.5921659696 -0.4005686818  2.9544617786  0.8735210768
##  [91]  1.3454264587 -0.7996673018  1.5395249009 -1.1883787297 -0.4751636790
##  [96]  0.0203498645 -1.2979330041  1.0713656055 -1.5555828257 -0.6579126314






#Part 2: Operations with vectors and matrices
#Mathematical operations between vectors (as well as matrices) happen element-wise. For example:

a <- c(1, 2, 3)
b <- c(4, 5, 6)
d <- a + b
d
## [1] 5 7 9
#Here, d is a vector with elements being equal to the summation of the respective elements of a and b.

#The same is true for other operations. For example:

a <- c(1,2,3)
b <- c(4,5,6)
d <- a*b
d
## [1]  4 10 18
a <- c(1,2,3)
b <- c(4,5,6)
d <- a^b
d
## [1]   1  32 729

#If for example we want to create a vector of ten elements with the ith
# element being ii, then

a <- 1:10
b <- 1:10
d <- a*b
d
##  [1]   1   4   9  16  25  36  49  64  81 100

#--------------------------------------------------------------------------------#
#Part 1: Creating Vectors
#Vectors is a foundamental and very important data structures in R. Therefore, it’s very important to know ways to generate them.

#Using the : operator
#This is useful if we want to generate vector of consecutive integers

a <- 1:10
b <- 50:70
d <- 70:50
a
##  [1]  1  2  3  4  5  6  7  8  9 10
b
##  [1] 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70
d
##  [1] 70 69 68 67 66 65 64 63 62 61 60 59 58 57 56 55 54 53 52 51 50
#Notice, that in the last case of vector d, the order of the elements is from the greatest to the smallest.

#Using the seq function
#The seq function allows to define the first and the last elements of a vector as well as the step OR the number of elelements. In the later case, the elements are distributed at equal distances. For example:

## a vector for which the first element is 1, the step is 2 and the last is NOT GREATER than 20. Thus here: 
a <- seq(from=1, to=20, by=2)
a
##  [1]  1  3  5  7  9 11 13 15 17 19
## a vector from 1 to 20 with 30 elements in total
b <- seq(from = 1, to=20, length.out = 30)
b
##  [1]  1.000000  1.655172  2.310345  2.965517  3.620690  4.275862  4.931034
##  [8]  5.586207  6.241379  6.896552  7.551724  8.206897  8.862069  9.517241
## [15] 10.172414 10.827586 11.482759 12.137931 12.793103 13.448276 14.103448
## [22] 14.758621 15.413793 16.068966 16.724138 17.379310 18.034483 18.689655
## [29] 19.344828 20.000000
#Using the rep function
#The rep function is very useful because it allows to generate vectors with repetitive elements. For example, if we want to generate a vector of 10 1s then:

## a vector of ten 1s
a <- rep(1, 10)
a
##  [1] 1 1 1 1 1 1 1 1 1 1
## a vector of ten "R"
a <- rep("R", 10)
a
##  [1] "R" "R" "R" "R" "R" "R" "R" "R" "R" "R"
## a vector of ten "WORDS"
a <- rep("WORDS", 10)
a
##  [1] "WORDS" "WORDS" "WORDS" "WORDS" "WORDS" "WORDS" "WORDS" "WORDS"
##  [9] "WORDS" "WORDS"
#Using the c() function
#The c function (from the word concatenate), allows to concatenate elements to a vector (or vectors to another vector). For example:

## a simple example
a <- c(1,2,3,4,5,6,6)

## a new vector from two other vectors
b <- c(c(1,2,3), c(34,5,6))
#It is of course possible to combine these methods

a <- 1:10
b <- rep(1:5, 4)
d <- c(a,b)
d
##  [1]  1  2  3  4  5  6  7  8  9 10  1  2  3  4  5  1  2  3  4  5  1  2  3
## [24]  4  5  1  2  3  4  5
#Vectors of random numbers
#Since R is a statistical language, it is very useful to be able to draw random numbers from distributions

#a vector of random elements from a uniform
## 100 random values from a uniform distribution with min:0 and max:10
a <- runif(100, 1, 10)
a
##   [1] 8.444331 4.879045 6.233414 2.148667 7.246836 6.454051 2.910835
##   [8] 8.820577 4.774129 6.535367 6.117422 6.130629 1.249406 6.401247
##  [15] 9.337773 5.649142 3.982129 8.555728 3.121173 7.688253 8.907029
##  [22] 4.557183 3.568504 2.103356 7.497209 4.652765 9.479190 8.462343
##  [29] 6.657565 5.263065 5.590213 6.501272 3.268033 6.860039 8.733813
##  [36] 2.559283 1.088549 3.111928 7.748140 3.503235 7.008526 1.692579
##  [43] 3.608574 9.371586 8.385626 2.580697 9.853283 2.233135 7.010064
##  [50] 9.401486 8.593154 5.235938 2.418401 2.544461 5.072425 9.817651
##  [57] 1.921284 1.909572 6.414346 9.267274 2.124362 9.554967 8.626550
##  [64] 1.230292 6.205839 2.509679 9.587234 7.188495 8.268828 6.530076
##  [71] 4.616666 3.049735 2.350251 5.737566 1.182483 5.990748 1.806069
##  [78] 9.611678 9.793944 3.199450 9.216586 3.822633 9.023409 3.051067
##  [85] 8.196802 1.659558 9.354643 8.443624 2.623373 6.044032 6.620612
##  [92] 6.916750 8.893715 7.307494 4.598399 8.094269 2.051436 4.914949
##  [99] 1.065808 5.454927
a vector of random values from a Gaussian distribution with mean 0 and variance being 1
a <- rnorm(100, 0, 1)
a
##   [1] -0.53180169  1.47357771 -1.29260609  0.61642095 -0.61926486
##   [6] -0.81908122 -0.30357929 -0.25037962 -0.55302262  1.69581395
##  [11]  0.64504865 -0.46468933  1.43111999 -0.83727716 -1.15911889
##  [16] -0.71988558 -1.29549666  0.65930801  1.73312496  1.30240243
##  [21]  0.86416889  3.57784473  1.85709327  0.56271396  0.49510954
##  [26] -0.49726317 -1.55979201  0.30278514 -0.32711513  0.54820943
##  [31]  0.70101236 -1.05513756 -0.64340282  0.43506756  1.33102245
##  [36]  0.82384158  1.31171070  1.12121884  0.44482299 -1.60720159
##  [41]  0.31112010 -0.45976606 -0.34951627 -1.23073440  1.44470932
##  [46] -1.02279711  0.74433349  1.27165979  0.37776969  0.96243470
##  [51] -0.47112515  1.41794058  0.45564411  0.42304487  0.27116125
##  [56]  0.12819769  0.72459963 -0.25831427 -1.23732378  0.14497209
##  [61]  1.25666167  1.10205660 -0.99839236  0.20934090 -0.91838769
##  [66] -2.03011078  1.05324937 -0.14538330  0.27083736  2.08628874
##  [71] -0.72981285  0.94409337  0.38824228  1.11839411  2.93109894
##  [76] -1.68473500 -0.98810788 -0.50056649  1.24220397 -0.97687030
##  [81]  0.96730106 -0.87373341 -1.00178756 -0.04379972 -2.61452528
##  [86]  1.50226917  0.26396248  0.07160694 -0.27387154  1.65248996
##  [91] -0.27822856  1.40147971  1.28874882 -1.30956512  1.51297045
##  [96]  0.77357950  0.67787305 -0.25526123  0.08368234  1.47194432
#Part 2: Operations with vectors and matrices
#Mathematical operations between vectors (as well as matrices) happen element-wise. For example:

a <- c(1, 2, 3)
b <- c(4, 5, 6)
d <- a + b
d
## [1] 5 7 9
#Here, d is a vector with elements being equal to the summation of the respective elements of a and b.

#The same is true for other operations. For example:

a <- c(1,2,3)
b <- c(4,5,6)
d <- a*b
d
## [1]  4 10 18
a <- c(1,2,3)
b <- c(4,5,6)
d <- a^b
d
## [1]   1  32 729
#If for example we want to create a vector of ten elements with the ith element being ii, then

a <- 1:10
b <- 1:10
d <- a*b
d
##  [1]   1   4   9  16  25  36  49  64  81 100
