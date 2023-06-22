#Until now, we have seen several data types in R. To see what type is your data, you can use the class() or the str() functions.

#For example:

vec <- c(1,5,7);mec <- 1:10;dec <- c("TRUE","TRUE","FALSE","FALSE","FALSE")
l<-list(vec,mec,dec)
#What class is my data? The class() and str() functions
#Too often are we faced with the problem of not realizing the class of the data we are handling. This is especially more troubling in the case of data frames and lists whose components may be of different classes. The class() function is called upon any data object and returns the data type. In the case of the above list l

class(l)
## [1] "list"
#class(l) returns “list” which is exactly what l is. Now what if we wanted to know what is the data type of each of the object in l? In this case we may use the str() function

str(l)
## List of 3
##  $ : num [1:3] 1 5 7
##  $ : int [1:10] 1 2 3 4 5 6 7 8 9 10
##  $ : chr [1:5] "TRUE" "TRUE" "FALSE" "FALSE" ...

#str(l) returns a much more detailed output that contains the data class of l (List), the number of objects in it (of 3) and the data type of each of the objects (num, int, character) alongside the first instances in each one. Notice how the last object of l has been assigned the type “character” (chr). We can easily coerce it back to logical and put it in the list with

l<-list(vec,mec,as.logical(dec))
str(l)
## List of 3
##  $ : num [1:3] 1 5 7
##  $ : int [1:10] 1 2 3 4 5 6 7 8 9 10
##  $ : logi [1:5] TRUE TRUE FALSE FALSE FALSE



#Subsetting
#Subsetting refers to the selection of parts of data from greater sets. Subsetting of data types is one very important aspect of the R environment in the sense that it can be performed with extreme precision and at great speeds. In this sense it constitutes one of R’s main advantages. Subsetting uses a number of special characters to perform various tasks, such as obtaining specific rows, columns, elements from all possible data types depending on the user’s choise. It can be roughly divided to:
#-structural subsetting, where data are subsetted based on the structure of the data type (e.g. the 3 first columns)
#-logical subsetting, where data are subsetted based on a logical restriction (e.g. all values that are not “NA”)
#-numerical subsetting, where data are subsetted based on numerical/categorical operation/control (e.g. all values >10 or all values equal to “FALSE”)

#Two subsetting operators [] and [[]]
#The single [] subsets everything (arrays, lists, matrices, vectors) and returns a datatype of the same type as the original object.

#On the contrast, the [[]] returns an object with the same datatype as the ‘inner’ object. An example will make things more clear.

l1 <- list(c(1,2,3), c(7,8,1))
l1[1]
## [[1]]
## [1] 1 2 3
## l1[1] is  a list because l1 is a list and we used the [] operator
class(l1[1])
## [1] "list"
l1[[1]]
## [1] 1 2 3
class(l1[[1]])
## [1] "numeric"
Of course we can combine [[]] and [] in lists.

l <- list(c(1,2,3), c(1,4,5))
l[[1]][2] ## this is 2
## [1] 2

#Why the class of an object is important?
#Several functions use arguments that are restricted to a specific class. For example, the function prcomp accepts only matrix and not a data.frame.

#Let’s see the classes of the main data objects

## matrix
m <- matrix(1,1,1)
class(m)
## [1] "matrix"
## list
l <- list(1)
class(l)
## [1] "list"
## vector
v <- c(121)
class(v)
## [1] "numeric"
##array
a <- array(1, dim=c(1,1,1))
class(a)
## [1] "array"
## data.frame
df <- as.data.frame(m)
class(df)
## [1] "data.frame"
#Note: there is no class vector. It appears as the type of data it handles. For example:

v1 <- c(1,2,3)
class(v1)
## [1] "numeric"
v2 <-c(T,F,T,F)
class(v2)
## [1] "logical"
v3 <- c("A", "B")
class(v3)
## [1] "character"


#The $ operator
#We have already seen that data.frames may have column names, lists may contain named objects. For example:

a <- data.frame(col1=c(1,2,3), col2=c(4,5,6))
a
##   col1 col2
## 1    1    4
## 2    2    5
## 3    3    6
a$col1
## [1] 1 2 3
a$col2
## [1] 4 5 6
a[,1]
## [1] 1 2 3
a[,2]
## [1] 4 5 6
#The same applies in lists

l <- list(first=c(1,2,3), sec="XXXX")
l$first
## [1] 1 2 3
l$sec
## [1] "XXXX"

#Structural Subsetting
#Makes use of the [], [[]] and $ operators, which with the clever combination of commas can provide absolute precision on the choise of data with only a few characters coding. [] may be used on any vector, factor, matrix or dataframe to subset it in one or two dimensions. In the case of vectors and factors there is only one dimension. Therefore, if we want the nth element of a vector v we simply put n within brackets

v<-1:10
x<-v[6]
#x now holds the sixth value of the vector v. R enumerates all data types starting from 1 (and not 0 like Perl or Python) so 6 will actually return the 6th value.

x<-v[6:8]
#will get a “slice” of v and store it in x which now becomes a vector itself carrying the 6th,7th and 8th elements of v. Remember how the “:” operator is used for ordered integers. Now what if we wanted some compartmentalized subsetting that does not follow a certain order

x<-v[c(1, 6:8, 11:15, 18, 20)]
#This intricate subsetting allows us to get the 1st, then the 6th-8th, then the 11th-15th, then the 18th and then the 20th elements of v and store them in a vector called x. Notice how the subsetting indices (the numbers in the parentheses) are a vector in themselves and thus they are introduced with the concatanate function c(). We could have greater control in this subsetting if we split the process in two

m<-1:100
ind<-c(1, 6:8, 11:15, 18, 20)
x<-m[ind]
#This first creates a vector called ind that carries the indices (the numbers of elements we want to obtain from v) and then passes it to v with [] to perform the subsetting. The exact same process stands for factors (which as we already know are categorical vectors). But what about matrices and dataframes? Here there are two dimensions on which we can subset (rows and columns). R uses the same operator [] but allows for two values separated by comma to provide information of rows and columns (in this order). Although both of them are not always needed (suppose you only want to subset columns but not rows) R needs to keep in mind we are treating two-dimensional datatypes so we need to use a comma inbetween. This will become less confusing with an example. Suppose we need to keep only the first and the third line from a matrix m. This is done with:

m<-mtcars
mm <-  m[,c(1,3)]
#Notice that the indices inside the brackets are of the form [ ,vector]. That is because the value before the comma is reserved for row subsetting. Since we don’t want to subset on the rows we leave this empty, but use the comma since this is compulsory. After the comma we simply provide the vector of indices we want to subset columns by (in this case c(1,3) for the 1st and the 3rd). In perfect symmetry subsetting on rows 10 through 20 would be performed with:

m<-mtcars
mm <- m[10:20,]
#As in this case the indices are serial we need not use a concatanated structure so 10:20 will do. Alternatively we could use c(10,11,12,13,14,15,16,17,18,19,20) but we don’t for obvious reasons. Bear in mind that in both cases above mm is a matrix (or data frame, depending on what m was) whose dimensions have now changed. If m was a MxN matrix, then mm is a Mx2 in the first case and 10xN in the second. In the case we subset on both rows and columns

m<-matrix(1:140, nrow=10, ncol=14)
mm <- m[c(1:5,8), c(10,11, 12:14)]
#mm is now a 6x5 matrix. Understandably we can return any single element of a matrix or data frame by providing its exact “coordinates” and thus

m[6,9]
## [1] 86
#will return the 6th element of the 9th column. In the case of data frames with named columns we can also use the $ operator to subset columns. Take the built-in R data frame called mtcars simply by typing

mtcars
##                      mpg cyl  disp  hp drat    wt  qsec vs am gear carb
## Mazda RX4           21.0   6 160.0 110 3.90 2.620 16.46  0  1    4    4
## Mazda RX4 Wag       21.0   6 160.0 110 3.90 2.875 17.02  0  1    4    4
## Datsun 710          22.8   4 108.0  93 3.85 2.320 18.61  1  1    4    1
## Hornet 4 Drive      21.4   6 258.0 110 3.08 3.215 19.44  1  0    3    1
## Hornet Sportabout   18.7   8 360.0 175 3.15 3.440 17.02  0  0    3    2
## Valiant             18.1   6 225.0 105 2.76 3.460 20.22  1  0    3    1
## Duster 360          14.3   8 360.0 245 3.21 3.570 15.84  0  0    3    4
## Merc 240D           24.4   4 146.7  62 3.69 3.190 20.00  1  0    4    2
## Merc 230            22.8   4 140.8  95 3.92 3.150 22.90  1  0    4    2
## Merc 280            19.2   6 167.6 123 3.92 3.440 18.30  1  0    4    4
## Merc 280C           17.8   6 167.6 123 3.92 3.440 18.90  1  0    4    4
## Merc 450SE          16.4   8 275.8 180 3.07 4.070 17.40  0  0    3    3
## Merc 450SL          17.3   8 275.8 180 3.07 3.730 17.60  0  0    3    3
## Merc 450SLC         15.2   8 275.8 180 3.07 3.780 18.00  0  0    3    3
## Cadillac Fleetwood  10.4   8 472.0 205 2.93 5.250 17.98  0  0    3    4
## Lincoln Continental 10.4   8 460.0 215 3.00 5.424 17.82  0  0    3    4
## Chrysler Imperial   14.7   8 440.0 230 3.23 5.345 17.42  0  0    3    4
## Fiat 128            32.4   4  78.7  66 4.08 2.200 19.47  1  1    4    1
## Honda Civic         30.4   4  75.7  52 4.93 1.615 18.52  1  1    4    2
## Toyota Corolla      33.9   4  71.1  65 4.22 1.835 19.90  1  1    4    1
## Toyota Corona       21.5   4 120.1  97 3.70 2.465 20.01  1  0    3    1
## Dodge Challenger    15.5   8 318.0 150 2.76 3.520 16.87  0  0    3    2
## AMC Javelin         15.2   8 304.0 150 3.15 3.435 17.30  0  0    3    2
## Camaro Z28          13.3   8 350.0 245 3.73 3.840 15.41  0  0    3    4
## Pontiac Firebird    19.2   8 400.0 175 3.08 3.845 17.05  0  0    3    2
## Fiat X1-9           27.3   4  79.0  66 4.08 1.935 18.90  1  1    4    1
## Porsche 914-2       26.0   4 120.3  91 4.43 2.140 16.70  0  1    5    2
## Lotus Europa        30.4   4  95.1 113 3.77 1.513 16.90  1  1    5    2
## Ford Pantera L      15.8   8 351.0 264 4.22 3.170 14.50  0  1    5    4
## Ferrari Dino        19.7   6 145.0 175 3.62 2.770 15.50  0  1    5    6
## Maserati Bora       15.0   8 301.0 335 3.54 3.570 14.60  0  1    5    8
## Volvo 142E          21.4   4 121.0 109 4.11 2.780 18.60  1  1    4    2
#This small dataframe contains makes of cars alongside their constructor specifications. In order to choose one specific specification simply ask for the name of the data frame followed by $ and the name of the column. For instance

mtcars$cyl
##  [1] 6 6 4 6 8 6 8 4 4 6 6 8 8 8 8 8 8 4 4 4 4 8 8 8 8 4 4 4 8 6 8 4
#will return the vector containing the cyl column. This is identical to calling mtcars[,2] asking thus for the second column of the data frame (not including the names of the cars). The [[]] and $ operators are mostly used in list context. Although, as we saw earlier, [n] can be used to invoke the nth element of a list, [[n]] does the same without getting back the reference (the name) of that element. This is sometimes desirable, especially when we want to pass data from a list to another function. Consider the list l we saw earlier

l<-list(vec,mec,dec)  
#where vec, mec and dec are three vector of different type and size. l can be created by keeping names for each of the three

l<-list(a=vec, b=mec, c=dec)
#now vec can be invoked with all of the following commands

l[1]
## $a
## [1] 1 5 7
l[[1]]
## [1] 1 5 7
l$a
## [1] 1 5 7
#There are subtle differences are in the form of the output. From top to bottom the output is stripped from (sometimes unneccesary) references.



#Logical Subsetting
#Makes use of logical operators (“&”, “!”, “|”, see more on that in Control Structures) Logical operators stand for AND (=“&”), OR (=“|”) and NOT (“!”). While the first two are more complex as “joining” operators and we will see more of them when we discuss control structures, NOT “!” can stand on its own as it signifies the negation of a statement. In this sense

!is.na(x)
##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#returns logical values (TRUE or FALSE) depending on whether any value in x is NOT “NA”. In this sense, !is.na() is the mirror-image of is.na(). Lets use this in a subset

y<- x[!is.na(x)]
#y is now a vector with all the elements of x that are NOT “NA”. There are two things to be careful of here. One, that the subsetting index refers to the subsetted data object (x). Notice how x is subsetted based on one of its own properties (how many of its elements are not “NA”). The other is that the output is a vector, regardless of the structure of x. Even if x was to be a matrix or a data frame, the logical subsetting will return the values reading by column starting from the top left.


#Numerical Subsetting
#Makes uses of the common numerical operators (>, <, >=, <=, “!=” and ==). You should already be familiar with “>”,“<” and “<=”, “=>” as “greater”, “smaller”, “smaller or equal” and “greater or equal”, but notice “!=” signifying “not equal to” and “==” for “equal to” (this is because in R as in most programming languages, we use “=” for assignment of variable values, (in R “->” and “=” are the same but we strongly recommend “->” to avoid confusion)). Numerical subsetting works exactly like logical. So assuming x is a data frame holding these values

x<- data.frame(a=1:4, b=c("Me","You","Him","Her"))
#Then, subsetting it numerically by asking only values equal to 2 would be

x[(x==2)]
## [1] "2"
#This returns the value “2” (only once in this case, but more times if 2 were to be found more than once). What if we wanted to get values that are greater than 2. We could try

x[(x>2)]
## Warning in Ops.factor(left, right): '>' not meaningful for factors
## [1] "3" "4" NA  NA  NA  NA
#which will give a > Warning message: In Ops.factor(left, right) : > not meaningful for factors

#What happened here? We got the two numerical values that fulfill our restriction (>2) but we got four NA values at the end and an error message that said our operation is not meaningful for factors. What R is trying to tell us is that a comparison >2 is not possible for categorical data like “Me” or “You”. In this sense the operation “Me”>2 returns nothing (NA). R prints this at the end of the vector but is kind enough to point it out to us. Numerical subsetting is not always numerical. It can be categorical as well with the use of the “==” and the “!=” operators. If we try

x[(x!="Me")]
## [1] "1"   "2"   "3"   "4"   "You" "Him" "Her"
#we get back all the elements of x that are not equal to “Me” regardless if they are numbers or characters. This is because R performs coercion of variables wherever this is possible before conducting the subsetting. This is very handy in the general cases of data comparisons that we discuss next in Control Structures. One very handy function for subsetting is which(). It can be used for both logical and numerical subsetting to produce a subset of indices fulfilling certain conditions.

x<-c(1,-1, 0, 3, -20, -2)
which(x>=1)->ind
x[ind]->new_x
#which() can also be used in logical context.


#An example. Subsetting for not NA values
#Lets now use subsetting in an example that is very useful (and also quite common). Removing NA values from a data frame. The task we want to complete is, given a data frame (or a matrix), extract the instances (rows) that do not have NA values (holes). As with many cases, R already has a built-in function to check for lines in a matrix that do not carry NA values. This function is complete.cases() and can be invoked on a matrix m like this

ind<-complete.cases(m)
#ind now is a vector that holds the line numbers of the data frame that fulfill the condition (not having a NA value). All we have to do now is to ask for a subset of m with the rows held in ind

mm<-m[ind,]
#And we are done! A clean data set.



#And finally. The subset function
#A number of things in R can be done with the use of predefined functions. Subsetting is not an exception. There is a specific subsetting function called subset(). subset() combines the use of all types of subsetting (numerical, logical and structural) in data frames with named columns. It also makes use of logical operators to combine subsetting commands in a single. An example may be seen with a default R dataset called “airquality”. “airquality” is structured data in a data frame concerning information on ozone, solar radiation, wind and temperature for a number of dates organized by month and day. To have a better view of its contents simply type

head(airquality)
##   Ozone Solar.R Wind Temp Month Day
## 1    41     190  7.4   67     5   1
## 2    36     118  8.0   72     5   2
## 3    12     149 12.6   74     5   3
## 4    18     313 11.5   62     5   4
## 5    NA      NA 14.3   56     5   5
## 6    28      NA 14.9   66     5   6
#Notice how there are holes in the data with a number of NA values. Suppose now that we want to obtain a slice of the data that we contains all dates with a temperature higher than 60 degrees and a wind of 10 knots or more. We would type

dates<-subset(airquality, Temp > 60 & Wind>=10)
#Observe how the function works. We call subset() on the dataframe airquality asking that a combined condition is fulfilled, so as Temp>60 AND Wind>=10. The logical “AND” is coded by the ambersand “&” symbol. Alternatively, had we wanted to keep the dates with either Temp>60 OR Wind>=10 we would have asked for

dates<-subset(airquality, Temp > 60 | Wind>=10)
#in which case the logical “OR” is coded with the bar symbol “|”. Think about the case where we would have wanted mutual exclusion of conditions (e.g. Temp>60 “AND NOT” Wind>=10). In this case we would have to think a bit more and code for an equivalent condition. That would be Temp>60 & Wind<10 (inverting the condition on wind). Finally, subset can also incorporate structural subsetting in the form of retaining specific parts of the data. In the case e.g. that we would have wanted to keep only the dates (that is month and day) fulfilling the above condition (that is dropping the meteorological data) all we need to do is to use the select argument of subset(). The command would now be:

dates<-subset(airquality, Temp > 60 & Wind<10, select = c(Day, Month))
head(dates)
##    Day Month
## 1    1     5
## 2    2     5
## 7    7     5
## 10  10     5
## 11  11     5
## 12  12     5
#Notice how we have complete freedom to manipulate the data frame in terms of ordering of its vectors. In this example we choose to show Days before Months although their order was the other way round in the initial data frame.



#Simple Functions and Control Structures
#One of R’s main powers is the great number of predefined functions. Taken together the set of built-in functions alongside those contained in various R packages constitute an extensive toolbox with which you can perform mathematical calculations, conduct complex statistical analysis and create elegant and highly informative graphs. A detailed summary of the complete “dictionary” of R functions is beyond the scope of these notes (and may well be beyond the scope of most of R manuals). Nonetheless you can find the index of R functions contained in the base “core” distribution here.
#For the purposes of our classes we will discuss only a subset of all available functions, focusing on basic mathematical operations, functions that deal with basic (and a little more advanced) statistics and plotting functions.



#An overview of operators in R
#We have already seen how mathematical calculations may be performed in R in more or less the way one uses a calculator. Remember that we can simply type the number of variables with operators inbetween and get the result by hitting return (enter)

a<-5
b<-6
a + b
## [1] 11
#The basic mathematical operators are recapped in the list below

#Arithmetic Operators * + addition * - subtraction * * multiplication * / division * ^ or ** exponentiation * x %% y modulus (x mod y) 5%%2 is 1 * x %/% y integer division 5%/%2 is 2 Out of those you may not be familiar with the last two. Modulus provides the remainder of a division between two integers, while integer division is the result of a division without proceeding further than the integral part of the quotient.


#Logical Operators
#Logical operators have been discussed previously in the chapters related to Subsetting. A formal definition of logical operators (or connectors, or connectives) is a set of symbols that may be used to connect two or more statements in a grammatically valid way. In the example discussed above for the subset() function these two sentences would be: a) “The temperature is greater than or equal to 60 degrees” and b) “The wind speed is greater than 10 knots”. Each of the two sequences contain some logical operators in itself. The temperature>=60 and the wind>10 are examples of numerical control operators. Combinations thereof consist of asking: * Both of them being true (A AND B) * Only one of them being true (A AND NOT B, B AND NOT A) * At least one of them being true (A OR B) All of the above can be coded with specific logical operators listed below. * < less than * <= less than or equal to * > greater than * >= greater than or equal to * == exactly equal to * != not equal to * !x Not x * x | y x OR y * x & y x AND y In the following chapters we will see how we go about using them and so getting used to working them on a routine basis


#R Functions
#Let’s now get to the core of the basic R function toolbox. As discussed previously R almost has a function for everything. A number of these functions are very useful for getting beginners accustomed to the way R works and to prepare them for writing their own functions (when they stop being beginners).
#Below we present some of the most important R functions for various objectives covering basic statistics, plotting and string manipulation.

#Numeric Functions
#Numeric functions include functions used to performed more advanced mathematical operations. Some of the most important are:

abs(x) : absolute value
#abs() is used to return the absolute value of a number (integer or not). This means that both abs(-5) abs(5) return the value 5.

sqrt(x) : square root sqrt() returns the square root of a number.

#ceiling(x), floor(x), trunc(x), round(x, digits=n), signif(x, digits=n) All of these functions treat real numbers in terms of rounding (that is ommission of decimal points). ceiling(x) rounds up the real number x to the closest integer that is greater than x (>x) while floor(x) does the opposite, rounding x to the closest integer that is smaller than x (<x). In this sense

ceiling(5.1)
## [1] 6
#returns 6, while

floor(5.98)
## [1] 5
#returns 5

#trunc(x) truncates all decimal points rounding the number to the closest integer in the same way floor() does it. round() and signif() both take a number of digits as an additional argument but slightly differ in that round(x, digits=n) does rounding in a way that keeps n decimal points while signif(x, digits=n) rounds to n total digits (not including the comma). Thus:

x<-45.6789
round(x, digits=3)
## [1] 45.679
signif(x, digits=3)
## [1] 45.7
#cos(x), sin(x), tan(x), acos(x), cosh(x), acosh(x) Refer to the corresponding trigonometric functions for cosine, sine, tangent, arc-cosine, hyperbolics etc.

#log(x) : natural logarithm, log10(x) : common logarithm
#log() is the natural logarithm while base-10 log is coded as log10(). Do you remember what you need to do to convert a natural log to, say, a base-2 logarithm? In case you don’t remember how to change the base you can always use logb(x, base=n)

#exp(x) : exponential of e remember generic exponentiation is coded either through a^x or a**x.


#Character Functions
#Character handling is not one of R’s major strong point and you can always work around data manipulation in character strings outside R with shell scripting or other scripting languages. Still, R provides a number of functions we can use to handle strings when we need to stay within the environment.

substr(x, start=n, stop=m)
#The substr() function works more or less in the way the function of the same name in Perl. In a string x it returns a substring starting from n and running through m. Remember that numbering in R starts from 1.

grep(pattern, x , ignore.case=FALSE, fixed=FALSE)
#Pattern matching in R is performed with grep() that searches for pattern in x. Using fixed=F assures that pattern search may be performed with a regular expression. Notice that the functions returns the matching indices, that is the result of the matching is the subset of elements of the vector x that match the pattern. A substitution function sub() is similar to grep but with the addition of a replacement string sub(pattern, replacement, x, ignore.case=F, fixed=F) sub(pattern, replacement, x, ignore.case =FALSE, fixed=FALSE) for instance

sub("\\s",".","Hello There") 
## [1] "Hello.There"
#Splitting of a string is performed with strsplit() strsplit(x, split) where split is the character(s) at which splitting takes place. strsplit() returns a character vector. As in the previous functions, fixed=T allows split to be a regular expression. For instance:

x<-"abc.d.eef.g"
strsplit(x,".", fixed=T)->d
d
## [[1]]
## [1] "abc" "d"   "eef" "g"
#The paste() function joins strings in a concatenation using sep as a separator

paste("x",1:3,sep="") 
## [1] "x1" "x2" "x3"
paste("x",1:3,sep="M") 
## [1] "xM1" "xM2" "xM3"
#Transforming strings to upper or lower case is explicitly done with toupper(x) and tolower(x).

#Other Useful Functions
#A number of very useful functions do not fall in a specifi category but will prove very handy once you start coding your own R scripts and functions. Suppose you need to generate a sequence of numbers with a fixed interval. Remember that this can be done with the “:” operator only for interval=1 but if we want another step we may use seq(from , to, by) If you code

x <- seq(1,20,3)
x
## [1]  1  4  7 10 13 16 19
#x becomes the vector of all values down to the largest that fulfils the condition < to-by Another type of sequence we may need to create is the repetition of elements. This is performed with the use of the rep() function rep(x, ntimes)

y <- rep(1:3, 2)
y 
## [1] 1 2 3 1 2 3
#Obtaining a random set of values from a greater set is done with sample(). sample(x, size, replace=F/T) sample() takes a subset of size=size from the vector x and returns it a smaller vector. If replace=T then the same element can be drawn more than once.

pretty(c(start,end), N) #returns a vector of equally spaced N values between start and end. Invaluable for simulating distributions, producing values for standard reference plots etc.

sort(x) #returns the vector x sorted in numerical order from the smallest to the largest element while

order(m[,c(i,j)]) orders a matrix m according first ot the values in the i-th and then in the j-th columns. Calling

m<-matrix(runif(100), nrow=50, ncol=2)
i<-1; j<-2
m[order(m[,i],m[,j]),]
##             [,1]       [,2]
##  [1,] 0.01681164 0.44279025
##  [2,] 0.02747127 0.33867720
##  [3,] 0.06231898 0.01018724
##  [4,] 0.11502784 0.45307732
##  [5,] 0.12628696 0.87329223
##  [6,] 0.14427334 0.98702452
##  [7,] 0.17735398 0.28601187
##  [8,] 0.18345233 0.09093880
##  [9,] 0.19853615 0.76766366
## [10,] 0.23524911 0.42016889
## [11,] 0.24598157 0.39011060
## [12,] 0.24664202 0.98766170
## [13,] 0.32960463 0.41007013
## [14,] 0.33213535 0.80586546
## [15,] 0.33648320 0.66095139
## [16,] 0.34764840 0.14601567
## [17,] 0.36429720 0.81331753
## [18,] 0.38955678 0.07215733
## [19,] 0.41853466 0.74105154
## [20,] 0.46363333 0.20020714
## [21,] 0.46474820 0.47704310
## [22,] 0.49868507 0.47445452
## [23,] 0.50358894 0.13286459
## [24,] 0.57433237 0.21894682
## [25,] 0.58258469 0.53161126
## [26,] 0.59675955 0.20117120
## [27,] 0.59712366 0.92219224
## [28,] 0.60729780 0.94107224
## [29,] 0.64163243 0.13822350
## [30,] 0.65285012 0.05774504
## [31,] 0.65616172 0.06775344
## [32,] 0.68971759 0.64666866
## [33,] 0.70531491 0.45574166
## [34,] 0.73645609 0.98786753
## [35,] 0.73983130 0.98634404
## [36,] 0.75310690 0.16646877
## [37,] 0.78081686 0.36546804
## [38,] 0.78204270 0.13208013
## [39,] 0.78260930 0.45575914
## [40,] 0.80338747 0.97839731
## [41,] 0.82253888 0.04585821
## [42,] 0.85882987 0.81417212
## [43,] 0.86706878 0.63837168
## [44,] 0.88925306 0.82480493
## [45,] 0.91049092 0.31693518
## [46,] 0.91474325 0.68112694
## [47,] 0.92281462 0.31350792
## [48,] 0.93960475 0.73392726
## [49,] 0.96470141 0.34136102
## [50,] 0.97641055 0.10618281
#will produce a re-ordered matrix according to the above ordering. (More about the very useful runif() function later on)
