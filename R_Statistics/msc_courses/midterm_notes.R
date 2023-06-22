## Random variable = Biological Data
## Concept of capital letter, means we cant know the exact value of a person from our data until we pick that person. That variable is specific to each person
##For a second random variable, if all students had 0 in their grades, then it doesnt matter how many hours they spent studying, so it is not a random variable.
##Combine 2 random variables into a new random variables. For example grade/hours

## ?kati --> Usage/Help
## help.search("kati") --> Help page
## x = 1:10 --> Vector of integers

##Scalars. R does not have scalars, but instead VECTOR (numeric vector), of length 1 if num
a = 1
b=2
a

class(a) #
typeof(a) # Double is the deffault type. This means its a float? and not an integer
length(a)

a = 1L #make the float into integer
b = 2L
typeof(a)

a/b #Double division
a%%b #Integer divison

b^3 or b**3 #Power

length(b**3)
typeof(b**3)

## String:
t = "ATCATTGGA"
class(t) #both are same
typeof(t)
length(t) #It shows 1 because its a vector of length 1.
nchar(t) #This show the length of the string

##Logical:
g = TRUE
r = F

class(r)
length(g) #Print 1 becuase its a vector
nchar(g)

## & (and, means both arguments TRUE), &&, | (or, when one is TRUE), ||


##
x = as.integer(3)
typeof(x)

as.numeric // as.string #Works for TRUE / FALSE


#### VECTORS #### 1 Dimensional constructions with 1 type of variable (only strings, only num

a = c(1,2,3) #All are doubles
a
length(a)

b = c("a", "fkafs", "fasrW") #Type = Character
b
length(b)

aa = c("a",2,3) #Will convert all to string.
aa

aa = 1:10
class(aa) #integer

bb = c(1,2,3,4,5,6,7,8,9,10)
class(bb)

##functions that return vectors

aa = seq(1,10,2) #(Start, End, Step)
aa

bb = seq(1,10,length.out=11)
bb


runif(10,10,20) #r = random, unif = uniform. (NumberCounts, Start, End). Randomly randomly uniformly distributed

rnorm(10,10,20) #Gaussian distribution (NumberCounts, Mean, SD) #Mean = Average, SD = How much they deviate from the mean

sample(c(TRUE,FALSE), 100, replace = TRUE) # Random from (vector, NumCount, with replacement)
sample(1:1000, 100, replace = TRUE) #Uniform distribution

rep(10, 5) #(Number i want to repeat, Times)
rep(c(10,2),each = 5) #Each will make every element repeat 5 times)
rep(c(10,2,3),c(5,10,3))


## Operations

c(1,2,3) * c(3,4,5) #Element-wise operations
c(1,2,3) + c(1,2) # It will add and when it is out of numbers it will cycle the small vec
c(1,2,3,4) + c(1,2) #No warning message, full cycle

(0:20)*2 #Best way of operations.

c(1,2, c(1,2)) #This concatenates the vectors
a = 1
a = c(a,2)

a = c() #Empty vector
a = c(a,2)

a = c()
for(i in 1:10){
    a = c(a,i*2)
}


a = vector(mode = "numeric", length=10)
for(i in 1:10){
    a[i] = i*2 #Better than appending
}


2^(0:10) #All the powers of 2
2^(10:0) #All the powers of 2 reversed



### MATRIX ### 2D : NxM (only 1 type of char)

M[2,2] #Get the element of a matrix
M[2,2] = 5 #Set the element of a matrix

M[,n] #get all the Nth column --> THIS IS A VECTOR : c(n1,n2,n3,n4,..)
M[n,] #get the Nth row of all columns --> THIS IS A VECTOR
M[,2:3] or M[,c(2,3)] # Get consecutive columns OR get all the columns in the vector
M[,seq(from=2,to=ncol(M), by=2)] #Get columns starting from 2, and get every 2 consecutive
M[,(1:5)*2] #Multiply it by 2
M[sample(1:nrow(M), 10, replace=TRUE),] # 10 random numbers from 1:nrow(M) with replace.
##It will give 10xncol(M), with the rows beeing the vector from sample

v = 1:88 #Vector 1:88
matrix(v, nrow=20,ncol=5) #This will fill a column first then go to the next. If its over, it will start over

matrix(1:100,nrow=20,ncol=5,byrow=TRUE) #This will fill row first then move to the next
matrix(0, nrow=10,ncol=4) #Make a zero table

vv=matrix(0:3,nrow = 3, ncol=3)
vv
vv>0 #Gives matrix with True / False
vv[vv>0] #Gives vectors with the values that fulfill requirements


### APPLY ###

a = list(a=c(1,2,3), b=runif(100, 1, 10), c=c(T,F,T))
head(a)

lapply(a, sum)#an thelw na einai pali lista.Pairnei orisma ath lista, kai sth synexeia th sunarthsh poy tha efarmosw se kathe stoixeio thw listas

sapply(a, sum) #epistrefei vector. An den mporei na to kanei vector epistrefei lista

lapply(a, function(x){return(c(x[1], x[2]))})#tha parei diadoxika a,b,c

m = sapply(a, function(x){return(c(x[1], x[2]))})#epistrefei pio aplh domh, dhladh vector h matrix analogws ti voleyei. Edw epistrefei pinaka

apply(m, 1, function(x){log(prod(x))}) #auto to efarmozei se kathe grammh (1)

apply(m, 2, function(x){log(prod(x))}) #auto to efarmozei se kathe sthlh (2)

logprod = function(x){log(prod(x))}
apply(m, 2, logprod)

logprod_2 = function(x, threshold){
    y=log(prod(x))
    ##y[is.inf(y)]
    return(y<threshold)
}
apply(m, 2, logprod_2, threshold=3) #alliws mporw apla na balw to 3

##### LISTS #####
##Store multiple values like (Vectors, Matrixes, DataFrame).
##We can store different types of data

## [] --> Slice operator and returns the same structure (so list)
## [[]] --> Returns whatever is in that position of the list

a = list(b = 1, d = c(1,2,3))
a
a[[1]] #Value of the first element of the matrix
a[[2]][1]
a$b #Value of the element b of the list
a$d

a[2] # THIS IS A LIST
a[2][[1]][1] #The second elements of the first list, then the first element of that and the first index of that
a[1:2] #This is a list with the first 2 elements of the list a


de = a[c(1,1,1,1)]
de
typeof(de) #It is a list
length(de) #It has length of 4
names(de) #Names can be the same

de[[4]]
de[[1]]=3 #Changes first element of the list to 3
de

str = "fasghjkl"
strsplit(str,"")##Create a list from string

str2 = c("fasfafffg","asfef")
strsplit(str2,"") ##Create a list from that string

#### PWM ####

## Biology Stuff ##

##1)How much RNA will be produced
##2)When transcription will start
##Protein (Transcription Factors) binds before TSS (transcription start site) to start trans
##Regulation is achieved by TF, and binds to specific sequence or similar to that.
##It binds to a motif (motif = sequences that kind of match). A motif is usually strict at a position (means it allows only a specific nucleotide) while at others it allows more



###STEPS###
##1) Create 4xl zero matrix (l = motif length). Add pseudocount to each position (10^-8)
##2) Fill each position with occurences of each nucleotide --> COUNT MATRIX
##3) Convert to frequencies: Divide each column with its OWN sum (because we might be missing some information - N)
##4) Find frequencies of nucleotides by what we would expect by chance --> 0,25 For each OR we count all the [AGCT] of the whole DNA seq -----> DIVIDE previous matrix by the frequency, and then we LOG the matrix.
## i) If the freq of a nucleotide is higher than what we expect, then the PSSM value is positive, OTHERWISE its negative.
## ii) The higher the frequency, the higher the PWM value
## iii) There is a maximum value : log(1+pseudocount/expected) // minimum log(pseudocount/ex)
### iv) Assume that different positions of the motif are independent. 
## v) P(ACGTTA | PWM) = P(A^1 | PWM) * P(C^2 | PWM) * ... P(A^6 | PWM)

##5) Find a score: Score(ACGTTA | PWM) = PWM(A,1) + PWM(C,2)+...+P(A,6). THE SUMS OF THE                                                                            RESPECTIVE PWM VALUES
##Higher score --> Higher probability to observe given DNA sequence for this PWM. (or how much it matches)

##6) Split the whole string into substrings of length = motif, and add the score of each position to that and SUM it. Greater score --> More plausible binding

##7) Compare the scores we get with cases where we know we have binding. This is difficult sice we dont know what kind of "affinity" to use // OR with cases we DONT have binding. For this case we can use strength of binding = 0.
##The score has a gaussian distribution. P(S >= 10,6 | NO BIND) prob to get score s with prob greater than 10,6 when i have no binding. That is really small. We would say it does not fit the distribution of no binding scores. So it would have a low pvalue and thus we cant say it fits the distribution. 
##The idea here is to calculate scores in cases with no binding and compare our scores. If they are very high, then we reject the H0 : No binding, and assume that my case is a case of binding (pvalue would be very low)

##Approach 1) We assume that if we construct artifically (simulate) a DNA then no binding will occur (or very few positions can be assigned to binding. (sample((A,C,G,T),...). In order to do that then we need to make sure we take into account the positions of preferrence of each nucleotide (for example maybe T wants to be next to A)

##Approach 2) The majority of real DNA has no binding. Thus, random locations represent no binding.

##After each approach, we need to get the distribution of scores. We generate that by calculating the frequency score.
##The hypothesis testing: we get a score and check where it is in the distribution of the scores. It resembles gaussian distribution, and in this case we want high score, so we find the probability to find that score above the P(0,05) threshold (upper.tail = True), where that is based on the null hypotehsis that there is no binding, and thus we rejected the null meaning we probably have binding.

random = c(0.2, 0.3, 0.2, 0.3)
pwm = matrix(sample(1:100, 40, TRUE), nrow=4) #create 4x10 motif
ppm = apply(pwm, 2, function(x){x/sum(x)})#pwm
pwm = log2(ppm/random)#pssm
alphabet = c("A","C","G","T")
row.names(pwm) = alphabet

dna = sample(alphabet, 20000, replace = TRUE) #Lets construct a random 20K length DNA, where we know we probably wont have a binding site
dna #The best idea is to construct based on nucleotides like ("AA","AC","AG",...,"TT"). We can take into account here the propabilities to find each dinucleotide


allScores = vector("numeric", (length(dna)-ncol(pwm)+1))
#This function is to find the score for each motif in that substring.
for(i in 1:(length(dna)-ncol(pwm)+1)){
    substring = dna[i:(i+ncol(pwm)-1)]
    score = 0
    for(j in 1:ncol(pwm)){
        score = score + pwm[substring[j],j]
    }
    allScores[i] = score
}
allScores
hist(allScores) #With a hist we can see all the distribution of scores. We notice that majority is negative

plot(density(allScores)) #best way for histogramms. We Find where the P(0.05) threshold is and we reject based on that.



sort(allScores) #WE NEED TO SORT FIRST
length(allScores)*0.05 #This is the number of values that are above the P(0.05)
length(allScores) - length(allScores)*0.05 #this will give us the index of that, so its 18991.45 --> 18992

threshIndex = ceiling(0.95 * length(allScores)) #get the index of the upper C.I
threshIndex #Ceiling means it gives the next integer

sort(allScores)[threshIndex] # --> This will give us that number. (C.I upper)

quantile(x = allScores, probs = 0.95) #This is probably the previous index. This gives us an AUTOMATIC way

####### P-values #######


## >How p-values are distributed when null hypothesis is correct?
## Usually we assume that H0: μ1 = μ2 (Population). Basicly we could see that both the distributions of those are gaussian and are the same exactly. To construct the t-statistic we take into account the mean/average of the 2 values (healthy-sick) and variances of the 2 values. The statistic is a value at 0.05?

## P-value: ITS THE DENSITY OF AREA ABOVE A CERTAIN THRESHOLD (T(critical)). ΕΜΒΑΔΟΝ!
##It follows a uniiform distribution when null hypothesis i correct ( pvalue ~ U(0,1) | H0 )

##If we construct a t-statistic given that the H0 is true, the t-statistic has a normal   distribution. H1: E[healthy] > E[sick] (E = expression)
## Higher t-statistic (or lower for lower tail) means lower p-value, so more extreme to find that value/pvalue .. --> P(x>= t1 | H0) = P1 // P(p<=P1 | H0) = P1 // P(pvalue <= x | H0) = x

m = matrix(rnorm(100000, 100, 10), 10000, 10)
#(100K values, mean = 100, sd = 10) Distributed in 10000 rows, and 10 columns)
#Since we create a matrix from values that come from same distribution, the H0 here is correct
labels = rep(c("h", "d"), each = 5)
colnames(m) = labels

mean(m[1,1:5])
mean(m[1,6:10])
dim(m)

myttest = function(mat, lab1, lab2){
    pvalue = t.test(mat[lab1], mat[lab2])$p.value
    return(pvalue)
}

pvalues = apply(m, 1, myttest, 1:5, 6:10) #Apply on the matrix m, 1 = row (2 = column), the function myttest, and split based on columns 1-5, 6-10
hist(pvalues)
sum(pvalues < 0.05)


############ Correcting for multiple tests ################
## t-test for thousands of genes. For every such test, we set significance (0.05 usually), then after we apply the t-test we get a p-value
## IF: p-value < 0,05 --> REJECT probably H0 // p-value > 0,05 CANT REJECT H0
##FP : Test rejects H0 but H0 is true
##FP RATE: When we reject the null when we shouldnt reject it
## When the null is TRUE then the probability to find a p-value that is smaller than a is a.

## Multiple hypothesis testing :
##Family wise Error Rate (FWER) --> Probability to make atleast 1 type 1 error/1 False Positive (FP) (when we reject null when we shouldnt)
##If null is correct then prob = a . NOT = 1-a
## FWER : a' = 1 - (1-a)^N <= a*m
##Probability for N experiments is (1-a)^N <-------------------- No TYPE 1 Error in any of the N independent tests
##Atleast 1 error : 1 - (1-a)^N --> Ekthetiko plot increasing towards 1



##>>>>>>>#Control the FP when we perform multiple experiments #<<<<<<<<<<<<<<##
## To make up for this we use an a.pe = a/N ---> a.pe*N = (a/N)*N >= a'
sum(pvalues<0.05/10000) #We use 10.000 for 10.000 rows basicly (genes). We get 0

##If p.values < a/N --> pvalue*N < a _______________ THIS IS THE CORRECTED (ADJUSTED) PVAL.
##Since this gives values > 1, and that is a prob, then we just report back 1.
##p.adjust, "bonferroni" (means we multiply pvalue with N experiments

pbonf = p.adjust(pvalues, method = "bonferroni")

sort(pbonf, decreasing = FALSE)[1] #The smallest pvalue is ~0.18, and even that doesnt reject H0 which is what we want since it is correct

m1a = matrix(rnorm(50000, 100, 10), 10000, 5)
m1b = matrix(rnorm(50000, 120, 10), 10000, 5)
m1 = cbind(m1a,m1b)

pvalues2 = apply(m1, 1, myttest, 1:5, 6:10)
hist(pvalues2) #Shows alot with low pvalue. Different from adjust

which(pvalues2<0.05/10000)#Shows the index of values that escape the thresh. Differs
which(p.adjust(pvalues2, method = "bonferroni") < 0.05) #VERY CONSERVATIVE

sum(p.adjust(pvalues2)<0.05) #We should have rejected it for all of them (10K) but with bonferroni we only rejected 7 of them
sum(p.adjust(pvalues2, method = "fdr") < 0.05) # ~6K genes where we reject H0

####### FOR BONFERRONI ####### : The threshold for all the pvalues here is 0.05/m
##First we sort
sortedpvalues = sort(pvalues2, decreasing = FALSE)

plot(1:length(sortedpvalues),sortedpvalues, pch=19, cex=0.4, ylim=c(0,0.1)) #pch = fill point, cex = diameter (40% smaller). Yaxis: 0,0.1
abline(h=0.05/10000, col = "red") #THIS IS FOR BONFERRONI
abline(b=0.05/10000,a=0,col="green") #b = slope, a = intercept. THIS IS FOR FDR

plot(1:length(sortedpvalues),sortedpvalues, pch=19, cex=0.4, ylim=c(0,0.00001))
abline(h=0.05/10000, col = "red")

########### FDR CORRECTION ########### : The threshold for the pvalues is increasing per index. So it will be for for(i in 1:N){ i*(a/m) } ---> Better for biology
## Sort p-values

##------------------------------------------------------## REMEMBER TO ADD T.TEST

t.test(x, y = NULL,
       alternative = c("two.sided", "less", "greater"),
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95, …)




myttest = function(v, lab1, lab2, alter = "two.sided"){
    myt = t.test(v[lab1], v[lab2], alternative = alter)
    return(myt$p.value)
}

pvalues <- apply(b,1,myttest,1:6,7:15)

which(pvalues < 0.05) #dinei indexes

geneNames[pvalues < 0.05]

d = b[pvalues < 0.05,]
dim(d)


heatmap.2(as.matrix(d[1000,]), Rowv=F,Colv=T,dendrogram = 'column', trace = 'none', scale='row')


#katwfli = 0,05/30_000



##alternative = "greater" (enallaktikh oti to prwto groyp einai megalutero apo to deutero groyp) // "less" (to prwto group einai mikrotero apo ti deytero group) // "two.sided" h sxesh einai diaforetikh p.x m1#m2.
##prwta mpainei to group poy exw valei egw prwta. Panta paei me vash to prwto

##t.test --> t = timh statistic // th katanomh thn kserei h R, opote ypologizei to pvalue me vash oti oi pragmatikh diafora metaksu tis meses times anaferetai sto plhthusmo, dhladh oti den einai ish me mhden, kathw m1#m2, kai einai eite to ena. Me p = 0.27 den einai duskolo na vrw diaforetikous mesoys oroys pou na antistoixoyn se t = 1.155

##False positive: antras thetiko test egumosynhw #pval deixnei megisth empistosynh? gia to FP
##False negative : gunaika eguos me arnhtiko test





## p-value for two-sided : The t(critical) = the t where we have statistic significance changes. For two-sided we want the total embadon to be the 5%, so it is 2*p(one-sided). For one sided, we have a "lower" t-critical, so we have a bigger embadon for one side. For the two-sided the p-value is bigger (twice)--> The one-sided is half of the two-tail P-value



################## LINEAR MODELS ##################
##What it does, is try to fit a single line to explain the data as best as possible.
##Matrix (columns individuals, rows genes)
##Columns : Individuals
##Rows : Genes

##We try to find genes which have statistically difference average values between the two classes
## Plot : Find Average of expression (dot) and the height is the difference.
##We see in the plot that dots represent individuals for 1 gene and their expression.
##We do this for every gene.


## Expression = bo + b1*F (constant + coef. of factor * effect of factor ). We try to explain the expression values in terms of a factor (F). 1 Condition is Healthy (0)/Sick (1).  We try to find genes affected by the condition (sick) and so we want to explain the expression levels through the condition. If we find genes that fit well in the model, then those genes are those beeing affected by beeing healthy/sick
##We apply this process for EVERY GENE SEPERATELY not all genes together

##Since the line is straight, then 2 consecutive points will always differ by the same amount as the next ones and so on by Δ(expr). We cant fit the line to every type of data, meaning they cant be explained by a straight line (maybe can be explained by exponential)

##Let F: age ------> CONTINUOUS VARIABLE
## E = b0 + b1*age --> If age = 0, E = b0 (intercept - average prediction for the expression levels for age = 0. What this models predicts, is that when increasing age by 1 point, the expression levels change by an amount equal to b1 (slope). b1 > 0: (/), b1<0: (\)

lm(E~Age) #Explain expression levels as a function of the Age. Performs statistical evaluation on the significance of age on the expression of the values





#># DISCRETE EXPLANATORY VARIABLES (_FACTORS_) ##
## Healthy/Sick or Smoke/Dont-Smoke or Drud/No-Drug 

## MEANS MODEL ## --> Explain expression level through the below function
## E = b1*x1 + b2*x2 (x : factor with 2 levels - ex. healthy/sick). We say individual has Healthy: x1 = 1 (when healthy), x2 = 0 (not sick) --> E = b1*1 + b2*0 --> E = b1                   (this is the mean expression values of all the individuals is around b1)         Sick:  x1 = 0 (not healthy), x2 = 1 (when sick) --> E = b2

##WE CARE ABOUT B2-B1 --> Δ = b2 - b1 --> H0: Δ=0 // H1:Δ#0. What we test is if the difference is different to zero  (the statistical significance) (of their difference)

## Gene 1
healthy = c(2.4, 7, 2.3, 3.1)
sick = c(21, 22.3, 16, 8.2)

expression = c(healthy, sick) #Concatenate vectors
expression

condition = factor(c(rep(0,4), rep(1,4)), labels = c("healthy", "sick")) ## We set the level of healthy to 0 and the level of sick to 1, and we assign to 0 the healthy, and the sick to 1. If we put rep(1,4), rep(0,4) --> We get sick first then healthy. 0: Basic level
condition

mymod = lm(expression~condition + 0) #Since we have MEANS MODEL we need +0. 
summary(mymod)
##Coefficients:
##conditionhealthy: Estimate: b1 p-value assigned to b1
##conditionsick: Estimate: b2 p-value assigned to b2
##The 2 p-values dont tell us if healthy is different than sick, but if the b1 and b2 are statistically different than 0. p-value < 0.05 --> The b is stat.sign different than 0.
##This tells us only if the estimations are different (the 2 means of healthy vs sick) and not if the b1 and b2 are statistically different from each other

mean(healthy)
mean(sick)

##Overall should be used with contrast model

## MEANS REFERENCE MODEL ##

##E = b0+b1*x x = 0 (healthy) --> E = b0 // x = 1 (sick) --> E = b0+b1
## The difference is now is b0 + b1 - b0 = b1
## If we estimate b1 and it is statistically significant than 0, then the 2 condition expressions for that gene are statistically significant different, since b1 represents the difference between the 2 conditions

mymod_ref = lm( expression ~ condition)
summary(mymod_ref)
##Coefficients:
##Intercept: b0 --> Represents mean value of reference level which is healthy
##conditionsick: Shows how much higher the mean expression level of that condition is (b1). The mean expression for condition level is Intercept + This (b0 + b1). If we get p<0.05 then we reject the H0, and we say that for this data the expression level for the sick is stat.sign. different than of the healthy


### MORE THAN 1 EXPLANATORY VARIABLES ### (no interaction between the variables)
## Let gender + health condition

## E = b0 + b1x +b2g (gender)
##Reference level needs 1 condition for each of the two factors.
##We will set male = 0, healthy = 0.
##> E = b0 (Male Healthy individual)
##> E = b0 + b1 (Male Sick individual). If we estimate b1 != 0, then it means beeing sick affects expression level of males.
##> E = b0 + b2 (Female Healthy individual)
##> E = b0 + b1 + b2 (Female Sick individual). For females difference is b0+b1+b2-b0-b2=b1

## The b1 is the same in both cases of male/female meaning that beeing sick is the same for both genders --> NO INTERACTION BETWEEN age + health condition

##Since we are interested in the effect of one of the 2 conditions on the expression level, GIVEN that we have in our sample both female + males. We add the gender condition variable to our script


healthy = c(2.4, 7, 2.3, 3.1,   6.4, 7.4, 8.1, 6.7)
sick = c(21, 22.3, 16, 8.2,  25.2, 24.1, 23.2, 24.7)

expression = c(healthy, sick) 
expression


condition = factor(c(rep(0,8), rep(1,8)), labels = c("healthy", "sick"))
condition

gender = factor(c(rep(0,4), rep(1,4), rep(0,4), rep(1,4)), labels = c("male", "female"))
gender

mymod_ni = lm(expression ~condition + gender) # ___NO INTERACTION___ BETWEEN THE 2 (+)
summary(mymod_ni)

##When we say no interaction we mean that the effect on expression level for beeing sick is the same for either male or female. 

##Intercept: b0 (male healthy). Prediction of this differs from the mean of males. Because we dont assume any interraction between the 2 conditions. In order to assume no interraction the model has to shift some values to fit the line. This might mean there is interraction 
##The estimates only represent mean values if there is truly no interaction between the 2 variables.
##Beeing sick has an increase of b1 = 15.16  (regardless of gender) compared to intercept
##Beeing female regardless sick or not has an effect of b2 = 5.4. Both of these are s.s.
##This shows strong effect on expression levels

mymod = lm(expression~condition)
summary(mymod)
##This shows the mean value of healthy people (5.4) regardless of gender


##### Interaction #####

## E = b0 + b1*x + b2g + b3*x*g
## Male Healthy : E = b0
## Male Sick : E = b0 + b1 # bi : The effect of beeing sick when male
## Female Healthy : E = b0 + b2 : The effect of beeing female healthy
## Female Sick : E = b0 + b2 + b3 : The effect of beeing female sick. It is negative when the female decrease, since without it the model would predict an increase.

## If we assume the effect of beeing sick is same for both male/female then we wouldnt be able to see the decrease.

healthy = c(2.4,7,2.3,3.1,   36.4, 37.4, 38.1, 36.7)
sick = c(21,22.3,16,8.2,   1.4,3.2,4.5,4.4)

expression = c(healthy,sick)
expression

condition = factor(c(rep(0,8), rep(1,8)), labels = c("healthy", "sick"))
condition

gender = factor(c(rep(0,4), rep(1,4), rep(0,4), rep(1,4)), labels = c("male","female"))
gender

mymod_i = lm(expression ~ condition*gender) # ___WITH INTERACTION___ (*)
summary(mymod_i)

##b0 (intercept),b1 (consick),b2 (for healthy genderfem, you go from b0 to b0+b2),        b3 (we would go from b0 to b0 + b1 for sick, then b0+b1+b2 for female and for sick female b0+b1+b2+b3 thats why negative, consick:genderfem)

m1a = matrix(rnorm(25000, 100,10), 1000, 25) #50.000 gonidia, mesh timh 100, sd = 10(?), 1000 gonidia, 25 sthles
m1b = matrix(rnorm(25000, 120,10), 1000, 25) #
m2a = matrix(rnorm(25000, 120,10), 1000, 25)
m2b = matrix(rnorm(25000, 100,10), 1000, 25)

m = cbind(m1a,m1b,m2a,m2b) # 1000 gonidia (rows) x 100 individuals (cols)
m


condition = factor(c(rep(0,50), rep(1,50)), labels = c("healthy", "sick"))
condition

gender = factor(c(rep(0,25), rep(1,25), rep(0,25), rep(1,25)), labels = c("male","female"))
gender

mylm = lm(m[1,] ~ condition*gender) # For 1 gene
summary(mylm)

names(summary(mylm))

summary(mylm)$coefficients[2,4]

pvalues = vector(mode="numeric", length = nrow(m))
pvalues

for(i in 1:nrow(m)){
    mylm = lm(m[i,] ~ condition*gender)
    mysum = summary(mylm) #second row, 4th column for the effect of healthy sick
    pvalues[i] = mysum$coefficients[2,4]
}

pvalues
corr.pval = p.adjust(pvalues, method = 'fdr')
corr.pval










## Proodos ##


##H diaspora otan einai megalh den shmainei pws kai megaluterh diafora stoys mesous orous tha mas dwsei me mikrotero p-values. gia paradeigma


y1 = c(rnorm(50,100,10), rnorm(50,120,10))
t.test(y1[1:50], y1[51:100])$pvalue
mean(y1[1:50]) - mean(y1[51:100])


y1 = c(rnorm(50,100,50), rnorm(50,120,50))
t.test(y1[1:50], y1[51:100])$pvalue
mean(y1[1:50]) - mean(y1[51:100])


## microarray : pws einai ta dedomena: Grammes - Gonidia, Sthles - Atoma, 2 kathgories // Koitaw katanomh timwn prin kanw analysh. Vriskw log timwn an einai ptwsh. Me boxplot prepei na einai normal h katanomh. // Stoxos analyshs: // pvalues: ti symbolizei, H0,H1.
##multiple testing correction: Giati xreizetai, me ti tropoys to kanoyme, poios einai o poio austhros klp. // kanw t.test (me apply, ti symbolizei to ~ sto ttest. H1 polu shmantikh (an einai greater,less or twosided). Kalw synarthsh, kai px ti alternative tha valw
##Grammika montela: Suntelesths ths aneksarthths metablhths x ( pairnei times 0,1). O suntelesths mas deixnei poso allazei otan metaboyme apo to 0 sto 1.
##Ti symbainei an den laboume ypopsh to variance. Prepei na to balw wste na brw gonidia    pou h ekfrash toys epirazetai apo astheneia/ygeia alla einai pio eykolo na ta doyme otan baloume mesa kai to fulo.
##Ti shmainoun oi suntelestes sthn ekfrash. E= a1*ygeia + a2*fulo + a3*fulo:uygeia. To a1 deixnei pws h ugeia allazei thn efkrash otan to fulo einai 0. to a2 to antristrofo. to a3 exei nohma mono otan einai kai ta 2 1, dhladh mas leei posh prepei na einai h allhlepidrash gia na parw th mia apo tis dyo kathgories, mas deixnei to pragmatiko ypsos dhladh pou tha htan autes oi times pragmatika. San vash (0,0) pairnw to Female:Healthy.

##To a2 einai h diafora metaksu twn meswn timwn gia ena condition (px ygeia)
##Xwris to a3, tha phgainame pros ta katw kathws einai arnhtiko to a1 kai to a2. Epeidh omws pame pros ta panw shmainei oti uparxei kapoia allhlepidrash eksisoy kai to a3. To a3 != 0 shmainei pws uparxei ontws kapoia allhlepidrash. Prepei na einai statistika shmantika (p<0.05) diaforetiko tou 0

