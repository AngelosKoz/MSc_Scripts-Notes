##---------- Exercise 2 ---------##

getwd()

t.test(Theoph[,5],mu=15)

##We get a t-statistic = -40.228, p-value < 2.2e-16, 95% confidence interval: (4.466749, 5.454160) and an estimated mean of x = 4.960455

##Our null hypothesis is that the mean of x = 15 (Ho = 15) and alternative hypothesis that it is not equal to 15 (H1 != 15). We get a p-value that is extremely low and therefore with a pretty high level of confidence (although never completely sure) we can assume the null hypothesis is rejected and that the true mean is not equal to 15. In order to not reject the null hypothesis, and assuming that  p-value was 2.2e-16, then we would need an alpha (α) < 2.2e-16, since in order to reject the null hypothesis we need the p-value to be lower than our significance level.


##---------- Exercise 3 ---------##
library(MASS)
attach(birthwt)

t.test(chem,mu = 0)

## 95% Confidence interval: (2.043523, 6.517311)

t.test(chem, mu = 2.043523) # We notice that the p-value for this and 
t.test(chem, mu = 6.517311) # for this is exactly 0.05, which means for any value that inside the 95% confidence interval, we probably would not reject the null hypothesis. So any values below 2.043523 and above 6.517311 would lead us to probably reject the null while the reverse would lead us to probably not reject the null

#For example 
t.test(chem, mu = 2.1) #p-value = 0.05558 > 0.05 (alpha) --> We probably dont reject the null 
t.test(chem, mu = 6.4) #p-value = 0.0622 > 0.05 (alpha) --> We probably dont reject the null


##---------- Exercise 4 ---------##
library(MASS)
library(MKinfer) #This is required for the perm.t.test command, so i included the results just in case
attach(birthwt)


#Lets create a function to make permutation easier for later use. We wont need to write this again later on, and we can just call it with permttest(pop1,pop2,R)
permttest = function(n1,n2,perm=999){
   
    stat = t.test(n1,n2)$statistic #Lets set the threshold for our statistics later
    all_stat = vector() #We initialize an empty vector to append our t statistics
    l1 = length(n1) #find the first sample size to use later for index
    l2_ind = l1 + 1 #I had issue using this directly in indexing so i set a "dummy" variable to use later
    l2 = length(n2) #We do the same for the second sample 
    pop_ind = l1+l2 #We combine those to find the total length (we could also just use a length(all_pop) later 
    pop = c(n1,n2)

    for(i in 1:perm){ #the perm here has a default of 999 but it can be changed to any value, even though it might take longer with bigger repetitions.
        pop_random = sample(pop) #lets shuffle all our sample size
        new_pop1 = pop_random[1:l1] #we now choose the same number of values as the initial sample 1, which is given by the length(n1). In our case this would give us 1 through 12.
        new_pop2 = pop_random[l2_ind:pop_ind] #The same for here, now instead we have length(n1)+1 which is the next index, until the end. That means this will give values from 13 through 189 in our case (which is 177, the initial sample 2 length)
        #Now that we have finished randomizing our samples, we can start doing the t.test
        all_stat = c(all_stat, abs(t.test(new_pop1,new_pop2)$statistic)) #now we can run the randomized samples t.test. I chose to save the statistic but we can also save pvalues if needed as seen at the start of function. We use absolute here since we are interested for two-sided test
    }
    #return(all_stat)
    return((sum(all_stat>=abs(stat))+1)/(perm+1))#I was not sure if i should add the initial statistic we found, so i pressumed since we divide by permutation repeats + 1 that i should add it, so instead i added 1, since we use the sum which is a number and we know that our statistic is the threshold so we can use 1. We use absolute here since we are interested for two-sided test
}

pop1 <- bwt[ht==1] #First lets set our sample 1
pop2 <- bwt[ht==0] #and our sample 2


t.test(pop1,pop2)$p.val #With no permuation. 95% C.I  (-1024.4717,153.6751) and p-val = 0.1332


permttest(pop1,pop2,999)
permttest(pop1,pop2,9999)

##With my function i get the exact same results with MKinfer, so i think the explanation mentioned below suffices and it includes some comments about confidence intervals.


##The test below is done with MKinfer.

perm.t.test(bwt[ht==1],bwt[ht==0], R=999) #With 999 repeats, one of the results generated 
#(Monte-Carlo) permutation p-value = 0.1451 
#permutation difference of means (SE) = -443.2735 (213.2768) 
#95 percent (Monte-Carlo) permutation percentile confidence interval:
# -882.67140  -15.81631

perm.t.test(bwt[ht==1],bwt[ht==0], R=9999) #With 9999 repeats, one of the results generated
#(Monte-Carlo) permutation p-value = 0.1306 
#permutation difference of means (SE) = -435.9929 (213.378) 
#95 percent (Monte-Carlo) permutation percentile confidence interval:
# -863.05953  -20.95508

## The explanations below concerns the MKinfer results.

##>We slight fluctuations of our p-values when using a permuation based t-test with 999 repeats, which is not that great but not insignificant either. The biggest change concerns our 95% confidence intervals where we notice a big difference from the results of no permuation when compared to either permuation results. Both the confidence intervals when using permuations look alike which i would say it the most noticable difference when doing a simple t.test. Another noticable conclusion is that while the p-value increased with 999 permuations, it decreases again for 9999 permuations, looking more like our no-permuation t.test and even lower than that even if that is by a small fraction, which would probably mean the more repeats we do the more sure we are of our conclusions, or rather it would lead me to believe we get closer to a more "objective" truth(if we can call it objective). But with p-values > 0.05 we probably cannot reject the null hypothesis which in this case is that the true difference in means is equal to 0 (Ho = 0 | H1 != 0).
##Overall we get bigger noticable fluctuations with 999 permutations while with 9999 we get more consistent results almost identical to the initial p-value, which is normal.

perm.t.test(bwt[ht==1],bwt[ht==0], R=99999) # I did one more for 99_999 repeats and we notice that the p-value stays almost the same and the same can be concluded for the 95% confidence interval


##---------- Exercise 5 ---------##

library(MASS)
library(MKinfer)
attach(birthwt)

class(deaths)
length(deaths)
deaths


x1 = deaths[1:36]
x2 = deaths[37:72]



t.test(x1,x2) # 95% C.I (-66.28383,501.89494) and p-value = 0.1307, while the values for x ( deaths ) are about 2165 in those years and for y (deaths) are about 1947 in those years


##With my function
permttest(x1,x2, 999)
permttest(x1,x2, 9999)

##As mentioned in exercise 4, the comments and explanation done for MKinfer also apply to my function since we get same results, so i kept them since they also include comments about confidence intervals

##The test below is done with MKinfer.

perm.t.test(x1,x2, R=999) #With 999 repeats, one of the results
#(Monte-Carlo) permutation p-value = 0.1491
# permutation difference of means (SE) = 213.1272 (143.6269)
# 95 percent (Monte-Carlo) permutation percentile confidence interval:
#  -76.71667 510.57500

perm.t.test(x1,x2, R=9999) #With 9999 repeats, one of the results
#(Monte-Carlo) permutation p-value = 0.1273
# permutation difference of means (SE) = 218.3698 (143.7449)
# 95 percent (Monte-Carlo) permutation percentile confidence interval:
#  -61.16944 497.18611


## The explanations below concerns the MKinfer results.

##The results here are almost similar to exercise 4, meaning we notice slight changes in p-value, with the same pattern as exercise 4, meaning slight noticable fluctuations for 999 repeats and a slight decrease for 9999 repeats, reaching a lower p-value than the initial t.test, but almost identical to it. While all 3 have p-values > 0.05 we probably cannot reject the null hypothesis that the deaths for those years are exactly the same with 95% confidence (Ho: meanx1 - meanx2 = 0 | H1 meanx1-meanx2 != 0). The biggest difference here is that the confidence intervals are about the same when comparing to exercise 4 were we had bigger fluctuations in the lowerer-upper parts. Instead here the confidence intervals are about the same, leaning towards a more "loose" C.I for 999 repeats and a more strict C.I for 9999 (a bit stricter than no-permuation t.test). Overall, again, more consistent results with more permutations




##---------- Exercise 6 ---------##

set.seed(1234) #set the seed so we get the same random numbers when re-running. We have to run this when re-generating x or y.
##In this case set.seed is not required since we want to run different outputs for x and y

##Lets create a function to automate the exercise
pvalrepeats = function(n){
    ##Preparing array names
    ##size=n
    r.names <- sprintf("p-values")
    c.names <- sprintf("Repeat%d",1:1000)
    m.names <- c("p-value_with_equal_variance", "p-value_with_different_variance", "p-value_with_permuation_(999)")
    ##Initialize empty vectors
    p_equal_var <- vector()
    p_diff_var <- vector()
    p_perm <- vector()
    for(i in 1:1000){ #We loop 1000 times each time creating a new randomized x and y meaning different p-values.
        x <- rexp(n,1) ##Generate random numbers from an exponential distribution with a rate of 1
        y <- rnorm(n,1,5)##Generate y from a normal distribution with mean = 1 and variance = 25. Since variance is 25, the standard deviation is sqrt(25)=5
        p_equal_var <- c(p_equal_var, t.test(x,y, var.equal = T)$p.value)#We append each new p-value to the appropriate vector
        p_diff_var <- c(p_diff_var, t.test(x,y, var.equal = F)$p.value)
        p_perm <- c(p_perm, permttest(x,y,999))#Here i use my function for the permutations and not MKinfer, and since it generates just the p-value we can append it to our vector.
    }
    allpval <- array(c(p_equal_var,p_diff_var,p_perm), dim = c(1,1000,3), dimnames = list(r.names, c.names, m.names))
    for(j in 1:3){
        png(paste(dimnames(allpval)[[3]][j],"(n=",n,")","hist.png",sep = ""),width=600, height=350) #Saves in current directory
        hist(allpval[,,j], main = paste("Histogram of",dimnames(allpval)[[3]][j]), xlab = paste("p-value ( n =",n, ")"))
        dev.off()
    }
    return(allpval)
}


##Now that we have made a function to automate our results, lets run this for different n.

#We can also automated the above function for all the sample sizes, but i also include the manual way to do them with steps for checking proportion for all the 2 t.test and permutation. The histograms are also saved automatically through the function.
list_with_all <- list()
ns = c(10,20,30,40,50,60)
ind = 1
for(sample in ns){
    result = pvalrepeats(sample)
    list_with_all[[ind]] <- result
    ind = ind+1
}

##Saddly my function takes a longer time to run than MKinfer, but we notice same results and histograms.

sum((list_with_all[[1]][,,1]<0.05)/1000)*100 #This way we can access the arrays within the list and check our results. The list has 6 arrays which can be acces like [[1:6]] and the result of each different t.test like [,,1:3] (in both cases the : just symbolises the list of numbers). The higher the pops become, the less fluctutations between the 3 different t.test approaches become.



##In general we notice in the histogramms that we have a distribution resembling uniform, which is more so when the number of samples increases (n). So when comparing the n=10 with n=60 we notice that the p-values for n=60 are more uniformly distributed and have closer values than those compares to n=10

    
##Below we can see the results for different sizes manually, reporting back the histograms (saved at current working directory as seen above in the function) and the proportion of pvalues that are below our significance level (in this case we use alpha = 0.05)



                                        #N=10
result = pvalrepeats(10)
##checking the proportion of p-values that are statistically significant (for α = 0,05)
sum((result[,,1]<0.05)/1000)*100 # Value in %
sum((result[,,2]<0.05)/1000)*100 # Value in %
sum((result[,,3]<0.05)/1000)*100 # Value in %


                                        #N=20

result = pvalrepeats(20)
##checking the proportion of p-values that are statistically significant (for α = 0,05)
sum((result[,,1]<0.05)/1000)*100 # Value in %
sum((result[,,2]<0.05)/1000)*100 # Value in %
sum((result[,,3]<0.05)/1000)*100 # Value in %


                                        #N=30


result = pvalrepeats(30)
##checking the proportion of p-values that are statistically significant (for α = 0,05)
sum((result[,,1]<0.05)/1000)*100 # Value in %
sum((result[,,2]<0.05)/1000)*100 # Value in %
sum((result[,,3]<0.05)/1000)*100 # Value in %


                                        #N=40


result = pvalrepeats(40)
##checking the proportion of p-values that are statistically significant (for α = 0,05)
sum((result[,,1]<0.05)/1000)*100 # Value in %
sum((result[,,2]<0.05)/1000)*100 # Value in %
sum((result[,,3]<0.05)/1000)*100 # Value in %


                                        #N=50


result = pvalrepeats(50)
##checking the proportion of p-values that are statistically significant (for α = 0,05)
sum((result[,,1]<0.05)/1000)*100 # Value in %
sum((result[,,2]<0.05)/1000)*100 # Value in %
sum((result[,,3]<0.05)/1000)*100 # Value in %


                                        #N=60


result = pvalrepeats(60)
##checking the proportion of p-values that are statistically significant (for α = 0,05)
sum((result[,,1]<0.05)/1000)*100 # Value in %
sum((result[,,2]<0.05)/1000)*100 # Value in %
sum((result[,,3]<0.05)/1000)*100 # Value in %





