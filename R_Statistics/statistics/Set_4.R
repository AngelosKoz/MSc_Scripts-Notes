## Exercise 1 ##

## We notice that we have categorical variables (Trauma: Yes/No, Depression: Yes/No) which leads us to use the G^2 test of independence. Our null hypothesis here is that the 2 variables (Trauma and Depression) are independent. The alternative hypothesis is that they are not.

allm <- rbind(n1,n2)

g2ind = function(n1,n2){
    allm <- rbind(n1,n2) #Row bind our data in a matrix
    n1.sum = sum(allm[1,]) #Find sum of first row
    n2.sum = sum(allm[2,]) #Second row
    n.1sum = sum(allm[,1]) #Find sum of first column
    n.2sum = sum(allm[,2]) #Second column
    ntot = n1.sum + n2.sum #Find that total size
    #Now we calculate all the error for each position 
    exp11 = (n1.sum*n.1sum)/ntot #[1,1]
#    print(error11)
    exp12 = (n1.sum*n.2sum)/ntot #[1,2]
#    print(error12)
    exp21 = (n2.sum*n.1sum)/ntot #[2,1]
#    print(exp21)
    exp22 = (n2.sum*n.2sum)/ntot #[2,2]
#    print(exp22)
    g2 = 2*( (n1[1]*log(n1[1]/exp11)) + (n1[2]*log(n1[2]/exp12)) + (n2[1]*log(n2[1]/exp21)) + (n2[2]*log(n2[2]/exp22)) )
    return(g2)
}

g2ind(c(50,10),c(3,5)) #This is to test if i get the same results as the example of the lectures.


n1 <- c(251, 4)
n2 <- c(131, 33)
nmat <- rbind(n1,n2)

g2ind(n1,n2)

## Next, we need to check if our statistic has a value that rejects our null hypothesis or not. An approach for this is to compare our statistic with chi-square distribution at some degrees of freedom.

dfchi = (nrow(nmat)-1)*(ncol(nmat)-1)

#Our degrees of freedom are 1, so for alpha = 0.05, we can go to the chi-square distribution table and look for X1,0.05 which is ~ 3.84. Since our statistic is much larger than the critical value of chi-square (44.36 > 3.84), we can say that we probably reject the null hypothesis that Trauma and Depression are independent meaning there probably is some sort of relationship between them.


## Exercise 2 ##
                 
library(MASS)

x <- Pima.te[,1:7]


## To do regression models we can use the lm() function (linear model). Since we want to compare 1 variable to the rest, we should think what approach to use. We will use the (+) inbetween the variables ( example1 ~ example2 + example3 +...) which considers that the variables have no interaction between them, otherwise we need to use the (*) ( example1 ~ example2 * example3 *...) which considers the variables interact with each other.


##Considering we take only the first 7 columns of our data, we can write the function as  E = b0*var0 + b1*var1 + b2*var2 + b3*var3 + b4*var4 + b5*var5 + b6*var6                   where we take for example the first variable (npreg, which is var0 and will always be 1) as response variable and compare it to the rest. The same model with different variable placements (and ofcourse different b) is used for the rest of the comparissons. In each case the b0 will be the Intercept, meaning when all other values are 0, that is where we expec the mean value of our response variable to be. All the variables take 0 or 1 values, so for example if var3 = 0, that would mean how my model is explained when we take out the variable 3, while the b are scalars that show the "impact" of that variable on the mean (whether it has a negative or positive value, or how big or small it is)

## The adjusted coefficient of determination (R-square) shows how well the data is explained by the linear model (how well the model is adjusted), or how much of the variance in the response variable is accounted for by the predictor variables (how close the data are)


## We will use an alpha = 0.05.


##1) Npreg as response ##

npreglm = lm( x[,1] ~., x[,2:7] ) #This is basicly what it does: lm(Pima.te[,1] ~ Pima.te[,2]+Pima.te[,3]+Pima.te[,4]+Pima.te[,5]+Pima.te[,6]+Pima.te[,7])
summary(npreglm)

## From the summary, we can see that the only statistically significant variable that seems to have an effect on npreg is age and for each age increase we expect npreg to increase bhy 0.2. A small but positive effect. The R-square adjusted seems to axplain a 44.64% of the variance.


npregres = resid(npreglm) # This will save the residuals to a variable for later tests (normality, homoscedasticity.


## Normality test
hist(npregres)

qqnorm(npregres) # We notice that alot of the residuals seem to follow normal, but then the shapiro test shows otherwise.

##Shapiro test tests normality with null hypothesis beeing H0: our data follow normal distribution.
shapiro.test(npregres) # p < 0.05, meaning we probably reject the null hypothesis that our data comes from normal distribution.


## Homoscedasticity check
plot(npregres) #We do not notice some kind of pattern

## I thought in case of finding some statistically significant variable, to do an anova test with the reduced model beeing that (or those) variables and see if indeed by adding that variable we get a better prediction model. The null hypothesis is that the effect of all the predictor variables is 0 (b1,b2,..,b6 = 0). Alternative that at least one is different than 0.
## F-test (anova)

modelnpreg = lm(npreg ~. , data = Pima.te[,1:7])
reducednpreg = lm(npreg ~. -age , data = Pima.te[,1:7])

anova(reducednpreg, modelnpreg) #pval < 0.05 meaning we probably reject the null hypothesis that non of the predictor variables have an effect on npreg, and since it was age that we removed, age has an impact indeed and we get a better model when included.




##2) Glu as response ##
glulm = lm( x[,2] ~., x[,c(1,3:7)])
summary(glulm)

## We notice that bmi,ped and age have a statistically signifiant effect on glu concentration, with ped having a very high positive effect, meaning for each ped glu increase by 14.5. The problem here is that our R-squared is very low and explains only 14,68% of the variance. This means probably other variables we have to take into account, or just not a linearly explainable model?

glures = resid(glulm)

## Normality test
hist(glures)

qqnorm(glures) #(same comments as before)

shapiro.test(glures) # p < 0.05 (same comments as before)

## Homoscedasticity check
plot(glures) #(same comments as before)

## F test

modelglu = lm(glu ~. , data = Pima.te[,1:7])
reducedglu = lm(glu ~. -bmi -ped -age, data = Pima.te[,1:7])

anova(reducedglu, modelglu) #Same comments as previous anova. This time we removed 3 variables. Also, doing a reduced model by removing each one individually, gives us a higher but still statistically significant p-value, except for bmi (non statistically significant).


##3) Bp as response ##
bplm = lm( x[,3] ~., x[,c(1:2,4:7)] )
summary(bplm)

##It seems bmi and age have a statistically signifiant effect on bp, meaning for increased bmi, the bp (bloodpressure) will increase by 0.6 and for each age by 0.4. Again we notice a very low percent of the variance beeing explained (R-squared adjusted = 0.2041)

bpres = resid(bplm)

## Normality test
hist(bpres) 

qqnorm(bpres) #(same comments as before)

shapiro.test(bpres) #(same comments as before)

## Homoscedasticity check
plot(bpres) #(same comments as before)


## F test

modelbp = lm(bp ~. , data = Pima.te[,1:7])
reducedbp = lm(bp ~. -bmi -age, data = Pima.te[,1:7])

anova(reducedbp, modelbp) #Again, removing bmi and age seems to have a statistically significant effect on the model.



##4) Skin as response ##
skinlm = lm( x[,4] ~., x[,c(1:3,5:7)] )
summary(skinlm)

##The only variable we notice that has a statistically significant effect on skin, is bmi and for each bmi increase, skin increases by 0.88 and with a decent percent of the variance explained (R-squared adjusted = 0.438)

skinres = resid(skinlm)

## Normality test
hist(skinres) 

qqnorm(skinres) 

shapiro.test(skinres) # p-value = 0.1 > 0.05. This time it seems that the residuals probably come from a normal distribution.

## Homoscedasticity check
plot(skinres) #(same comments as before)

## F test

modelskin = lm(skin ~. , data = Pima.te[,1:7])
reducedskin = lm(skin ~. -bmi, data = Pima.te[,1:7])

anova(reducedskin, modelskin) # Same comments as npreg (the first one)



##5) BMI as response ##
bmilm = lm( x[,5] ~., x[,c(1:4,6:7)] )
summary(bmilm)

##Bp, skin and glu seem to have a statistically significant effect on bmi, with glu having a small effect of 0.02 increase of bmi per glu, while each increase of bp increases bmi by 0.13 and skin by 0.44. So far the highest explained variance of the data, close to 50% (R-squared - 0.49).

bmires = resid(bmilm)

## Normality test
hist(bmires)

qqnorm(bmires)

shapiro.test(bmires) #(p<0.05, not from normal distribution)

## Homoscedasticity check
plot(bmires) #(same comments as before)

## F test

modelbmi = lm(bmi ~. , data = Pima.te[,1:7])
reducedbmi = lm(bmi ~. -glu -bp -skin, data = Pima.te[,1:7])

anova(reducedbmi, modelbmi) #Same comments as before.



##6) Ped as response ##
pedlm = lm( x[,6] ~., x[,c(1:5,7)] )
summary(pedlm)

## The only statistically significant variable that seems to effect ped is glu by a very small amount (0.002 increase of ped per glu) and at same time it seems that ped cannot be explained by a linear model since the adjusted R-squared = 0.06 (0.6% explained).

pedres = resid(pedlm)

## Normality test
hist(pedres)

qqnorm(pedres) #We can clearly see that the data are curved and dont fit a straight line

shapiro.test(pedres) #(p<0.05, probably not from normal distribution)

## Homoscedasticity check
plot(pedres) #(same comments as before)


## F test

modelped = lm(ped ~. , data = Pima.te[,1:7])
reducedped = lm(ped ~. -glu, data = Pima.te[,1:7])

anova(reducedped, modelped)#A very high compared to the rest p-value, but still statistically significant, meaning glu probably has some kind of effect on ped




##7) Age as response ##
agelm = lm( x[,7] ~., x[,c(1:6)] )
summary(agelm)

##As expected, we notice npreg,glu and bp having a positive statistically significant effect on age. Meaning that for each npreg increase, age incrases by 2. For glu 0.04 and for bp 0.17. All in all higher values from npreg glu and bp are correlated with higher age (probably)

ageres = resid(agelm)

## Normality test
hist(ageres)

qqnorm(ageres) #We can clearly see a curve

shapiro.test(ageres) #The residuals do not follow a normal distribution (p<0.05)

## Homoscedasticity check
plot(ageres) #(same comments as before)

## F test

modelage = lm(age ~. , data = Pima.te[,1:7])
reducedage = lm(age ~. -npreg -glu -bp, data = Pima.te[,1:7])
anova(reducedage, modelage)

reducedage = lm(age ~. -npreg, data = Pima.te[,1:7])
anova(reducedage, modelage)

reducedage = lm(age ~. -glu, data = Pima.te[,1:7])
anova(reducedage, modelage)

reducedage = lm(age ~. -bp, data = Pima.te[,1:7])
anova(reducedage, modelage)

## All of these variables seem to have a statistically significant effect on the prediction of age whether individually removed or alltogether, with npreg having the strongest effect (which makes sense since npreg = number of pregnancies)


## Overall and as expected we notice same type of correlations as we saw in the previous set of exercises.
