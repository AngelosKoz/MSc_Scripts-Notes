## >> Exercise 1 << ##

library(MASS)

x <- Pima.te[1:7] #lets rename the data for easier use, while keeping only quantitative columns

## Now all we have to do is run cor() on our data and will give back a correlation matrix.
## The matrix is square symmetrical and has its diagonal values equal to 1 since they are compared with theirselves (meaning directly proportional).

pearson_cor <- cor(x) #pearson is the default
pearson_cor

spearman_cor <- cor(x, method = c("spearman"))
spearman_cor

## We notice a few differences when looking at the pearson and spearman correlation matrix.Looking at our data we could say the we can notice a few outliers within the small data matrix, which usually affects Pearson correlation. Overall, we notice the only negative correlation is for npreg-bim but with a very low correlation coefficient. For the rest we notice slight positive correlation, with the highest beeign for bmi-skin and for age-npg (with about 0.66 for each). Noticable is also correlation for skin-glu, bmi-glu, bmi-bp, ped-glu, age-glu and age-bp with about 0.22-0.32.
##Since the default for our analysis is pearson, i cointinue my code using pearson.

## >Fischer transformation to get test statistic

## 1)With atanh
## Atanh will turn our matrix into z-score. This converts the data (correlation coefficients) to normal distribution data. The infinite is due to the fact the comparison is with itself. Info from: https://blogs.sas.com/content/iml/2017/09/20/fishers-transformation-correlation.html

z_score <- atanh(pearson_cor)
z_score


## 2)With psych (same results just mentioned it as information)
library(psych)
z_score <- fisherz(pearson_cor)
z_score




## Now that we have the z scores, we need to calculate the standard error and then the statistic. The standard error or z for fisher transformation is ~ 1/sqrt(N-3) where N is the size of the sample.

se_pima <- 1/sqrt(nrow(x)-3)
se_pima # se~0.0551

##Our null hypothesis here is H0: Correlation coefficients = 0 (means no correlation at all) and H1 != 0.
##To find the statistic, we need to divide the difference of our sample observations (our data) and the null hypothesis, whih in our case is that there is no correlation, thus cor.coef = 0.

statistic_fisher <- (z_score - 0)/se_pima
statistic_fisher # We notice that the statistic is about the same as the correlation matrix


## Since we now know the statistic we can calculate the pvalue based on the formula 2 * P(Z>=z | H0 is true). We just need to calculate the lower or upper tail and multiply by 2 as indicated in the formula. For correlation we calculate with 2 free variables, thus df = N-2, where N is the size of the sample. We can use the pt() function to calculate the upper/lower tail. From: https://bookdown.org/mslacour87/ch4_correlation_and_regression/Ch.4-Correlation-and-Regression.html

p_val <- 2*pt(abs(statistic_fisher), df = nrow(x)-2, lower.tail = F) #lower tail is false since we calculate one and multiply by 2.
p_val #We see that the diagonal of the matrix has pvalue = 0, which is normal since that is the maximum for the correlation coefficient. Means there is a direct proportional correlation. The same would be also be observed if we had -1 (inverse proportional).

dim(p_val)#7x7 matrix, so if we subtract the diagonal with is 7, we get a matrix of 42 entries, which half of them are identical, meaning we will have a total of 21 possible outcomes.

(sum(p_val<0.05)-7)/2 #In total we have 14/21 outcomes with a statistically significant pvalue with 95% confidence interval.

##And now to add the p in parenthesis, we convert to dataframe which are easier to work with.
pval_df <- as.data.frame(p_val) 
statistic_df <- as.data.frame(statistic_fisher)

lenp = dim(pval_df)[1]
for(i in 1:lenp){#we just need to loop over all the elements and copy contents from one to another
    for(j in 1:lenp){#it is a square matrix so we dont need to worry about diameters here
    statistic_df[i,j] <- paste( statistic_df[i,j]," (p=", pval_df[i,j],")", sep="")
    }
}

pearson_cor #correlation matrix

statistic_df #statistic and pvalue in parenthesis (as mentioned in exercise)


##Now lets do some plots.
## 1) Scatterplot for p-values - Correlation Coefficients

for(i in 1:lenp){
    png(paste(dimnames(p_val)[[1]][i],"_cor-pval_","scatter.png",sep = ""),width=700, height=400)
    plot(p_val[i,], pearson_cor[i,], main=paste(dimnames(p_val)[[1]][i],"as reference"),sep=" ", xlab="p-value", ylab="correlation")
    text(p_val[i,], pearson_cor[i,], row.names(p_val), cex=0.8,pos = 4, col="red")
    abline(v=0.05, col="red", lty="dashed") # Added a 0.05 line for taste
    dev.off()
}

## Comments: In general we notice as expected, that the higher (or lower, but we dont have any in this case) the correlation, the lower the p-value is as seen in the plots. That is expected since our H0 is that the correlation is equal to 0, so higher absolute deviation will lead to more statistically significant results (p-value). What we mean with higher correlation is that a change in one of the variables leads to change of the other variable. Higher correlation higher changes. The comments below concern 95% Confidence Interval, and i have denoted a red dotted line which shows the value for p = 0.05.
## i) Npreg vs Rest: Statistical significant with age (higher) and bp 
## ii) Glu vs Rest: Statistical significant with all except npreg
## iii) Bp vs Rest: Statistical significant with all except ped
## iv) Skin vs Rest: Statistical significant with all (higher correlation with bmi) except npreg+age
## v) B.M.I vs Rest: Statistical significant with all (higher correlation with skin) except age + npreg
## vi) Ped vs Rest: Statistical significant with all except npreg and bp
## vii) Age vs Rest: Statistical significant with all (higher correlation with nperg) except skin+bmi

##For the data.


lendata = length(x)
for(i in 1:lendata){
    for(j in 1:lendata){
        if (j <= i){
            next
        }
        xl = colnames(x)[i]
        yl = colnames(x)[j]
        png(paste(xl,"_",yl,"_scatter.png",sep = ""),width=700, height=400)
        plot(x[,i], x[,j], col = rep(1:2),main = paste(xl,"_",yl,"_scatter.png",sep = ""), xlab = xl, ylab = yl)
        dev.off()
    }
}

#We notice that the values for variables with statistically significant correlations have data that resemble each other, almost as if displaced upper or lower, more so for higher correlations.

## 2) Boxplot. For the boxplot i tried with and without the outliers (basicly the 1's that come with when correlating the variables with themselves, and as expected there was a noticable difference. To do so, i just replace the 1's with NA which is not taken into account at the boxplot.

fixed_cor <- as.data.frame(pearson_cor)
fixed_cor[fixed_cor==1] <- NA
fixed_cor

## Without removing 1's
png("cor_boxplot.png")
boxplot(pearson_cor)
dev.off()


## After removing 1's
png("fix_cor_boxplot.png")
boxplot(fixed_cor, ylim=c(-0.1,1))
dev.off()


## Comments: We notice that our boxes have bigger interquantile range (Q1-Q3), but the median does not change considerably. We notice that the outliers this time are the statistically significant variables that are correlated with the given variable. We also notice the negative correlation as mentioned before. Bigger interquantile range means more variable correlations with other variables, while the mean seems to be in most cases on the verge of the Q1 or Q3.

## 3) Histogram

png("cor_hist.png")
hist(pearson_cor)
dev.off()

png("fix_cor_hist.png")
hist(as.matrix(fixed_cor), xlim = c(-0.1,1))
dev.off()

## Comments: We can see the occurences of correlation values, most beeing below 0.4, but still at a statistically significant level

## 4) Correlation plots: http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
library(corrplot)

png("corrplot.png")

corrplot(pearson_cor, type="upper", order="hclust", 
         p.mat = p_val, sig.level = 0.05, insig = "blank") #blank leave non statistically significant correlations as blank


dev.off()

## 5) Heatmaps (same as above)
palette = colorRampPalette(c("green", "white", "red")) (20)

png("heatmap.png")

heatmap(x = pearson_cor, col = palette, symm = TRUE)

dev.off()

## Comments: Both the corrplot and the heatmap show almost identical things. For the corrplot we see only the values that are statistically significant when paired with one another. The deep blue circles are when variables are pitted with themselves thus leading to correlation equal to 1. On the side we see the graph of [-1,1] showing the color of correlation. This is a very nice and clear way to see visually the correlations. Same thing applies to the heatmap, where we can also see some type of clusters forming to the higly correlated variables. In this case high correlation is red.


##Using Hmisc library, we can get the correlation matrix (same as psych and atanh), which also calculates the p-value on its own. (An alternative way)
library(Hmisc)

hm_cor <- rcorr(as.matrix(x), type = c("pearson"))
hm_cor
hm_cor$r#([1]) #Correlation matrix

hm_cor$P#([3]) #P-value of correlations

corrplot(hm_cor$r, type="upper", order="hclust", 
         p.mat = hm_cor$P, sig.level = 0.05, insig = "blank")




## >> Exercise 2 << ##

## Lets recreate the steps we used on exercise 1 to make a function for it.
fisherpval = function(dataf,h0){ #added the h0 parameter to change it if we are not looking for correlation = 0
    cor_df = as.data.frame(cor(dataf))
    fz_score = atanh(cor_df)
    se_df = 1/sqrt(nrow(dataf)-3)
    statistic_df = as.matrix(fz_score - h0)/se_df #it needs as matrix to calculate pvalue
    p_val <- max(2*pt(abs(statistic_df), df = nrow(dataf)-2, lower.tail = F)) #We are using max since we know that there will be in total 4 values in our matrix, since are comparing 2 variables. Both pairs are identical, with one beeing 0 since it comes from correlation = 1. Thus since pval >= 0, the max will give us the pvalue of intereset
    return(p_val)
}

#Same function but stops at statistic (for use at the permutation filtering)
corstat = function(dataf,h0){
    cor_df = as.data.frame(cor(dataf))
    fz_score = atanh(cor_df)
    se_df = 1/sqrt(nrow(dataf)-3)
    statistic_df = as.matrix(fz_score - h0)/se_df
    return(min(statistic_df)) #We can use minimum since again we have pairs but this time the variables with themselves return +Infinite.
}


##I chose not to over-comment on this permutation again, since its the same core as the one from set 2. The change here is that we want to keep one of the variables the same, while changing the other.

corpermttest = function(n1,n2,perm=999,h0=0){   
    my_df = cbind(as.data.frame(n1), as.data.frame(n2)) #With this we ensure the data does not become a problem
    corr_stat_thresh = corstat(my_df,h0)
    all_stat = vector()
    
    for(i in 1:perm){
        pop_rand <- as.data.frame(n2[sample(nrow(n2)),]) #Lets shuffle the second variable
        new_pop = cbind(as.data.frame(n1),pop_rand) #And then bind n1 (unchanged) with the shuffled n2        
        new_stat = corstat(new_pop,h0)
        all_stat = c(all_stat, abs(new_stat))
    }
    return((sum(all_stat>=abs(corr_stat_thresh))+1)/(perm+1))
}

#Again we have the same core function as set 2 so i chose not to overcomment.
corpvalrepeats = function(n,h0=0){
    r.names <- sprintf("p-values")
    c.names <- sprintf("Repeat%d",1:1000)
    m.names <- c("p-value_with_asymptotic", "p-value_with_permuation_(999)")
    p_cor_fish <- vector()
    p_cor_perm <- vector()
    for(i in 1:1000){ 
        x <- as.data.frame(rbeta(n,3,4))
        y <- as.data.frame(rexp(n,1))
        temp <- cbind(x,y) #This is for fisher functions
        p_cor_fish <- c(p_cor_fish, fisherpval(temp,h0))
        p_cor_perm <- c(p_cor_perm, corpermttest(x,y))
        }
    allpval <- array(c(p_cor_fish,p_cor_perm), dim = c(1,1000,2), dimnames = list(r.names, c.names, m.names))
    for(j in 1:2){
        png(paste(dimnames(allpval)[[3]][j],"(n=",n,")","hist.png",sep = ""),width=600, height=350)
        hist(allpval[,,j], main = paste("Histogram of",dimnames(allpval)[[3]][j]), xlab = paste("p-value ( n =",n, ")"))
        dev.off()
    }
    return(allpval)
}


list_with_all <- list()
ns = c(10,15,20,30,40,50)
ind = 1
for(sample in ns){
    result = corpvalrepeats(sample)
    list_with_all[[ind]] <- result
    ind = ind+1
}


## With this we can see the proportion of each test with different n in a fancy and easy way.

for(arr in 1:6){ #We have 6 different populations, so we can access each one of them from the list like so.
    print(paste("For n=",ns[arr],sep=""))
    for(tp in 1:2){ #Since we have 2 different p-values we can access them like so
        if( tp == 1 ){
            print("1)The proportion for asymptotic pvalue is:")
        }
        if( tp == 2 ){
            print("2)The proportion for pvalue with permutation (999) is:")
        }
        print(sum((list_with_all[[arr]][,,tp]<0.05)/1000*100))
    }
}
                


## From the histograms we notice that with lower n, without permuation we have higher fluctuations of p-values, but as the n increases, we notice the tendency towards uniform distribution. As far as permutations go, we notice the uniform distribution between all chosen n. This means its more "stable" and is not affected by population number (n)
