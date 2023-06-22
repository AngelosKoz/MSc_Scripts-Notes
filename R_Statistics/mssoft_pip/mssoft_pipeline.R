## Line = Genome
##Column = Polymorphism


## 1 = Mutation (On the query when compared to the reference)
## 0 = No Mutation ( - // - )
##Removed columns with all(1) // all(0)

##library("readr") read_lines()

#obs <- read.table("ms_obs_final.out", sep = "", header = F)
#file_str <- paste(readLines("ms_obs_final.out"), collapse="\n")


## First we will read the files as strings.
obs <- readLines("ms_obs_final.out")
sim <- readLines("ms_sim_final.out")


##Now we will create a function that will take 3 parameters. 1) File, 2)Number of repeats and 3)The length of the original (observed) data, since the rest of the simulated are based on the original length of the observed, with the only difference beeing the polymorphic positions
getstats = function(file,reps=1,obslen=50){
    skipn = 1 #This will be used to skip the new-line character included when reading the files.
    k_pwd = vector() #Initialize 3 empty vectors for all our statistics.
    w_pp = vector()
    tajimad = vector()
    for(rep in 1:reps){ #We will use this loop to iterrate through the simulated data file
        setrep = file[(skipn+((obslen*skipn)-obslen)):((skipn+(obslen*skipn))-1)] #This is a correlation i found that works to keep one replicated (simulated) data at a time while skipping new line. This works only if the format is the same as simulated with new lines.
        genome = length(setrep) #This will be used for iterration in our for loops (for query)
        counter_genome = length(setrep) #While this will be used to iterrate through the rest of the rows (genomes)
        k_rep = 0 #Initialize the k value
        comparisons = 0 #This will basicly keep track of the comibnation (N 2).
        alpha1 = 0 #Initialize the alphas values
        alpha2 = 0
        for(i in 1:(genome-1)){ #Now the idea is to create a loop where we will use each of the genomes (rows) as a query and compare with the rest of the genomes (rows) in a pairwise mannerr. To do that we will go in an orderly manner and use as query each genome (row) and compare with the ones below that, to skip repeated calculations.
            query = as.numeric(strsplit(setrep[[i]], "")[[1]]) #With this we can load the string as a numeric vector and make the calculations. We assign each query (to compare) with each loop.
            alpha1 = alpha1 + (1/i) #Calculate the alphas
            alpha2 = alpha2 + (1/(i^2))
            for(j in 1:(counter_genome-1)){ #Now we iterrated through the rest of the rows to make pairwise comparisons.
                compare = as.numeric(strsplit(setrep[[i+j]], "")[[1]]) #With the i+j we can make sure we will not make comparisons with query itself and start the pairwise difference starting the row below query.
                k_rep = k_rep + sum((query - compare)==abs(1)) #calculate the numerator of k.
                comparisons = comparisons + 1 #Add 1 each time a pairwise difference is completed
            }
            counter_genome = counter_genome - 1 #With this we make sure that the second for loop (the rows we compare query with) are always 1 less to make up for the last comparison.
        }
        #Now this part is done after each replicate (or the original observed) is iterrated through. We calculate all of the statistics and values in one go and append them in the initialized empty vectors .
        s_temp = length(query)
        k_temp = k_rep/comparisons
        k_pwd = c(k_pwd, k_temp)
        w_temp = s_temp/alpha1
        w_pp = c(w_pp, w_temp)
        b1 = (obslen+1)/(3*(obslen-1))
        b2 = (2*( (obslen^2) + obslen + 3))/( (9*obslen)*(obslen-1))
        c1 = b1-(1/alpha1)
        c2 = b2 - (( obslen+2 )/(alpha1*obslen)) + (alpha2/(alpha1^2))
        e1 = c1/alpha1
        e2 = c2/((alpha1^2)+alpha2)
        taji_temp = (k_temp - w_temp) / (sqrt((e1*s_temp) + (e2*s_temp)*(s_temp-1)))
        tajimad = c(tajimad, taji_temp)
        skipn = skipn+1
    }
    statistics = rbind(k_pwd, w_pp, tajimad) #Row bind all of the statistics to keep at once
    return(statistics)
}


start_time <- Sys.time()
statsobs <- getstats(obs)
statssim <- getstats(sim, obslen=50, reps=10000)
statistics <- cbind(statsobs,statssim) #Add both observed and simulated data in a single matrix array
end_time <- Sys.time()
end_time - start_time


#Now we will do some visual modifications to our matrix array (change columns and rows names).
rnames <- c("k_pwd", "w_pp", "Tajima_D")
cnames <- sprintf("Sim_data%d",0:10000)
cnames[1] <- sprintf("Observed")
rownames(statistics) <- c("k_pwd", "w_pp", "Tajima_D")
colnames(statistics) <- c(cnames)

statdf <- t(statistics) #Transpose the matrix array so each column will be a statistic.
#head(statdf)

## ---------------------------------- ##
##Calculating the normalized counts for k statistic

k_mean <- mean(statdf[,"k_pwd"][-1])
k_var <- var(statdf[,"k_pwd"][-1])
k_sd <- sd(statdf[,"k_pwd"][-1])
statdf <- cbind(statdf, as.matrix(apply(statdf[,"k_pwd", drop = FALSE], 1, function(k) (k - k_mean)/k_sd))) #We calculate and add as a new column the normalized counts in one go
colnames(statdf)[4] <- "k_pwd_norm"

##Calculating the normalized counts for w statistic

w_mean <- mean(statdf[,"w_pp"][-1])
w_var <- var(statdf[,"w_pp"][-1]) #This is what was originally used in the script
w_sd <- sd(statdf[,"w_pp"][-1]) # <<------------
statdf <- cbind(statdf, as.matrix(apply(statdf[,"w_pp", drop = FALSE], 1, function(w) (w - w_mean)/w_sd)))
colnames(statdf)[5] <- "w_pp_norm"

##Calculating the normalized counts for Tajima's D statistic.

tajid_mean <- mean(statdf[,"Tajima_D"][-1])
tajid_var <- var(statdf[,"Tajima_D"][-1])
tajid_sd <- sd(statdf[,"Tajima_D"][-1])
statdf <- cbind(statdf, as.matrix(apply(statdf[,"Tajima_D", drop = FALSE], 1, function(td) (td - tajid_mean)/tajid_sd)))
colnames(statdf)[6] <- "Tajima_D_norm"

##ADD SQUARE ROOT

#Now we will make a function to calculate the euclidian distance.
eucld = function(st){
    obsnormal = st[1,4:6] #We want the normalized counts so we only use the last 3 columns, while keeping the first row only (obesrved)
    simnormal = st[-1,4:6] #The same logic, but here we remove the first row which corresponds to the observed
    dist = (t(apply(simnormal,1,function(x) obsnormal - x)))^2 #We will calculate the denominator that is inside the root with this for each statistic and use the power of 2 on those
    distances = apply(dist,1,function(x) sqrt(sum(x)))
    head(distances)
    return(distances)
}        
    

distances <- eucld(statdf)
smalldistances <- distances[head(order(distances),500)] #Now that we have calculated the distances, we use the order function, which sorts from low to high and keeps only the index. With head we can get the first 500 values from the output of order and use them as index.


smalldistances

##For our last step we need to load the parameters and keep the same indexes
params <- read.delim("pars_final.txt", header = F)
params <- params[,1]
parameters <- params[head(order(distances),500)]

meanpars <- mean(parameters) #Calculate mean
medianpars <- median(parameters) #Calculate median

meanpars
medianpars

plot(hist(parameters))

plot(density(parameters))



###
plot(hist(smalldistances))
plot(density(smalldistances))

