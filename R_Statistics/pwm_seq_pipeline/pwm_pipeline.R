getwd()
library(seqLogo) 
library(gplots)
library(stringdist)


# ----------- Bash data manipulations ------------ #

#sort all_genes.txt -R | head -n 1000 > random_1000.genes #Sort all the genes randomly and pick the first 1000

#sed -i '1s/^/Gene_Stable_ID\n/' random_1000.genes #Insert a first row header
#sed 's/ /_/g' gluconeogenesis_gname.txt #Remove blank characters
#sed 's/ /_/g' glycolysis_gname.txt      #--------- // ----------




#--------------- Question 1 --------------#
##An aproach to find commons with bash##
#cat glycolysis_gname.txt gluconeogenesis_gname.txt > gly_glu_gnames.txt # merge files together
#sed -i '/Gene\ stable\ ID/d' gly_glu_gnames.txt #remove the gene stable id
#sort gly_glu_gnames.txt | uniq -d #print only diplicates lines


#Before we even test this, we expect to find no commons when taking genes in the class of tens/dozens from a pool of 69_000 genes. While doing preliminary test with bash showed 14 commons between the 2 pathways. Now lets check this with phyper.

gly67 <- read.delim("glycolysis_gname.txt", header=T)
colnames(gly67)[1] <- "Gene_ID"

glu35 <- read.delim("gluconeogenesis_gname.txt", header=T)
colnames(glu35)[1] <- "Gene_ID"

pathways_common  <- merge(gly67, glu35, by = "Gene_ID")

nrow(pathways_common)  # The common genes: 14

#So in order to check the probability of getting common genes, we could aproach this with dhyper which would give us the probability of getting from 0 to N (where N would be the total genes). We could check all the probabilities from 0 to common as so:

q_1  <- phyper(nrow(pathways_common)-1, 35, 69299-35, 67, lower.tail=FALSE)
q_1

##As expected the probability (pval) is very low since we have a huge pool and very small samples. This means our initial hypothesis that we expect to find more than 14 is rejected, which is also why we reduce the common by 1 (common - 1). This means that we have significant enrichment (in this case commons) meaning there is a very low probability of seeing equal to or bigger overlaps (commons)
##It does't matter if we swap 67 with 35, we get the same probability.


# ------------------ # Starting script to get questions 2/3 # -------------- #


#background = 1000 random genes from 69_299
#foreground = my selected path(s) genes

# Make a function to read the data i will use (both foreground and background) into a specific form (here we use list)
concatenate_data = function(filepath) {
  con = file(filepath, "r")
  data = list()
  while ( TRUE ) {
    line = readLines(con, n = 1) #read the text line by line (n=1)
    if ( length(line) == 0 ) { #If we reach the last line of the file we break from the infinited loop (since while true runs forever)
      break
    }
    if(length(grep(">", line)) > 0){ #we check lines for the ">" which indicates the gene name, and the line below it are the gene sequence. In this case we dont take the gene itself we take 500 nucleotides upstream from the gene.
        name = gsub("^>(\\w+)", replacement="\\1", x=line) #Since the ensembl gene name starts with > we replace it to keep the name only. this would be helpful if we had even more information we did not need
        data[[name]] = "" #we store the above-mentioned ENS name in the list as empty in order to add the sequences that follow it in the lines below, until the next name pops up. Here we also stored the sequence that is 1 line below our ENS ID
    }
    else{
      #line = gsub("\r?\n|\r", "", line)
      ##print(line)
      data[[name]] = paste(data[[name]], line, sep="") #As mentioned here we store all the lines that have sequence in order to end up with a single concatenated sequence in one line
    }
  }
  close(con) #closes the opened file
  return(data) #and returns my data
}


#---------------------------------------------#


kmer=8 #set the desired motif length

#Making a function to find and count all of the substrings

processSequences = function(seqList, len=8){
  x = sapply(seqList, function(i){ #Sapply works like a loop, applying a kind of function to all the elements (i - is a single element) of a data type. Here we work on a list.
    v = strsplit(i,"")[[1]] #We split the values (sequence) of each element (i which is gene name) into single nucleotides using nothing ("") as seperator. v is a vector so we use the [[1]] to get the single nucleotides stored there
    sapply(1:(length(v)-len+1), function(j){paste(v[j:(j+len-1)], collapse="")}) #Here, for each nucletodide (element of the vector) we take the whole length of the vector and remove the kmer length and add 1. Basicly this is seq - motif +1 which is all the possible places for a given kmer/motif, and with the sapply we take a step at a time finding those motifs.
  })
  table(x)#this returns all the possible substrings and counts the substrings so we know how many of each k-mer exist (in our case it will count 8-mers). They are returned as vectors with the name beeing the morif and the value is the number of occurences (counts)
}


#---------------------------------------------#


# Make a function to find p-value based on hypergeometric distribution of most common 8-mers Foreground vs Background. Basicly we want to see if we get a motif that is over-represented in the foreground when compared to the background.

getProb = function(foreground, background){ #create an empty vector to add the pvalues later
  probs = vector("numeric", length=length(foreground))
  sumforground = sum(foreground)
  sumbackground = sum(background)
  for(i in 1:length(foreground)){#For every foreground motif, it will test if it exists in the motifs of the background
    bcounts = 0
    if( names(foreground)[i] %in% names(background)){
      bcounts = background[[names(foreground)[i]]]#if it exists, it will save the background counts of that motif
      prob = phyper(q=foreground[i]-1,  m = bcounts, n = sumbackground - bcounts, k = sumforground, lower.tail = FALSE)#and here it will calculate the probability, since we have the foreground motif counts. (m = backgrounds counts (white) // n = size of background - motif counts (black) // k = total motif counts in foreground). We use False since we want counts that are equal to, or greater tha background.
      probs[i] = prob
    }else{#This contemplates for a motif not existing in the background, which is why we set it to 1
      bcounts = 1
      prob = phyper(q=foreground[i]-1,  m = bcounts, n = sumbackground - bcounts, k = sumforground, lower.tail = FALSE)
      probs[i] = prob
    }
  }
  names(probs) = names(foreground)
  return(sort(probs, decreasing=FALSE))#sort from small to great (since we care for lowest p-values) Lowest means it is the most over-represented in the foreground when compared to background
}

#---------------------------------------------#


##For the random genes (background dataset)## (we dont need the sample() from tha part of the code since we have already selected 1000 of our own)

random_genes_data = concatenate_data("random_genes.fa") # class() = list # random_genes_data$ENSG00000291294 

background_data =  processSequences(random_genes_data, kmer) #len = kmer

#names(background_data) # this shows all the k-mers # 61_225 different k-mers
#background_Counts = sum(background_data) #this counts the occurences of all the k-mers
#background_Counts #493_000 occurences of k-mers



##For glycolysis (foreground dataset)##

glycolysis_data = concatenate_data("glycolysis_67.fa")

foreground_glycolysis = processSequences(glycolysis_data, kmer)

#names(foreground_glycolysis) # 21_595 different 8-mers
#glycolysis_Counts = sum(foreground_glycolysis) # 33031



## Hypergeometric, find p value (glycolysis)

hyperprobs_glycolysis = getProb(foreground_glycolysis, background_data)

#save(hyperprobs_glycolysis, file="hyperprobs_glycolysis.RData")

#sort(foreground_glycolysis, decreasing=TRUE)[1] #GGGGCGGG most common with 31 counts





##For gluconeogenesis (foreground dataset)##

gluconeogenesis_data = concatenate_data("gluconeogenesis_35.fa")

foreground_gluconeogenesis = processSequences(gluconeogenesis_data, kmer)

#names(foreground_gluconeogenesis) #12_727 different 8-mers
#gluconeogenesis_Counts = sum(foreground_gluconeogenesis) #17_255 occurences of the 8-mers


##Hypergeometric, find p value (gluconeogenesis)

hyperprobs_gluconeogenesis = getProb(foreground_gluconeogenesis, background_data)

#save(hyperprobs_gluconeogenesis, file="hyperprobs_gluconeogenesis.RData")

#sort(foreground_gluconeogenesis, decreasing=TRUE)[1] #AAAAAAAA most common with 26 counts




# ---------- Question 2 --------- #

filter_glycolysis  <- hyperprobs_glycolysis[hyperprobs_glycolysis < 0.001]
length(filter_glycolysis)# 468 substrings from glycolysis with pval < 0,001


filter_gluconeogenesis  <- hyperprobs_gluconeogenesis[hyperprobs_gluconeogenesis < 0.001]
length(filter_gluconeogenesis)# 171 substrings from gluconeogenesis with pval < 0,001


common_filtered <- filter_glycolysis[names(filter_glycolysis) %in% names(filter_gluconeogenesis)] #here we use the %in% and we use the names of the vector to check for common motifs
common_filtered
length(common_filtered) #23 common motifs. 


# -------------------------------------------------------- #


###Building the PWM ###


## Make a function to find similar motifs (hamming distance <= 2) to the most overrepresented string found from hyper (based on pvalue). That would be for example hyper_result[1]
##We use the library stringdist for this

getAllInstances = function(candidate, foreground, threshold){#Candidate = lowest p.val motif, foreground = all of the motifs, threshold = 3 because we use smaller than below). We can try to set this to <= 2 instead??
  allnames = names(foreground)#All the motif names (not the counts)
  motifstrings = c()
  for(i in 1:length(allnames)){
    if( stringdist(candidate, allnames[i], method = "hamming") < threshold){#we compare the candidate motif with all the rest of the motifs based on hamming distance. 
      motifstrings = c(motifstrings, rep(allnames[i], foreground[i]))#If we get a match based on hamming distance it will take that specific motif (i) and find it's counts in  the foreground
    }
  }
  return(motifstrings)
}


#---------------------------------------------#


## Make a function to find the occurences of all 4 nucleotides to use it in the PWM later

getACGT = function(dataset, alphabet=c("A", "C", "G", "T")){
  counts = vector("numeric", length=length(alphabet))
  names(counts) = alphabet
  for(i in 1:length(names(dataset))){
    nam = strsplit(names(dataset)[i], "")[[1]]
    ntall = rep(nam, each=dataset[i])
    ntcounts = table(factor(ntall, levels=alphabet))
    counts[names(ntcounts)] = counts[names(ntcounts)] + ntcounts 
  }
  counts/sum(counts)
}

basefreqs = getACGT(background_data)
basefreqs #A=0.2559671 | C=0.2487546 | G=0.2475139 | T=0.2477644

#---------------------------------------------#

## Function to build the PWM 


getPWM = function(stringMotifs, length=6, alphabet =c("A", "C", "G", "T"),  freqs = rep(0.25, 4)){
  pfm = matrix(0, nrow=4, ncol=length)
  row.names(pfm) = alphabet
  for(i in 1:length(stringMotifs)){
    v = unlist(strsplit(stringMotifs[i], ""))
    for(j in 1:length(v)){
      pfm[v[j], j] = pfm[v[j], j] + 1
    }
  }
  ppm = pfm/colSums(pfm)
  pwm = pwm = log2((ppm+1e-4)/freqs)
  return(list(pwm=pwm, ppm=ppm))
}

#---------------------------------------------#


                                        # For glycolysis

names(hyperprobs_glycolysis)[1] #top candidate is AAACCCGG

hyperprobs_glycolysis <- sort(hyperprobs_glycolysis,decreasing=FALSE)#Sort so we can get the first since we have ties.


motifs_glycolysis = getAllInstances(names(hyperprobs_glycolysis)[1], foreground_glycolysis, 3) #find similar motifs to the top candidate ( names(hyperprobs_glycolysis)[1] ). We could use a different number other than 3 (for example 2.5). 3 is more flexible i supposed

ps_glycolysis = getPWM(stringMotifs = motifs_glycolysis, length=kmer, freqs = basefreqs)#Build the position weight matrix and visualize it with seqlogo using the matrix of frequencies (ppm) which is the most common nucleotide represented in each position (sum of each position should be 1)

png(file="~/Desktop/Master/R_masterclass/scripts/seqlogo_glycolysis.png",width=600, height=350)

seqLogo(ps_glycolysis$ppm)

dev.off()



                                        # For gluconeogenesis

names(hyperprobs_gluconeogenesis)[1] #top candidate is AAGTCGGG

hyperprobs_gluconeogenesis <- sort(hyperprobs_gluconeogenesis,decreasing=FALSE)#Sort so we can get the first since we have ties.


motifs_gluconeogenesis = getAllInstances(names(hyperprobs_gluconeogenesis)[1], foreground_gluconeogenesis, 3)

ps_gluconeogenesis = getPWM(stringMotifs = motifs_gluconeogenesis, length=kmer, freqs = basefreqs)

png(file="~/Desktop/Master/R_masterclass/scripts/seqlogo_gluconeogenesis.png",width=600, height=350)

seqLogo(ps_gluconeogenesis$ppm)

dev.off()

#---------------------------------------------#


## Scan the sequences to get the maximum locations

#---------------------------------------------#

# Make a function to scan sequences (8mers) and get a score based on the $pwm which shows each position score (basicly how often we see a specific nucleotide), which can be seen in the seqlogo
getScoreSimple  = function(vstring, pwm){
  score = 0
  v = strsplit(vstring, "")[[1]]
  scores = vector("numeric", length=length(v)-ncol(pwm)+1)
  for(i in 1:(length(v)-ncol(pwm)+1)){
    score = 0
    for(j in 1:ncol(pwm)){
      letter = v[i+j-1]
      score = score + pwm[letter, j] 
    }
    scores[i] = score
  }
  return(scores)
}



#---------------------------------------------#

##Making a heatmap for glycolysis

res_glycolysis = t(sapply(glycolysis_data, getScoreSimple, ps_glycolysis$pwm)) #we apply the function to the original sequence file (after we have concatenated the sequences for each gene) based on the $pwm and we transpode the data so we have rows with the genes (67) and columns the sequence (here we get 493 columns)


logicalres_glycolysis_strict = matrix(as.numeric(res_glycolysis > 3), nrow=nrow(res_glycolysis))#we keep scores that sum above 3, we expect less results in the heatmap

logicalres_glycolysis_soft = matrix(as.numeric(res_glycolysis > -4), nrow=nrow(res_glycolysis))#and scores that sum above -4, we expect more results in the heatmap





png(file="~/Desktop/Master/R_masterclass/scripts/heat_glycolysis_cutoff(3).png",width=600, height=350)

heatmap.2(logicalres_glycolysis_strict, dendrogram='none', Rowv=F, Colv=F, trace = 'none')

dev.off()



png(file="~/Desktop/Master/R_masterclass/scripts/heat_glycolysis_cutoff(-4).png",width=600, height=350)

heatmap.2(logicalres_glycolysis_soft, dendrogram='none', Rowv=F, Colv=F, trace = 'none')

dev.off()




#---------------------------------------------#


##Making a heatmap for gluconeogenesis

res_gluconeogenesis = t(sapply(gluconeogenesis_data, getScoreSimple, ps_gluconeogenesis$pwm))

logicalres_gluconeogenesis_strict = matrix(as.numeric(res_gluconeogenesis > 3), nrow=nrow(res_gluconeogenesis))

logicalres_gluconeogenesis_soft = matrix(as.numeric(res_gluconeogenesis > -4), nrow=nrow(res_gluconeogenesis))




png(file="~/Desktop/Master/R_masterclass/scripts/heat_gluco_cutoff(3).png",width=600, height=350)

heatmap.2(logicalres_gluconeogenesis_strict, dendrogram='none', Rowv=F, Colv=F, trace = 'none')

dev.off()




png(file="~/Desktop/Master/R_masterclass/scripts/heat_gluco_cutoff(-4).png",width=600, height=350)

heatmap.2(logicalres_gluconeogenesis_soft, dendrogram='none', Rowv=F, Colv=F, trace = 'none')

dev.off()


## --------------- Question 3 --------------- #
##Q3.2 For tasks 2 and 3 we simply find the maximum score for each row (gene) with a max()
glycolysis_maxscores <- data.frame(apply(res_glycolysis, 1, max))
names(glycolysis_maxscores)[1] <- "max_motif_score_gly"
glycolysis_maxscores


##Q3.3
gluconeogenesis_maxscores <- data.frame(apply(res_gluconeogenesis, 1, max))
names(gluconeogenesis_maxscores)[1] <- "max_motif_score_glu"
gluconeogenesis_maxscores


##Q3.4 For task 4 and 5 we will use the inverted data with pwm as seen below, and to check the differences we could do a simple subtraction to see how far off the 2 scores are. I did an extra step and merged those differences to check how they compare between the common genes

res_gluconeogenesis_by_glycolysis = t(sapply(gluconeogenesis_data, getScoreSimple, ps_glycolysis$pwm))


gluconeogenesis_by_glycolysis_maxscores <- data.frame(apply(res_gluconeogenesis_by_glycolysis, 1, max))
names(gluconeogenesis_by_glycolysis_maxscores)[1] <- "max_motif_score_glu_gly"
gluconeogenesis_by_glycolysis_maxscores

q3_4 <- cbind(gluconeogenesis_maxscores, gluconeogenesis_by_glycolysis_maxscores)
q3_4$diff1 <- q3_4[,1] - q3_4[,2]
q3_4$genes <- rownames(q3_4)
colnames(q3_4)
q3_4_sort <- q3_4[order(q3_4$diff1, decreasing = TRUE),]
q3_4_sort#We dont see great difference beside the top and bottom score (which could be inverted if we swap subtraction order) while the rest seem to have the same maximum scores (thus a lower difference)

##Q3.5

res_glycolysis_by_gluconeogenesis = t(sapply(glycolysis_data, getScoreSimple, ps_gluconeogenesis$pwm))


glycolysis_by_gluconeogenesis_maxscores <- data.frame(apply(res_glycolysis_by_gluconeogenesis, 1, max))
names(glycolysis_by_gluconeogenesis_maxscores)[1] <- "max_motif_score_gly_glu"
glycolysis_by_gluconeogenesis_maxscores


q3_5 <- cbind(glycolysis_maxscores, glycolysis_by_gluconeogenesis_maxscores)
q3_5$diff2 <- q3_5[,1] - q3_5[,2]
q3_5$genes <- rownames(q3_5)
colnames(q3_5)
q3_5_sort <- q3_5[order(q3_5$diff2, decreasing = TRUE),]
q3_5_sort#We get the same conclusion as Q3.4 here, and an explanation is probably that indeed some motifs might be more present in these 2 pathways since they are kind of similar (having also 14 common genes) and that the motifs dont vary greatly 


comparison <- merge(q3_4[,c(3,4)],q3_5[,c(3,4)], by="genes")
comparison <- comparison[order(comparison$diff2, decreasing = TRUE),]
comparison


#We can see that the differences when scanning with different pwm (gly with glu and glu with gly) we get the exact same difference (one beeing negative the other positive but we could just swap the subtraction order and it will be the same) besides 1 gene ENSG00000108515 (ENO3). This probably means that the best scores for those genes are the same which is to be expected. The difference with ENO3 might be some kind of variant?
