# Question 1
# eixa 8ema me to grapsimo diavasma arxeiwn. Gia na mh xasw monades
# egrapsa me to xeri tis times kai synexisa thn askhsh
cells <- c(7.3, 7.7, 8.2, 7.6, 7.9, 7.6, 7.2, 7.2, 7.9, 7, 6.9, 7.7)

m <- mean(cells) # mean value

var(cells) # variance

s <- sd(cells) # standard deviation
s

100*sd(cells)/mean(cells) # syntelesths metavlhtothtas 

N <- 12 # mhkos dedomenwn
conf.interval <- qnorm(0.95)*s/sqrt(N) # confidence interval
low_lim = m - conf.interval # lower limit of conf.interval
high_lim = m + conf.interval # upper limit of conf.interval
# print values
conf.interval
low_lim 
high_lim


# Question 3
pnorm(25, 30, 15/sqrt(25), lower.tail = T) # answer = 0.04779035

