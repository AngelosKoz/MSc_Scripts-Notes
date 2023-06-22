D7/7

setwd("I:/Bstat/RR")
Temp<-read.table("crabs.txt",header=T,dec=",")
t.test(Temp$T, mu=24.3, alternative ="two.sided")
m<-mean(Temp$T)
s<-sd(Temp$T)
n<-length(Temp$T)
t<-(m-24.3)/(s/sqrt(n))
p1<-pt(abs(t), df=n-1, lower.tail = F)
cat("p-value =", 2*p1)
t.value1<-qt(0.025,24,lower.tail=F)
t.value1
t.value2<-qt(0.005,24,lower.tail=F)
t.value2
power.t.test(n=25,delta=1,sd=s,sig.level=0.05,power=NULL,type="one.sample", alternative = "two.sided")
power.t.test(n=NULL,delta=1,sd=s,sig.level=0.05,power=0.9,type="one.sample",alternative="two.sided")
power.t.test(n=NULL,delta=1,sd=s,sig.level=0.01,power=0.9,type="one.sample",alternative="two.sided")
power.t.test(n=NULL,delta=1,sd=s,sig.level=0.05,power=0.8,type="one.sample",alternative="two.sided")
power.t.test(n=NULL,delta=1,sd=s,sig.level=0.05,power=0.9,type="one.sample",alternative="two.sided")
power.t.test(n=25,delta=1,sd=s,sig.level=0.01,power=NULL,type="one.sample", alternative = "two.sided")
s2<-var(Temp$T)
chi2<-(n-1)*s2/2
chi2
p2<-pchisq(chi2,df=n-1,lower.tail=T)
cat("p-value =", p2)
x.value<-qchisq(0.95, 24, lower.tail = F)
x.value

D8/3

setwd("I:/Bstat/RR")
data1 <- data.frame(sex=c("male","male","male","male","male","male","male","female","female","female","female","female","female"),
cholest =c(221,219,230,220,222,224,227, 223,221,230,224,223,231))
data1
names(data1)
attach(data1)
tapply(cholest, sex, length)
tapply(cholest, sex, mean)
tapply(cholest, sex, var)
tapply(cholest, sex, sd)
tapply(cholest, sex, summary)
cholest_male <- cholest[sex=="male"]
cholest_female <- cholest[sex=="female"]
par(mfrow=c(1,2))
qqnorm(cholest_male, main="Normal Q-Q Plot for males")
qqline(cholest_male)
qqnorm(cholest_female, main="Normal Q-Q Plot for females")
qqline(cholest_female)
par(mfrow=c(1,1))
shapiro.test(cholest_male)
shapiro.test(cholest_female)
var.test(cholest_male , cholest_female )
library(car)
leveneTest(cholest, sex, center="mean")

#D8/12

n <- 60; x <- 36; p0 <- 0.5
z <- (x/n-p0)/sqrt(p0*(1-p0)/n)
cat("?? ???????? ?????? ?????????????????????? z=", z ,"\n")
z.value <- qnorm(0.025, 0, 1, lower.tail = F)
p.value <- 2* pnorm(abs(z), lower.tail=FALSE)
n1 <- 2030; x1 <- 84; p1 <- x1/n1
n2 <- 2051; x2 <- 56; p2 <- x2/n2
phat <- (x1+x2)/(n1+n2)
s <- sqrt(phat*(1-phat)*(1/n1+1/n2))
z.stat <- (p1-p2)/s
cat("?? ???????? ?????? ?????????????????????? z=", z.stat ,"\n")
z.value <- qnorm(0.05, 0, 1, lower.tail = F)
p.value <- pnorm(z.stat, lower.tail=FALSE)



#D8/19

#PARADEIGMA 3// H0 oti pithanothta tha einai idia na epileksoun kathe trofh dhladh 
#H0: P1=1/5,P2=1/5,P3=1/5,P4=1/5,P5=1/5
#H1: Pi$1/5 gia kapoio i----xrhsimopoiw x² kalhs prosarmoghs
#????=35*1/5=7 (???? ?????????? ???????? ???????? ?????????????????? ???? ???????????? ???? ???????? 7 ???? ???????? ??????????)



obs<-c(8,13,6,6,2)#parathroumenes diafores)
n<-sum(obs)
p0<-c(1/5,1/5,1/5,1/5,1/5)#pithanothtes gia to kathena)
expect<-n*p0 #anamenomenes)
chi2<-sum(((obs-expect)^2)/expect)
#akolouthei x2 katanomh me k-1 b.e dhladh 5-1=4 b.e
#thelw posotiaia shmeia ara quantile
chi.value<-qchisq(0.05,4,lower.tail=FALSE)#einai ligo pio mikro ara den anoikei sthn perioxh aporipshs poy einai oi times panw apo thn timh auth
#ara den aporriptetai se epipedo shmantikothtas 95%
chi2>chi.value
#mas leei an einai megalutero h oxi 
chisq.test opoy dinw x=dianysma me parathroymenew, p=p0 dhladh dianusma me pithanothtes, obs=observed

