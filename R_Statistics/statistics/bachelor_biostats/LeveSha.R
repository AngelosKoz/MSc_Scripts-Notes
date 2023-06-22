D9/3

setwd("I:/Bstat/RR")
dataStrodium <- read.table("Anova1.txt", header=T, dec = ",")#anova deixnei an oi diafores einai shmantikes
attach(dataStrodium)
names(dataStrodium)
tapply(STRODIUM, AREA, mean) #sunarthsh san trito orisma thn efarmozei se kathe epipedo tou paragonta aerea, ara vriskei tis m.t edw gia kathe perioxh
tapply(STRODIUM,AREA, var) #eikona gia diaspores twn 5 deigmatwn apo tis perioxes
boxplot(STRODIUM~AREA)#na dw an exw paratupa shmeia kai an moiazoun katanomes/an exw isa h diamesos sth mesh orthogwnioy /// den vlepw megalh diafora stis summetries edw
S_A<- STRODIUM[AREA=="A"] #prepei na to kanw gia kathe deigma ksexwrista dioti theloyme to kathe deigma na einai aneksarthto kai theloyme KATHE deigma na exei proelthei apo K.K plithismo
S_B<- STRODIUM[AREA=="B"]#pairnei times opoy h perioxh einai A,B,C klp
S_C<- STRODIUM[AREA=="C"]
S_D<- STRODIUM[AREA=="D"]
S_E<- STRODIUM[AREA=="E"]
shapiro.test(S_A)#prepei na ftiaksw 5 diaforetika dianysmata me tis metrhseis tou strodium gia kathe perioxh gia na ta diavasei to test mas
shapiro.test(S_B)#vlepw an akolouthei K.K 
shapiro.test(S_C)#mhdenikh upothesh oti exoyn proelthei apo K.K ara ta dedomena mas prosarmozontai sthn KK 
shapiro.test(S_D)#megalo p-value, dexomai ara apokliseis apo kanonikathta den einai shmantikes
shapiro.test(S_E)
library(car)#levene einai gia diaspores!mporei na kanei anova an exoyme kentrarei ws pros mesh timh h diaspora einai gia mhdenikh upothesh diasporwn
leveneTest(STRODIUM~AREA)#deixnei an oi diafores einai shmantika diaforetikes
#median apo mono toy kai mpainei otan exoyme apokliseis apo kanonikothta, to F einai apo thn anova, kai to pvalue einai h shmantikothta toy statistikou, ara elegxos edeiksa oti den uparxei shmantika diaretikh diafora dioti
#exw timh statistikou ... kai to pvalue
ANOVA1 <- aov(STRODIUM ~ AREA) # to "," einai gia aneksarthta deigmata, otan h metablhth mas einai omadopoihmenh ws pros ta epipeda toy paragonta tote vazw "~"
summary(ANOVA1)#residuals mesa sthn omada, Fvalue=MSB/MSE (sfalma sto paranomasth), 3,95e-12--poly mikros arithmos
#mikro pvalue ara shmantikes diafores
#otan den exw omoiogeneia kanw bf h welch
library(onewaytests)
bf.test(STRODIUM ~ AREA, data=dataStrodium, alpha = 0.05,)
welch.test(STRODIUM ~ AREA, data=dataStrodium, alpha = 0.05,)
TukeyHSD(ANOVA1)#διαφορα στισ μτ, διαφορα στα Δ.Ε και το pvalue adjustive/prosarmozei h methodos to pvalue gia periorismo sfalmatos typoy 1 otan exoyme na kanoyme me pollaples sugkriseis.. 
#tis sugkrinw me epipedo shmantikothta, edw to 0,05 kai koitazw ana 2. an einai mikrotero toy 0.05 uparxei diafora, ektos apo C-B , D-B kai D-C
plot(TukeyHSD(ANOVA1))
pairwise.t.test(STRODIUM, AREA, p.adjust = "bonferroni")#to 1 den einai pote 1 einai stroggulemeno, omoiws me ta 0
library(DescTools)
DunnettTest(STRODIUM ~ AREA, data=dataStrodium)#to control mpainei to prwto







bird.pyr<-c(1.11,1.23,0.91,0.95,0.99,1.08,1.18,1.29,1.12,0.88)
bird.peu<-c(1.60,2.17,1.85,1.99,1.74,1.54,1.86,1.87,2.04,1.70)
bird.pso<-c(0.42,0.93,0.77,0.37,0.50,0.48,0.68,0.62,0.67,1.03)
#μοιράζονται το ίδιο περιβάλλον
mean(bird.pyr)
mean(bird.peu)
mean(bird.pso)
var(bird.pyr)
var(bird.peu)
var(bird.pso)
boxplot(bird.pyr)
boxplot(bird.peu)
boxplot(bird.pso)

setwd("I:/Bstat/RR")
birds<-read.table("Birds.txt", header=T, dec = ",")
attach(birds)
names(birds)
tapply(TIME,BIRD,mean)
tapply(TIME,BIRD,var)
boxplot(TIME~BIRD)
bird.pyr<-TIME[BIRD=="Pyr"]
bird.peu<-TIME[BIRD=="Peu"]
bird.pso<-TIME[BIRD=="Pso"]
shapiro.test(bird.pyr)
shapiro.test(bird.peu)
shapiro.test(bird.pso)
library(car)
leveneTest(TIME~BIRD)
anova<-aov(TIME~BIRD)
summary(anova)


setwd("I:/Bstat/RR")
L <- read.table("fyta.txt", header=T)  #δημιουργουμε ενα data frame με τα δεδομενα μας
L
names (L)
#Θεωρούμε οτι κάθε μέτρηση ειναι ανεξάρτητη και οτι η δειγματολιψεια έγινε τυχαία.
#Ετσι λοιπόν έχουμε σαν Η0 υπόθεση οτι η αναλογίας των τεσσάρων αυτών φαινοτύπων ειναι 9/3/3/1.
#Για να ελέγξουμε πόσο πιθανό ειναι να έχουν προέλθει παρατηρούμενες αναλογίες απο εναν πληθυσμο με αναλογία φαινοτύπων 9/3/3/1 θα χρειαστέι να εφαρμόσουμε ενα χ2 τεστ καλής προσαρμογής.
observed <- c(152, 39, 53, 6)#δημιουργουμε τωρα ενα διανυσμα με τις παρατηρησιμες τιμές που έχουμε στο δείγμα μας.
n <- sum(observed)# τεστάρουμε οτι το σύνολο τους ειναι 250, το οπόιο ειναι
n
p0 <- c(9/16, 3/16, 3/16, 1/16)# η αναλογια των αναμενώμενων τιμών.
expect <- n*p0
expect
chi2 <- sum((observed-expect)^2/expect)
chi2 #8.972444
chi.value <- qchisq(0.05, 3, lower.tail = F) # οι βαθμοί ελευθερίας της κατανομής ειναι 3 
chi.value #7.814728
#συνεπώς εχουμε περιοχή απόρριψης της τιμής του στατιστικού< του chi.value. αρα η υπόθεση απορριπτεται σε επιπεδο σημαντικοτητας 0.05


