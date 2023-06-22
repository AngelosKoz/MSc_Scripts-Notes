#analush diaspores me duo paragontes//prwth periptwsh paradeigma selida 7
#paragontas me dyo epipeda (control xwris ormonh kai h allh omada me ormonh)
#theloyme na doyme an oi paragontes epidroune sth sygkentrwsh asvestioy sta poylia
#5 metrhseis apo kathe sundiasmo (ars+thil control/ars+thil me ormonh)
#kathe parathrhsh grafetai san sundiasmos twn orwn a gia epipedo epidrashs a kai b gia b
#epeidh exoyme 1> metrhseis se kathe --- mporw na dw an exw allhlepidrash
#deuterh periptwsh  h epidrash tou enos paragonta den einai statherh se ola ta epipeda touy alloy paragonta , an einai statheri tote shmainei oti den exw allhlepidrash ara kanoyme prwta ton elegxo gia thn allhlepidrash kai
#-an broume oti den exw shmantikh allhlepidrash tote proxwrw na kanw ton elegxo gia th shmantikothta toy kathe paragonta , an h epistash tous einai statistika shmantikh
#otan den exw genika shmantikh allhlpidrash dhladh otan thelw na elegxw an upraxei diafora metaksu a/th tote pairnw oles tis metrhseis ars kai oles twn thil kai agnow ta alla 2 epipeda toy apargonta, kathws mas edeikse oti
#-den eksartatai apo thn paroysia touy alloy paragonta


data <- read.table("Anova2.txt", header=T)#metablhth apokrishs h allanhni kai theloyme na doyme an uparxei diafora metsku eidwn kai fullwn sth sygk allanhnhs
attach(data)#2 paragontes na elegksoyme an allhlepidroyn kai an oxi tote kathe ena ksexwrista
names(data)
tapply(alanine,list(sex,species),mean)
tapply(alanine,list(sex,species),var)
tapply(alanine,species,mean)
tapply(alanine,sex,mean)
barplot(tapply(alanine,list(sex,species), mean), beside=T, legend =T, ylim=c(0,30)) #tash ta arsenika na exoyn megaluterh sugkentrwsh
boxplot(alanine~ sex*species) #prwth eikona gia tis katanomes twn 6 deigmatwn akolouthoyn k.k me ish diaspora // to * einai oti tha parei ola ta epipeda mazi olous dhladh tous sundiasmous
f_s1<- alanine[sex=="f" & species=='s1'] #xwrizw pali ta dedomena gia na dw to kathena
f_s2<- alanine[sex=="f" & species=='s2'] 
f_s3<- alanine[sex=="f" & species=='s3'] 
m_s1<- alanine[sex=="m" & species=='s1'] 
m_s2<- alanine[sex=="m" & species=='s2'] 
m_s3<- alanine[sex=="m" & species=='s3'] 
par(mfrow=c(3,2))
qqnorm(f_s1, main="Normal Q-Q Plot for females of speceis 1")
qqline(f_s1)
qqnorm(m_s1, main="Normal Q-Q Plot for males of speceis 1")
qqline(m_s1)
qqnorm(f_s2, main="Normal Q-Q Plot for females of speceis 2")
qqline(f_s2)
qqnorm(m_s2, main="Normal Q-Q Plot for males of speceis 2")
qqline(m_s2)
qqnorm(f_s3, main="Normal Q-Q Plot for females of speceis 3")
qqline(f_s3)
qqnorm(m_s3, main="Normal Q-Q Plot for males of speceis 3")
qqline(m_s3)
par(mfrow=c(1,1))
tapply(data,shapiro.test)
shapiro.test(f_s1) #den kanw ola ta dedomena shapiro, (den exw proupothesh oti allanhni akoloythei kanonikh) alla oti kathe ena apo ta 6 exoyn proelthei apo k.k.
shapiro.test(f_s2)
shapiro.test(f_s3)
shapiro.test(m_s1)
shapiro.test(m_s2)
shapiro.test(m_s3)
library(car)
leveneTest(alanine ~ sex*species)
interaction.plot(sex, species, alanine) #idia aukshsh diathreitai den einai akribws parallhles alla einai sxedon parallhles ta euthigramma tmhmata ara oi endeikseis einai oti mallon DEN einai shmantikh h 
#-allhlepidrash --shmantikh allhkepidrash tha htan an eixa meiwsh se sxesh thhlikoy kai ars//h epidrash eidous px 1  den tha htan idia sta fula
interaction.plot(species, sex, alanine)
model1 <- aov(alanine ~ sex + species + sex:species)
summary(model1) #vlepoyme oti sto sex uparxei shmantikh diafora metaksu ars/thil kai omoiws sto eidos
model2 <- aov(alanine ~ sex*species)
summary(model2)
TukeyHSD(model1) #metaksu 2 kai 3 den uparxei shmantikh diafora

#ta block omadopoioyn me idia xarakthristika gia mia sunthiki (p.x sunthikes periballontos gia na megalwsoyn zwa opws einai h thermokrasia,ugrasia klp//exw se ena tetragwno aspoyme sugekrimeno xwma poy exei diaforetika futa na
#megalwnoyn se idies synthikes,diaforetikes metaksu twn block. sto ena px upshlh thermo sto allo block xamhlh//mas endiaferei h diaita kai oxi sta block)
setwd("I:/Bstat/RR")
data <- read.table("RBD_ANOVA.txt", header=T)#se kathe sxediasmo epipedou tou block exoyume mono 1 timh / h0 mesh aukshsh varous idia gia oles tis diaites, h1 oti h mesh aukshsh varous den einai idia
attach(data)#to montelo 2 orous a,b poy einai oi epidraseis paragontwn (pera apo th mesh timh)
names(data)
model <- aov(WEIGHT ~ DIET + BLOCK) #eksarthmenh metavlhth to weight, ws prws tous duo paragontes athristika
summary(model)#mesa athrismata prokuptoun ahtrismata/bathmoi.eleutherias,
#milame gia montelo 3 anova (1 paragontas einai fixed(autos poy elegxoyme-diaites)enw o paragontas block einai o tuxaios paragontas)
#otan h drash enos paragonta den einai idia se ola ta epipeda toy paragonta tote den kanoyme elegxo block, dioti exoyme mono 1 parathrhsh se kathe sundiasmo kelioy
#F=11.82, kai p.value polu mikro, ara h0 aporr, ara exoume shmantikh diafora anamesa stris tesseris diaites
interaction.plot(BLOCK, DIET, WEIGHT)#ta zig zag tha eprepe na einai parallhla an htan isa
#an apomonwsoume 1 kai 3 fainetai na mhn exoyme allhlepidrash, alla an ta paroyme ola mazi oi grammes temnontai metaksu tous pou shmainei oti mallon exoyme allhlepidrash
#einai elegxos pollaplwn sugrisewn kai den kanoume TukeySHD (allhlepidrash)
#sta block 1 2 kai 3 ta apotelesmata ths diaita menei sxedon idia enw exoyme aukshsh kai ptwsh sta block 4 5
#enw sth diaita 1 sto bloc k1 einai xamhla, auksanei kai meta peftei ara den idia se ola ta block 
#fwto sth d2 ,blepoyme oti exei kalutero "metabolismo, kai diatireitai h aukshsh-meiwsh stis diaites
#enw h D3 den paramenei statherh to epipedo ths se ola ta epipeda tou paragonta, tote leme oti oi paragontes allhlepidroyn metaksu tous OXI TA EPIPEDA (oi paragontes block kai epipeda)
#o skopos pou ginetai autop einai dioti theloyme na meiwsouyme thn anekshghth metablhtothta (residuals)
models <- aov(WEIGHT ~ DIET)
models
summary(models)


