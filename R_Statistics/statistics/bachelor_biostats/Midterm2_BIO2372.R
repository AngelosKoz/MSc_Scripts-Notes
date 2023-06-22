#Askhsh 1 Aggelos Kozonakis bio2372
setwd("I:/Bstat/RR")
plant<-read.table("question1.txt", header=T, dec=",") #tha ginei elegxos gia diafora plhthusmiakwn meswn timwn. opou H0:μfse=μfsp, oti dhladh den exoyme diafora sto meso epipedo molubdoy kai H1:μfse diaforetiko apo μfsp
attach(plant) #gia na ginei o elegxos, prepei prwta na doyme ean kathe ena apo ta futa mas proerxetai apo kanonikh katanomh kai sth synexeia na doyme th sxesh sth diaspora toys
tapply(Pb,SPECIES,mean) 
tapply(Pb,SPECIES,var)
fse<-Pb[SPECIES=="Fse"] #vazoyme ta duo futa se ksexwristes metablhtes wste na ginei ksexwrista gia to kathe ena o elegxos kanonikothtas
fsp<-Pb[SPECIES=="Fsp"]
par(mfrow=c(1,2)) #ayto to vhma den xreiazetai parolauta mas dinei mia idea pws einai katanemhmena se sxesh me th grammh
qqnorm(fse, main="Normal Q-Q Plot for Fse")
qqline(fse)
qqnorm(fsp, main="Normal Q-Q Plot for Fsp")
qqline(fsp)
par(mfrow=c(1,1))
shapiro.test(fse) #me to shapiro tha doyme twra thn kanonikothta. Pairnoyme p-value=0.8766 to opoio einai megalutero apo to a=0.05 poy einai to epipedo shmantikothtas mas
shapiro.test(fsp) #omoiws p.value=0.9932 poy einai pali megalytero apo to a=0.05 poy shmainei oti kai ta duo proerxontai apo kanonika katanemhmeno plhthusmo

library(car)
leveneTest(Pb,SPECIES,center="mean") #twra tha doume th sxesh twn diasporwn wste na prosarmosoume analoga to t.test mas kai blepoyme oti F=1.3906, p.value=0.2612 pou einai megalytero apo to a=0.05 (epipedo shmantikothtas) ara
#mporoyme na poyme pws pithanws den uparxei statistika shmantikh diafora opote mporoyme na poyme pws einai ises oi diaspores
t.test(fse,fsp, paired=F,var.equal=T,conf.level=0.95) #ginetai twra to t.test me aneksarthta deigmata, ises diaspores, kai se epipedo shmantikothtas a=0.05
#pairnoyme t=2.3181,df=12,p.value=0.0389 pou shmainei oti h mhdenikh mas upothesh H0 aporriptetai se epipedo shmantikothtas a=0.05, dhladh oti yparxei statistika shmantikh diafora sto meso epipedo molubdoy metaksu fse kai fsp

#Askhsh 2
amino<-read.table("question2.txt", header=T,dec=",") #theloyme na doyme an se epipedo shmantikothtas a=0.05, uparxei statistika shmantikh diafora sth mesh sugkentrwsh MDA gia ta tria aminoksea. O elegxos autos tha ginei me 
#anova kathws exoyme 3 metablhtes, kai efoson plhrountai oi proupotheseis (exoyn epilagei tuxaia, prepei kathena na akoloythei kanonikh katanomh kai na exoyn ises diaspores) pame na kanoyme thn anova
attach(amino) #H0 einai oti h mesh sugkentrwsh MDA einai idia sta pontikia gia kathe aminoksu (mCarnosine=mHistidine=mImidazole to m einai ellhniko m apla gia na mhn uparksei thema to afhsa etsi), enw H1 oti h mesh sygkentrwsh
#MDA diaferei.
anova<-aov(MDA~AMINO)
summary(anova) #pairnoyme F2,12=20.65 kai p.value=0.00013, poy einai poly mikrotero apo to epipedo shmantikothtas a=0.05 sunepws mporoyme na poyme oti h H0 mas aporriptetai dhladh den einai idia h sugkentrwsh MDA ara twra 
#prepei na doyme me to test pollaplwn sygkrisewn poio diaferei me poio (tukey)
TukeyHSD(anova) #apo tis times vlepoyme pws h mhdenikes upotheseis pou aporriptontai einai metaksu Histidine-Carnosine (p.val=0.0002147) kai Imidazole-Histidine (p.val=0.0005696)
#enw den vlepoyme statistika shmantikes diafores metaksu Imidazole-Carnosine (p.val=0.8144940)

#Askhsh 3
wat<-read.table("question3.txt", header=T,dec=",")
names(wat)#kathws den mporw na dw to plot tha sas grapsw opws mporw th diadikasia.Gia na ginei analush grammikhs palindromhshs prepei prwta na kanoyme to diagramma diaxushs. Me to plot tha doyme tis times tis driasthriothtas
#toy neroy gia diaforetikes sugkentrwseis sakxarozhs, kai me vash ta shmeia ayta, tha pame sth sunexeia me thn sunarthsh lm
water<-c(0,5,10,15,20,25,30,35,40,45,50,55,60)
y<-c(0.998,0.997,0.991,0.988,0.985,0.982,0.978,0.97,0.97,0.963,0.957,0.953,0.951)
plot(water,y, xlab="Συγκέντρωση Σακχαρόζης (w/v)", ylab = "Δραστηριότητα Νερού") #Error in plot(Water, Y, xlab = "S<U+03C5><U+03B3><U+03BA><U+03AD><U+03BD>t<U+03C1><U+03C9>s<U+03B7> Sa<U+03BA><U+03C7>a<U+03C1><U+03CC><U+03B6><U+03B7><U+03C2> (w/v)",  : object 'Water' not found
lm(y ~ water) #tha sas steilw mia fwtografia na deite kai eseis an ekana kati lathos, tha sas steilw kai to arxeio poy onomasa htan auto me tis paules
fit <- lm(y ~ water) # tha baloyme thn eutheia elaxistwn tetragwnwn ths driasthriothtas nerou epi ths sugkentrwshs sakxarozhs
plot(water,y)#blepoyme oti den uparxoyn shmantika "paratupa" shmeia einai ola konta sthn eutheia elxaistwn tetragwnwn
abline(fit) #H0, a=0 kai enallaktikh H1 oti a#0.(t=1039.79,p=2e=16 ara aporriptetai) kai omoiws gia th shmantikothta ths eksiswshs elaxistwn tetragwnwn H0:b=0 kai H1:b#0 (t=-30.48,p=5.6e-12 kai omoiws aporriptetai) 
#kai stis duo periptwseis tha sugkrinoyme to p.value me to epipedo shmantikothtas a=0.05 alla kai stis duo einai mikrotero ara aporriptontai kai oi dyo
summary(fit) #apo edw tha pername tis times statistikou kai ta p.value
confint(fit, level = 0.95) #apo edw fainetai to diasthma empistosunhs gia th stathera kai thn klhsh eutheias
