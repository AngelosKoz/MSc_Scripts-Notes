#thetikh susxetish otan auksanei to ena auksanei to allo
#statistikes methodoi gia na doume an uparxei sxesh metaksu 2 metavlhtwn.
#Susxetish: prosdiorizei an uparxei GRAMMIKH sxesh metaksu 2 metablhtwn alla oxi ton tropo eksarthshs
#vriskoyme to suntelesth gia na doyme an h mia einai apotelesma ths allhs aitio-apotelesma
#Palindromhsh: x elegxomenh y eksarthmenh-agnwsth metablhth
#prwta kanw diagramma diaxyshs oste na mas apikonisei an h sxesh metaksu twn 2 metablhtwn einai grammikh
#Suntelesths pierson//parametrikos suntelesths posotikopoiei grammikh sxesh metaksu twn 2 metablhtwn den apodiknyei sxesh aitioy-apotelesmatos
#H timh pou dinei to r einai [-1,1] kai ektimaei plhthusmiako sunt.susx: to ??//??=?? den exw grammikh sxesh
#Elegxos shmantikothtasHo: ??=0,H1 ??/-0// Gia r=0.8 exw sysxetish kai me p-value vlepw an einai statistika shmantiko
#X kai Y prepei na einai Kanonika Katanemhmenes
#times sta akra = thetikh grammikh susxetish = auksanei h mia oso meiwnei? h allh
#Akraies times H den plhreitai Kanonikothta, xrhsimopoiw mh parametrikous suntelestes Spearman/Kendel tau rank
#Anakaluptoyn monotones sxeseis-susxetiseis den xrhsimopoioun ta x kai y ara xanoyn se isxy
#
Paradeigma 1
#gia pearson prepei na kanw kai elegxous

#D10/8

setwd("I:/Bstat/RR")
data<- read.table("Correlation.txt", header=T, dec=",")
attach(data)
data
names(data)
plot(Tail_length, Wing_length, xlab="Tail length (cm)", ylab = "Wing length (cm)", xlim=c(7,8.5), ylim=c(10,11.5))
#fainetai na einai grammikh kai psaxnw twra to vathmo, dhladh upologizw pearson + meta 
#boxplot gia outliers
par(mfrow=c(1,2))
boxplot(Wing_length, ylab = "Wing length (cm)", ylim=c(10,11.5))
boxplot(Tail_length, ylab="Tail length (cm)", ylim=c(7,8.5))
par(mfrow=c(1,1))
shapiro.test(Wing_length) #deixnei oti einai grammikes ara den aporiptw
shapiro.test(Tail_length)#omoiws
#mhdenikh upothesh oti ??=0//b.e to n-2//p-value h pithanothta to statistiko na parei timh sta 2 akra ths katanomhs apo auth poy edwse to deigma maw, edw to 5,58
cor(Wing_length, Tail_length) #pearson, me upsulh timh ara eimaste ok Efoson to statistiko poy xrhsimopoioyme gia elegxo akolouthei student katanomh kai to pvalue proerxetai apo upologismo Statistiko + pithanothta statistikoy
#na parei times megaluteres apo 5.589
#krisimo shmeio:
qt(0.025,10,lower.tail = F)
qt(0.025,10,lower.tail = T)
#Vrisketai to statistiko 5.58 ektos toy diasthmatos ara aporriptetai
cor(Wing_length, Tail_length, method="spearman") # Spearman
cor(Wing_length, Tail_length, method="kendall") # Kendall ties=epanalamvanomenes times oti den einai akribhs o upologismos alla me epanaliptikh diergasia//ta deigmata exoyn 2ples h 3ples times mesa sta eigmata
cor.test(Wing_length, Tail_length, method = "s")
cor.test(Wing_length, Tail_length, method = "k")
#gia ??=0.5?
#Dipleuros elegxos :
z<-sqrt(12-3)*((1/2*log10(1.87/0.13))-1/2*log10(1.5/0.5))
z
z.val<-qnorm(0.025,0,1,lower.tail=F) #0.025 dioti einai dipleuros kai oxi 0.05 (dipleuros einai to ??#??0)
z.val
#Vriskw z=2,35 kai zval vlepw oti diaferei shmantika kai einai logiko na einai diaforetiko apo to 0.5

#y=a+bx+e to e einai to sfalma dhladh diafores metaksy parathroymenos apo problepomenes
#E(Y)=a+bx
#E(e)=0
#yi kapelo einai ektimomenh, kai yi einai h timh poy pairnw kai etsi prokyptei y(kapelo)=a(kapelo)+b(kapelo_x <-eutheia elaxistoy tetragwnoy kai einai ektimhsh toy monteloy grammikhw palindromhshs)

#D10/13

#Kane diagramma diaxushs//Statistika//ektimhsh parametros kai statistiko elegxo shmantikothtas me H0 b=0 kai H1 b#0 . Ean den aporipsw mhdenikh upothesh tote exw E(Y)=a+bx kai tote den einai shmantikh 
Temp <- c(0, 4, 10, 15, 21, 29, 36, 51, 65)
Y <- c(67.4, 72.0, 74.3, 80.6, 83.7, 92.9, 95.4, 107.6, 120.6)
plot(Temp, Y, xlab="?????????????????????? (oC)", ylab = "??????????????????????")#thetikh klhsh dioti blepoyme anebainei to ena oso to allo
lm(Y ~ Temp)#mas dinei suntelestes/ektimhseis// to 0.8021 einai h klish Kai mexri stigmhs exw kanei diagramma diaxushs
fit <- lm(Y ~ Temp)
plot(Temp,Y)
abline(fit) #eutheia elaxistwn tetragwnwn//to abline vazei thn eutheia auth
summary(fit) #interecept einai parametros a//S.E//t.value koita LAB 10 H0, t7(to 7 einai bathmoi eleutherias)
#me to a poy prokuptei h eutheia den pernaei apo arxes aksonwn//to b deixnei th shmantikothta ths eksiswshs kai h timh deixnei pws diaferei apo to 0 satistika shmantika
confint(fit, level = 0.95)
#problepseis xrhsimopoiwntas thn eutheia 
predict(fit,list(Temp =18))#times ths aneksarthths metablhths = 18 gia na doyme thn problepsh
predict(fit,list(Temp =c(5, 18, 25)))#apla exw dianusma problepsewn
#an isxue h mhdenikh ypothesh h eutheia tha pernouse apo akswna x
#psaxnw eksarthsh ths y apo x, mas endiaferei h klhsh to b, an den diaferei shmantika apo to 0 tote h eutheia den exei grammikh sysxetish

#Eisagwgh grafhma diasporas sthn excel kai blepw th sxesh tous kai meta pataw prosthiki grammh tashs panw sta shmeia poy exw//probolh ...gia r kai mas emfanizei to a