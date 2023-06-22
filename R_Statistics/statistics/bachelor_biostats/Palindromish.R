#Palindromish
#μονο στη πανω για γραμμικη
#πρεπει να κανω διαγραμμα διαχυσης πριν κανω παλινδρομιση
#αν δεν ξεφευγουν απο ζωνη [-2,2] δηλαδη το 95%
#ολα ισα διασπορα

#D10/13 

Temp <- c(0, 4, 10, 15, 21, 29, 36, 51, 65)
Y <- c(67.4, 72.0, 74.3, 80.6, 83.7, 92.9, 95.4, 107.6, 120.6)
plot(Temp, Y, xlab="Θερμοκρασία (oC)", ylab = "Διαλυτότητα")
lm(Y ~ Temp)
fit <- lm(Y ~ Temp)
plot(Temp,Y)
abline(fit)
summary(fit) #t.value timh statistikoy, kai dipla to pvalue
confint(fit, level = 0.95)
anova(fit) #(temp = regression, residuals = sfalmata, )anova + ttest isodunamoi otan milame gia apllh gram palin dhladh 1 eks metablhth
predict(fit,list(Temp =18))
predict(fit,list(Temp =c(5, 18, 25)))
predict(fit, list(Temp=18), interval="confidence", level=0.95)
predict(fit, list(Temp=18), interval="prediction", level=0.95)
plot(fitted(fit), rstandard(fit), ylim=c(-2,2), xlab="predicted values", ylab="standardized residuals")
abline(h=0)
#οταν θελω να σχεδιασω στην ιδια γραφ παρασταση πολλες γραμμες, χρησιμοποιω την συναρτηση abline
#residual standardized error mas deixnei thn r2
#theloyme paromoies apostaseis, h na auksanoyn sfalmata 
#apo residuals shapiro??? gia na dw isws kanonikothta
#kuriws theloyme omoiogeneia diasporwn


#D10/14

Area <- c(208, 4.33, 3.83, 3.34, 0.62, 0.53, 0.92, 0.39, 0.018, 0.016, 0.003, 0.002)
Species <- c(42, 21, 20, 16, 15, 14, 14,10, 4, 7, 7, 6)
plot(Area, Species)
LogArea <- log10(Area)
LogSpecies <- log10(Species)
plot(LogArea, LogSpecies)
lm(LogSpecies ~ LogArea)#το πιο σημαντικο, και μετα να το βαλω σε fit
fit <- lm(LogSpecies ~ LogArea)
plot(LogArea,LogSpecies)
abline(fit)
summary(fit) #ελεγχος σηνμαντικοτητας για β χρησιμοποιω το 8,369 και pval ,αρα β=0 απορριπτεται//το ιτνερσεπτ ειναι το α/σταθερα μας///
#το λογκ area einai η μεταβλητη//ελεγχος γινεται με βαση τη κλιση
anova(fit)
confint(fit)#D.E
plot(fitted(fit), rstandard(fit), ylim=c(-2.5,2.5), xlab="predicted values", ylab="standardized residuals")
abline(h=0)

2 graf parastaseiw

x<-c(0.2,1.4,3.1)
y<-c(0.08,2,4)
x1<-seq(0,4,0.1)# apo 0 mexri 4 me bhmata 0.1
y1<-x1 #gia ayta ta x ypologizw ta y
y2<-x1^2 #%
plot(xr,yr)
lines(x1,y1) #vazw y=x
lines(x1,y2) #vazw y=x^2


Temp<-c(-18,-15,-10,-5,0,5,7,10,15,19)
Oxy<-c(5.2,4.7,4.5,3.9,3.4,3.1,3.0,2.7,2.5,1.8)
plot(Temp, Oxy, xlab="Temperature (oC)", ylab = "Oxygen(mg/g/h)")
fit <- lm(Oxy ~ Temp)
plot(Temp,Oxy)
abline(fit)
summary(fit)
anova(fit)
predict(fit, list(Temp=12), interval="confidence", level=0.95)
predict(fit, list(Temp=12), interval="prediction", level=0.95)
xx1<-predict(fit, list(Temp=c(-18,-15,-10,-5,0,5,7,10,15,19)), interval="confidence", level=0.95)
xx2<-predict(fit, list(Temp=c(-18,-15,-10,-5,0,5,7,10,15,19)), interval="prediction", level=0.95)
x1<-c(4.897012, 5.256196, 4.662120, 4.981533,4.267348, 4.527048,3.865490, 4.079649,3.451395, 3.644486,3.021204, 3.225420,2.844934, 3.061987,2.577106, 2.820261,2.124434, 2.423674,1.758823, 2.109880)
x2<-c(4.722983, 5.430225,4.477878, 5.165775,4.066055, 4.728341,3.649675, 4.295463,3.228385, 3.867496,2.802032, 3.444591,2.630084, 3.276837,2.370696, 3.026670, 1.934672, 2.613437,1.582777, 2.285926)
plot(Temp,Oxy,xlab="Temperature (oC)", ylab="Oxygen(mg/g/h)")
abline(fit)
arrows(Temp, xx1[,2], Temp, xx1[,3], length=0.05, angle=90, code=3, col="firebrick 2")
arrows(Temp, xx2[,2], Temp, xx2[,3], length=0.05, angle=90, code=3, col="steelblue 4")
points(y=xx1[,1],x=Temp, col="darkgreen")
points(y=xx2[,2],x=Temp, col="darkorange2")
abline(fit)


plot(Temp,Oxy, pch=19, main="Oxygen vs Temperature measurements", xlab="Temperature (Celsius)", ylab="Oxygen", ylim=c(min(xx2[,2]),max(xx2[,3])))
arrows(Temp, xx1[,2], Temp, xx1[,3], length=0.05, angle=90, code=3, col="firebrick 2", lwd=2)
lines(Temp,xx1[,2],lty=2, col="firebrick 2",lwd=2)
lines(Temp,xx1[,3],lty=2, col="firebrick 2",lwd=2)
arrows(Temp, xx2[,2], Temp, xx2[,3], length=0.05, angle=90, code=3, col="steelblue 4", lwd=2)
lines(Temp,xx2[,2],lty=2, col="steelblue 4",lwd=2)
lines(Temp,xx2[,3],lty=2, col="steelblue 4",lwd=2)

plot(Temp,Oxy,pch=19, xlab="Temperature (Celsius)", ylab="Oxygen")
abline(fit,lwd=4, col="darkgreen")
lines(Temp,xx1[,2], col="firebrick 1",lwd=2)
lines(Temp,xx1[,3], col="firebrick 1",lwd=2)
lines(Temp,xx2[,2], col="cyan",lwd=2)
lines(Temp,xx2[,3], col="cyan",lwd=2)
arrows(Temp, xx1[,2], Temp, xx1[,3], length=0.05, angle=90, code=3, col="darkred", lwd=2)
arrows(Temp, xx2[,2], Temp, xx2[,3], length=0.05, angle=90, code=3, col="navy", lwd=2)

install.packages("car")

