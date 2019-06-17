### testfile for learning various clustering packages
##nor1mix package


library(nor1mix)


#several mixtures

plot(MW.nm2) #skew
plot(MW.nm3) #strong skew
plot(MW.nm4) #kurtotic
plot(MW.nm5) #outlier
plot(MW.nm6) #bimodal
plot(MW.nm7) #separated
plot(MW.nm8) #asym bim
plot(MW.nm9) #trimodal
plot(MW.nm10)#claw
plot(MW.nm11)#double claw
plot(MW.nm12)#asym claw
plot(MW.nm13)#as do claw
plot(MW.nm14)#smooth comb
plot(MW.nm15)#disc comb
plot(MW.nm16)#dist bim


#rnorMix function

x1 <- rnorMix(500, MW.nm1)
hist(x1)

x10 <- rnorMix(5000, MW.nm10)
hist(x10) #not visible
hist(x10, breaks = 100)#visible
lines(density(x10), lwd= 2)#not strongly visible

x10 <- rnorMix(50000, MW.nm10)
lines(density(x10), lwd= 2)#slow improvement for exp increase in data


# clus2normix?

# dnorMix

dnorMixL(MW.nm10) #returns about 500 density values(why 500?)
dnorMixL(MW.nm10, n = 600) #ah

x10 <- rnorMix(5000, MW.nm10)
dnorMix(MW.nm10, seq(-3,3, length=200)) #deprecated use of fctn use dnorMixL

absc <- seq(-3,3,length=511)
ffs <- dpnorMix(absc, MW.nm10) #list with "d" & "p" values
plot(absc, ffs$p) # for example


# llnorMix



