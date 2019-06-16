### testfile for learning various clustering packages
##nor1mix package

library(nor1mix)


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


x1 <- rnorMix(500, MW.nm1)
hist(x1)

x10 <- rnorMix(5000, MW.nm10)
hist(x10) #not visible
hist(x10, breaks = 100)#visible
lines(density(x10), lwd= 2)#not strongly visible

x10 <- rnorMix(50000, MW.nm10)
lines(density(x10), lwd= 2)#slow improvement for exp increase in data


