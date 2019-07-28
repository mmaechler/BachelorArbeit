#### testing the various forms of calculating log likelihood of norMmix obj's

devtools::load_all("~/ethz/BA/norMmix")
library(MASS)

tr <- "clr1"
mo <- "EVV"
par. <- par212
x <- rnorMmix(MW212, n=5000)
p <- 2;k <- 2

# testing norMmix function for speed
system.time( for (i in 1:10000) {llnorMmix(par., x, p, k, trafo=tr, model=mo)} )
#    user  system elapsed 
#   2.140   0.004   2.144 

# compare to function using mvtnorm function
system.time( for (i in 1:10000) {llmvtnorm(par., x, p, k, trafo=tr, model=mo)} )
#    user  system elapsed 
#   5.228   0.004   5.231 

norMmixMLE( x, p=2, k=2, trafo="clr1", model="EVV" )

par. <- par213
x <- rnorMmix(MW213, n=10000)
norMmixMLE( x, p=2, k=2, trafo="clr1", model="VVV" )

par. <- par211
x <- rnorMmix(MW211)
norMmixMLE( x, p=2, k=2, trafo="clr1", model="VEE" )

p <- 2;k <- 1
par. <- par2nm1
x <- rnorMmix(MW2nm1)

norMmixMLE( x, p=2, k=1, trafo="clr1", model="VEE" )





#### testing MLE function

x <- rnorMmix(MW26, n=1000)

for (i in 1:5) {
	f <- norMmixMLE(x,p=2, k=i, trafo="clr1", model="EEI", maxiter=200, trace=1)
	res[i] <- f$optr$value
}

plot(-res, type="l")
plot(f$norMmix)
points(x)
g <- norMmixMLE(x,2,2,trafo="clr1",model="EEI",maxiter=100)
metric.norMmix(MW26,g$norMmix)
plot(MW26)
plot(g$norMmix, newWindow=FALSE)



x <- rnorMmix(MW213, n=12000)
h <- norMmixMLE(x,2,2, trafo="clr1", model="VVV", maxiter=300)
plot(MW213)
points(x)
plot(h$norMmix, newWindow=FALSE)
metric.norMmix(MW213,h$norMmix)

res <- 2:5

for (i in 10^(2:5)) {
	x <- rnorMmix(MW213, n=i)
	h <- norMmixMLE(x,2,2, trafo="clr1", model="VVV", maxiter=300,trace=1)
	res[i] <- (metric.norMmix(MW213,h$norMmix)$penalty)
}



