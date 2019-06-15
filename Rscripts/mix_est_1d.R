####  Fitting  1D  normal mixture models
####  Comparing  EM and direct Maximum Likelihood
####
#### --- this is a *subset* of MM's ./n1m-estim-2.R  -- for Lecture
##                                  ~~~~~~~~~~~~~~~

library(nor1mix)

plot(MW.nm9, lty=2, col = "blue", p.norm=FALSE, p.comp=TRUE)
MW.nm9 ## <- see the parameters (mu, sig2, w)

set.seed(17)
x9 <- rnorMix(5000, MW.nm9)
lines(density(x9), lwd=1.8)# "clearly" 3 components
## ok, so I hope any "normix" estimator should be able to get them ...

###----------------------------- MM's own approach ----------

## MM: I really think ML (without EM) and restrictions on sigma_k will work...
## --- goal: using optim() or mle() ==> Now implemented in
## ===>  nor1mix :: norMixMLE( ..... )
##                  ~~~~~~~~~           --> ~/R/Pkgs/nor1mix/R/norMixEM.R

### manual likelihood...


if(FALSE)##' now part of  'nor1mix' package (==> ~/R/Pkgs/nor1mix/R/norMixEM.R )
estep.nm <- function(x, obj, par)
{
    ## Purpose: 1 E-step for data 'x' and  obj = <norMix> object
    ##                               *or*  par = our parametrization
    ...........
}
## mstep.nm <- function(x, z)
## emstep.nm <- function(x, obj)


iniParNorMix0 <- function(x, m)
{
    ## Purpose: All equal components initial parameter for normal mixtures
    ## This is *TOO* simple (-> see
    ## --------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 17 Dec 2007, 18:47
    stopifnot(is.numeric(x), length(x) >= 1,
              is.numeric(m), length(m) == 1, m == as.integer(m), m >= 1)
    mu <- mean(x)
    s <- sd(x)
    nM2par(norMix(rep.int(mu,m), sigma = rep.int(s, m)))
}

# The trivial parametrization: pi_j = 1/K;  mu_j == all = global mean, etc
p0 <- iniParNorMix0(x9, 3)
par2norMix(p0) # looks "ok" (as expected)

plot(par2norMix(p0))
## log-Likelihoods :
llnorMix(p0, x9) # -8236.142
## at true parameter
llnorMix(nM2par(MW.nm9), x9) # -7778.055  (much better than "trivial initialization")

z0 <- estep.nm(x9, par=p0)
str(z0)
summary(c(z0)) ## all identical == 1/3
## .. i.e. starting with 'p0' -- EM always remains there...
## ..  ==> actually trivial if you think about it !!

## ==> something better :

##' Random Start at  m *different* observations:
iniParNorMix <- function(x, m)
{
    ## Purpose: equal *size* components, with means at random  x[i]
    ## This is *TOO* simple (-> see
    ## --------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 17 Dec 2007, 18:47
    stopifnot(is.numeric(x), length(x) >= 1,
              length(ux <- unique(x)) >= m, ## more unique obs. than components!
              is.numeric(m), length(m) == 1, m == as.integer(m), m >= 1)
    mu.s <- sample(ux, m) # take 'm' *different* x[i] as mu's
    s <- sd(x) / m
    nM2par(norMix(mu.s, sigma = rep.int(s, m))) # and pi. = weights = 1/m
}

set.seed(17)# as it is random --- NB: different sample() in R >= 3.6.0 !!!!
p0 <- iniParNorMix(x9, 3)
par2norMix(p0) # looks "ok", i.e., as expected

plot(par2norMix(p0), p.norm=TRUE, p.comp=TRUE)
## log-Likelihoods :
llnorMix(p0, x9) # -31833.75 ## much worse than trivial!
## at true parameter
llnorMix(nM2par(MW.nm9), x9) # -7778.055

## 1 E-step :
z0 <- estep.nm(x9, par=p0)
summary(c(z0)) ## now diverse in [0,1]

## Now an  M-step
p1 <- mstep.nm(x9, z=z0)
str(p1)
## or, the E- &  M- steps combined -- and working norMix objects
(P1 <- emstep.nm(x9, par2norMix(p0)))
plot(P1, p.comp=TRUE)
lines(hist(x9, plot=FALSE), freq=FALSE, col=adjustcolor("gray", 0.2))

## Another E- &  M- step :
(P2 <- emstep.nm(x9, P1))
plot(P2, p.comp=TRUE)
lines(P1, col="gray")

### Now, instead of guessing 'theta' (and then start with E-step),
## we guess the initial  Z  and then start with M- step rather than E-step
##    ---------------------

require(cluster)##--> Use cluster::clara() for initial 'Z' !

set.seed(1)
system.time(clp0 <- clara(x9, k=3, samples=1000, medoids.x = FALSE, rngR=TRUE))
## ~ 1 second
if(FALSE)## Note: could well use  pam() - it's relatively fast nowadays:
system.time(pam0 <- pam(x9, k=3, cluster.only = TRUE))
## 21.445 {nb-mm} // 26.973 {lynne, 2007} // 15.5 {lynne, 2009-11} // 8.4 {nb-mm4}
## 5.0 sec {lynne, ~2017}

clp0
## the 'clustering' component gives cluster/group labels  {1, 2, 3} :
table(clp0$clustering)  ## and split(..) builds a list of the 3 groups:
str( x9.Grps <- split(x9, clp0$clustering) )
nm1 <- clus2norMix(clp0$clustering, x9) ## <-- M-step from cluster:
nm1
plot(nm1, p.comp=TRUE) # and show the three clusters via color:
for(j in 1:3) rug(x9.Grps[[j]], col = j)

z2 <- estep.nm(x9, obj = nm1)
colSums(z2) # three groups of still about 1/3 each

## Try a few clara() calls and see the variation {!! ~= 1 sec/clara on nb-mm !!}
set.seed(11)
##for(i in 1:500) {
for(i in 1:100) {
    clp. <- clara(x9, k=3, samples=200, medoids.x = FALSE, rngR=TRUE)
    cat("."); if(i %% 50 == 0) cat("",i,"\n")
    (if(i == 1) plot else lines)(clus2norMix(clp.$clustering, x9))
}
## in general *quite* similar !

## now visualize a few E-M steps:
p.EMsteps <- function(nSteps, x, nmObj,
                      main = "EM for 1-dim normal mixtures",
                      sleep = nSteps <= 40, sec.sleep = 0.2, col2 = "tomato",
                      showLik = nSteps <= 80,
                      plotIni = FALSE, # plot the initial param. ?
                      doNote = TRUE)
{
    ## Purpose: show progress of EM - steps -- only few steps
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 26 Nov 2009
    stopifnot(is.norMix(nmObj), is.numeric(x), nSteps >= 1)
    nm. <- nmObj
    ll.em <- numeric(nSteps)
    for(k in 1:nSteps) {
        z. <- estep.nm(x9, obj = nm.)
        nm. <- do.call(norMix, c(mstep.nm(x=x9, z=z.), name=""))
        if(k == 1) { ## only after the first E+M - step
            op <- par(mar = c(0,0,0, 2.) + par("mar")); on.exit(par(op))
            plot(nm., p.norm=FALSE, main = main)
            if(plotIni) lines(nmObj, p.norm=FALSE)
            rug(x) ; mtext(paste("n = ", length(x)),
                           line = -1., adj=.99)
            nm2 <<- nm. ## save the 2nd
            if(showLik) {
                yVal <- approxfun(x=c(nSteps,1), y=par("usr")[3:4])
                xV <- { u <- par("usr")[1:2]; u[2] + .01*(u[2]-u[1]) }
            }
        } else lines(nm., p.norm=FALSE, col = rgb(0,0,0, 0.2))
        print(ll <- llnorMix(nM2par(nm.), x=x9))
        ll.em[k] <- ll
        if(showLik)
            text(xV, yVal(k), sprintf("%3d: %5.1f", k, ll),
                 cex = 0.8, adj = 0, xpd = NA)
        if(sleep) Sys.sleep(sec.sleep)
    }
    ## plot the last one more visibly
    lines(nm., p.norm=FALSE, col = col2, lwd = 2)
    if(doNote)
        mtext(paste(nSteps," E-, M-steps seem to have converged nicely .."),
              col=col2)
    ## return final normal mixture {but invisibly}:
    invisible(list(nm = nm., logLik = ll.em))
}

## start from clara() initialization  'nm1' :
p.EMsteps(20, x9, nm1, plotIni=TRUE)# 20 steps: final logLik ("ll"): -7779.5
## hmm, let's try more :
p.EMsteps(30, x9, nm1)


## aah, it seems to go on
r <- p.EMsteps(50, x9, nm1)
r <- p.EMsteps(200, x9, nm1)
## .. aah: needs *many* iterations -> finally gets to "good" solution
## 200 steps: final logLik ("ll"): -7771.226
nm. <- nm.200 <- r$nm # to keep

### hmm, converges to the same solution with a ``two wide'' middle peak
lines(density(x9, bw="SJ"), lwd=2, col="purple2")# "clearly" 3 components
##
lines(nm1, col = "brown")# the clara start (!)

## How did the log likelihood increase (it *does* !  famous EM theorem")
plot(r$logLik, type = "l", xlab = "iter", ylab = "log Likelihood",
     main = "EM - iterations --> Likelihood increases")
text(1:5, r$logLik[1:5])

## [ ............... ]





##---------------- Continue working with 'x9' example ---
tol <- 1e-6
oo <- norMix(mu = c(1,2,5))
llnorMix(nM2par(MW.nm9), x=x9)# -7778.055  (the "true" parameters)
llnorMix(nM2par(oo),  x=x9)# -14831.36 **FAR* away

if(FALSE) { ## for all these EM converges to pathological {2-comp.mix}:

oo <- norMix(mu = c(1,1,20))#-> EM "converges" to
##              mu         sig2            w
## [1,] 0.01692748 1.578665e+00 5.000000e-01
## [2,] 0.01692748 1.578665e+00 5.000000e-01
## [3,] 3.03457160 1.829534e-14 1.547347e-75

oo <- norMix(mu = c(-1,1,20))##-> conv to
## r; ll.em[i-1] # :
##             mu          sig2            w
## [1,] -1.062622  4.780516e-01 4.909407e-01
## [2,]  1.058053  4.322154e-01 5.090593e-01
## [3,]  3.034572 1.917842e-144 2.298869e-66
## Log.lik -7803.433

oo <- norMix(mu = c(-1,0,2))##  --- now converges fine
}

if(.Platform$OS.type == "unix")## use *FAST* device, not anti-aliased
x11(type="Xlib")

nSteps <- 2000 ## needs *many* iterations
ll.em <- rel.dp1 <- rel.dpI <- rep(NA_real_, nSteps)
r <- oo ; p <- nM2par(r)
for(i in 1:nSteps) { ## now with 'tol' effectively a while() loop
    p.o <- p
    plot(r <- emstep.nm(x9, r), p.comp=TRUE, p.norm=FALSE)
    ll.em[i] <- sum(dnorMix(x9, r, log=TRUE))
    ## our parametrization:
    p <- nM2par(r)
    dp <- abs(p - p.o)
    ## L1 and L_Inf relative errors -- w/ respect to *parameters*
    np <- sum(abs(p))
    rel.dp1[i]  <- sum(dp) / np
    rel.dpI[i]  <- max(dp) / (np/length(p))
    mtext(sprintf("it. = %4d; log.lik= %10.4f; rel.d(p)_1, _Inf = %10.4e,%10.4e",
                  i, ll.em[i], rel.dp1[i], rel.dpI[i]))
    if(max(rel.dp1[i], rel.dpI[i]) < tol) {
        ## Convergence w/ respect of relative *parameter* changes
        ##  (in both  L_1 and L_{\infty} norms) :
        cat("converged after",i,"iterations\n")
        break
    }
    if(i < 20) Sys.sleep(0.4)
}
##--> converged after 1283 iterations (!)

nm. <- nm.conv <- r # keep final
length(ll.em) <- length(rel.dp1) <- length(rel.dpI) <- i
ll.em <- c(sum(dnorMix(x9, oo, log=TRUE)), ll.em) # prepend LL(<start>)
d.ll <- abs(diff(ll.em)/((ll.em[-(i+1)]+ll.em[-1])/2))

if(names(dev.cur()) == "X11") ## back to use a "nice" (rather than fast) device:
    x11()# -> 'X11cairo'

plot(ll.em, type = "l", xlab = "iter", ylab = "log Likelihood",
     ylim = c(-8400, max(ll.em)),
     ## ylim = c(-7900, max(ll.em)),
     ## ylim = c(-7820, max(ll.em)),
     main = "EM - iterations")
text(1:9, ll.em[1:9])

par(new=TRUE)
plot(rel.dp1, type="l", ann=FALSE,axes=FALSE, log = "y", col="blue3")
lines(rel.dpI, col="blue4")
axis(4); mtext(expression(group("||", p[i] - p[i-1], "||")[1]), 4)
## The log.likelihood converges "much earlier":
lines(d.ll, col="green2")

dmat <- cbind(rel.dp1, rel.dpI, d.ll)
mcol <- c("blue2", "blue3", "green3")
matplot(dmat, type = "l", log = "y", ylim = c(1e-10, max(dmat)),
        xlab = "EM iterations", col = mcol, yaxt = "n"); eaxis(2)
legend("topright",
       expression(group("||", p[i] - p[i-1], "||")[1],
                  group("||", p[i] - p[i-1], "||")[infinity],
                  Delta(log(lik.(theta)))),
       col=mcol, lty=1:3)
1

###--- INSTEAD:  Maximum Likelihood  directly ----  using optim() or mle() --

neglogl <- function(par) -llnorMix(par, x=x9)
neglogl(nM2par(nm1)) # 7941 or  a bit better, ok

## [ ............... ]


or <- optim(nM2par(nm1), neglogl, method = "BFGS",
            control = list(maxit=6000, trace=TRUE, REPORT=1))
## quite (but not extremely) fast ..
## but  Rprof() gives that  84% of CPU is spent in dnorm(.)

## NOTA BENE:  Nowadays the *same* via
orM <- norMixMLE(x9, m=3, trace=4)

## When clara() had 1000 samples:
## initial  value 7882.153847
## iter  10 value 7780.962937
## iter  20 value 7771.274283
## final  value 7771.207469
## converged

## Start from the EM - solution:
or. <- optim(nM2par(nm.), neglogl, method = "BFGS",
            control = list(maxit=6000, trace=TRUE, REPORT=4))
## potentially is improved as well!! (depending where you stopped above)
## ---> 7771.207

.llnorMix(or$par,x9) ##  -7771.207  --- ok
(ob <- sort(par2norMix(or$par, name = "from optim(*, \"BFSG\")")))
plot(ob, p.comp=TRUE)## wow! much better [than the "wrong" mclust one]
lines(nm., col="blue")## not distinguishable
lines(density(x9, bw="SJ"), lwd=2, col="purple2")
## WOW!  This is MUCH better than the
## ---   __wrongly-converged__ EM-result (as from 'mclust' or 'flexmix')

### further examples
plot(nm3 <- norMix(c(1, -1), sig2 = c(.25,1), w= c(.2, .8)))

plot(nm4 <- norMix(c(0, 2),  sig2 = c(3,  1), w= c(.3, .7)), p.comp=TRUE)

x.2 <- rnorMix(1000, nm4)

plot(nm4, col="light blue", p.comp=TRUE, p.norm=FALSE)
lines(density(x.2))

## [ ............... ]
