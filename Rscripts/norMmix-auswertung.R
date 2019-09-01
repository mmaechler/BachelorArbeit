## 

devtools::load_all("~/ethz/BA/norMmix")
options(error=recover)
options(error=NULL)
set.seed(2019)
library(MASS)


models <- c("EII","VII","EEI","VEI","EVI",
	    "VVI","EEE","VEE","EVV","VVV")


n <- 2000
x <- rnorMmix(n=n,MW213)

res <- {}

for (i in models) {
	res[i] <- norMmixMLE(x,2,2,trafo="clr1",model=i)$optr$value
}

plot(-res,type="l")

res <- {}

for (k in 1:5) {
	res[k] <- norMmixMLE(x,2,k,trafo="clr1",model="EVV")$optr$value
}

bic <- matrix(0,5,length(models))
colnames(bic) <- models

aic <- matrix(0,5,length(models))
colnames(aic) <- models

p <- ncol(x)

for (i in models) {
	for (k in 1:5) {
		ans <- norMmixMLE(x,p,k,trafo="clr1",model=i,maxiter=200)
		bic[k,i] <- ans$parlen*log(n) +2*ans$optr$value
		aic[k,i] <- ans$parlen*2 +2*ans$optr$value
	}
}


overmod <- function(x,k) {
	bic <- {}
	aic <- {}
	n <- nrow(x)
	p <- ncol(x)
	for (i in models) {
		ans <- norMmixMLE(x,p,k,trafo="clr1",model=i)
		bic[i] <- ans$parlen*log(n)-2*ans$optr$value
		aic[i] <- ans$parlen*2-2*ans$optr$value
	}
	list(bic=bic, aic=aic)
}



####
#-------------------------------------------------------------------------------
####
# work on 2019-08-03
####

x <- rnorMmix(n=1000,MW213)

ans <- norMmixMLE(x,2,2,trafo="clr1",model="EVV")

det(ans$norMmix$Sigma[,,1])
# [1] 2.137754
det(ans$norMmix$Sigma[,,2])
# [1] 3.10791
## not the same...

## fixed some stuff with j-sum(j)/p instead of j-sum(j) 

x <- rnorMmix(n=1000,MW213)

ans <- norMmixMLE(x,2,2,trafo="clr1",model="EVV")

det(ans$norMmix$Sigma[,,1])
# [1] 2.552217
det(ans$norMmix$Sigma[,,2])
# [1] 2.552217
## finally produces equal alpha!

ans <- fit.norMmix(x,k=1:10,models=1:10, trafo="clr1")

ans$BIC
#          EII       VII       EEI       VEI       EVI       VVI
# 1  16593.547 16593.547 16606.824 16606.824 16606.824 16606.824
# 2  10992.291 10968.495 10306.654 10203.404 10274.819 10167.386
# 3  10688.971 10632.723 10296.067  9899.632 10284.475  9858.474
# 4  10292.573 10278.963  9638.813  9576.585  9637.941  9573.378
# 5  10177.416 10242.754  9284.409  9274.085  9297.188  9376.938
# 6  10334.018  9660.044  9137.138  9058.347  9169.358  9092.142
# 7   9583.349  9469.232  8983.291  8893.875  9339.560  8946.258
# 8   9604.084  9251.874  8894.268  8780.539  8958.188  8848.534
# 9   9624.810  9127.971  8756.540  8709.597  8840.023  8787.550
# 10  9538.400  9044.314  8724.225  8660.516  8723.664  8764.543
#          EEE       VEE       EVV       VVV
# 1  12730.596 12730.596 12730.596 12730.596
# 2   8458.926  8455.978  8080.928  8077.116
# 3   8307.503  8307.306  8110.035  8100.802
# 4   8216.067  8227.658  8146.920  8136.717
# 5   8158.428  8181.352  8177.743  8183.882
# 6   8139.580  8156.873  8201.515  8221.520
# 7   8136.885  8186.882  8243.219  8256.903
# 8   8155.261  8194.675  8279.479  8295.496
# 9   8175.112  8212.897  8323.535  8337.049
# 10  8189.729  8226.334  8347.291  8415.326

ans$topbic
#      row col                     
# [1,] "2" "EVV" "8080.92815908273"
# [2,] "2" "VVV" "8077.11578062968"
# [3,] "3" "VVV" "8100.80239791501"

## actually found correct model and cluster (2,VVV)!

## how does it look?
val <- norMmixMLE(x,2,trafo="clr1",model="VVV")

plot(val$norMmix)
points(x)
## good

## how about wrong model?
val <- norMmixMLE(x,2,trafo="clr1",model="EVV")

plot(val$norMmix)
# no still wrong

####
#-------------------------------------------------------------------------------
####
# work on 2019-08-04
####

## correction to above, it does work det() reveals that is is in fact of equal
## volume.
## seems the issue with unif alpha has been fixed
## tests:

x26 <- rnorMmix(n=1000, MW26)

val <- norMmixMLE(x26,2,trafo="clr1",model="EEI")
# det(Sig) are equal
plot(val$norMmix)


## try again the bic over models:

ans <- fit.norMmix(x26,k=1:10, models=1:10,trafo="clr1")
ans$topbic
#      row col                     
# [1,] "2" "EII" "9458.05706553518"
# [2,] "2" "VII" "9464.18990862962"
# [3,] "2" "EEI" "9465.60963158039"
## found correct one in top 3 (2, EEI)
MW26$model
# [1] "EEI"

## what if we grossly overestimate cluster number?
val <- norMmixMLE(x26,20,trafo="clr1", model="VVV")

## but abominably slow
## how does it look
matplot(ans$BIC, type="l")
## looks reasonably clear despite very overlapping clusters

## test speed

x26 <- rnorMmix(n=400, MW26)
## a bit smaller
timeMLE <- matrix(0,10,10)

for (i in 1:10) {
    for(j in 1:10) {
        timeMLE[i,j] <- as.numeric(system.time(norMmixMLE(x26,i,trafo="clr1", model=models[j]))[2])
    }
}

colMeans(timeMLE)
#  [1] 0.0084 0.0180 0.0132 0.0212 0.0220 0.0296 0.0128 0.0200
#  [9] 0.0308 0.0380
rowMeans(timeMLE)
#  [1] 0.0000 0.0016 0.0016 0.0040 0.0124 0.0220 0.0260 0.0304
#  [9] 0.0492 0.0668

## it seems models do not worsen time as much as cluster number
## optimization should start with cluster number
##


####
#-------------------------------------------------------------------------------
####
# work on 2019-08-05

## implemented use of mclusts init system using model based agglomerative 
## hierarchical clustering (MBAHC)
## can be chosen trough ini=c("cla", "mcl") in norMmixMLE()

## see what initial values are
x <- rnorMmix(1000, MW213)
invisible(norMmixMLE(x, 2, trafo="clr1", ini="cla", model="VVV", maxiter=10))
# initial  value 4007.604167 
# final  value 4007.604167 
# converged

invisible(norMmixMLE(x, 2, trafo="clr1", ini="mcl", model="VVV", maxiter=10))
# initial  value 4007.604167 
# final  value 4007.604167 
# converged

## seems more or less equal
## try false clusters:
norMmixMLE(x, 3, trafo="clr1", ini="cla", model="VVV", maxiter=10)
# initial  value 4027.656459 
# iter   5 value 4013.299059
# iter  10 value 4008.969466
# final  value 4008.969466 
# stopped after 10 iterations
norMmixMLE(x, 3, trafo="clr1", ini="mcl", model="VVV", maxiter=10)
# initial  value 4045.824343 
# iter   5 value 4017.915335
# iter  10 value 4009.437017
# final  value 4009.437017 
# stopped after 10 iterations

## clara seems better in this instance early on.

## after plotting 3 clusters for MW213(which has 2 clusters) once w/ mcl and 
## once w/ cla, it seems that there is a difference in how it handles the 
## third, superfluous, cluster.

## mclust does some strange clustering

x213 <- rnorMmix(1000, MW213)
mclclus <- hcVVV(x213)
indexmcl <- hclass(mclclus, 3)

plot(x213)
points(x213[indexmcl==1,], col="green")
points(x213[indexmcl==2,], col="blue")
points(x213[indexmcl==3,], col="red")



## subtract BIC of norMmixMLE and mclust

models <- c("EII","VII","EEI","VEI","EVI",
	    "VVI","EEE","VEE","EVV","VVV")

ansnMm <- fit.norMmix(x, k=1:10, models=1:10, ini="cla", trafo="clr1")

library(mclust)
ansmcl <- Mclust(x, G=1:10, modelNames=models)

ansnMm$BIC + ansmcl$BIC
#              EII           VII         EEI         VEI
# 1   1.469743e-09  1.469743e-09    6.907755    6.907755
# 2   3.819878e-11  3.427960e-06    6.907756    6.907755
# 3  -1.978644e+02 -2.133989e+01 -429.590749    6.563318
# 4  -1.247004e-01 -2.707077e-01    6.651546    6.539718
# 5  -9.796136e+01  2.317187e+01 -195.944574    6.437470
# 6   2.757753e+02  2.161082e+02    6.515021    6.425753
# 7   1.204404e+02  2.629785e+02    6.478995    6.536986
# 8  -2.473476e+01 -7.746697e+01 -251.231479 -143.361807
# 9  -8.762882e+01 -2.393603e+01 -116.271937  -94.064328
# 10  1.292093e+01  1.851451e+01  -45.518955   -1.816445
#            EVI        VVI        EEE        VEE       EVV
# 1     6.907755   6.907757   6.907755  6.9077553  6.907755
# 2    13.815516  13.815511   6.907755  6.9077527 13.815523
# 3  -426.204675 -66.076576 -46.424579 14.6932031 20.232303
# 4    27.223089  27.182020   6.730684  6.4667086 17.760054
# 5  -166.531539  40.934680 -28.577340  5.4530410 17.082762
# 6    41.076641  41.206814   6.564753  5.6334878 29.868263
# 7    32.475350 -29.953848 -23.814524 -8.1668977 25.241187
# 8  -258.580695 -51.968273 -27.531406  2.0798031 35.657169
# 9    61.911261  63.153852   5.710438  0.4673593 46.489965
# 10   68.798296  70.438876   2.676379 -6.1780501 69.015885
#          VVV
# 1   6.907755
# 2  13.815511
# 3  12.891298
# 4  21.838922
# 5  36.181956
# 6  33.913796
# 7  37.347069
# 8  48.423062
# 9  41.723388
# 10 45.575421

## if >0 => nMm > mclust and vice versa
## there seems to be a systematic errorconnected to the value 6.907755



diffbic <- ansnMm$BIC + ansmcl$BIC

mean(diffbic)
# [1] -7.588035

## mean would suggest that nMm is better, but a histogram shows, most diffs are
## a bit over 0
hist(c(diffbic))


ansnMmmcl <- fit.norMmix(x, k=1:10, models=1:10, ini="mcl", trafo="clr1")

#Error in nMm2par(obj = nMm.temp, trafo = trafo, model = model) :
#  4 diffbic <- ansnMm$BIC + ansmcl$BIC                                              |  isTRUE(all(apply(sig, 3, function(j) (ldl(j)$Diag >= 0)))) is not TR
#  5                                                                                 |UE




## different mixture, smaller size.
## does the systematic pattern reoccur?

x26 <- rnorMmix(500, MW26)

ansnMm2 <- fit.norMmix(x26, k=1:10, models=1:10, ini="cla", trafo="clr1")

ansmcl2 <- Mclust(x26, G=1:10, modelNames=models)

asdf <- ansnMm2$BIC + ansmcl2$BIC   # note mclust has bic defined as -BIC
# Bayesian Information Criterion (BIC): 
#              EII           VII        EEI          VEI
# 1   1.455192e-11  1.455192e-11  6.2146081   6.21460810
# 2  -1.395888e-01 -4.756255e-01  5.5012411   3.07846393
# 3  -1.323026e+00 -1.147379e+01  4.7637577  -5.20547432
# 4  -3.618632e+00 -6.123394e+00  5.2040117   0.66130297
# 5  -4.924523e+00 -1.169028e+01  3.3295968  -4.14854975
# 6  -2.824544e+00 -1.264269e+00  3.3806121  -5.46576008
# 7  -3.541775e+00  1.419451e+00  2.5362436   6.13049247
# 8  -5.041197e+00 -2.735047e-01  4.7643657  -0.05276408
# 9  -2.966653e+00 -3.971131e+01 -0.6441811  -9.15207041
# 10 -7.692870e+00 -1.122794e+01  2.8993839 -31.39121540
#          EVI       VVI       EEE        VEE       EVV       VVV
# 1   6.214608  6.214608  6.214608   6.214608  6.214608  6.214608
# 2  11.934082  9.631218  5.029121   2.663087 11.286666  8.834553
# 3  23.042092 14.411384  5.341183  -1.551446 13.557000 10.458182
# 4  19.177244  7.179716  5.398366 -16.677092 20.193401 -7.081574
# 5  23.924815 12.839222  5.434241   1.692398 27.991784 22.585312
# 6  34.881048 27.344113  3.512624  -2.375613 26.986966 25.253365
# 7  38.580564 37.169461 -4.700895  -1.186825 32.048934 20.212452
# 8  40.480655 38.245426  4.380756  -3.710059 46.525256 34.334146
# 9  48.278538 39.196200 -4.997565 -38.804163 53.815070 37.696185
# 10 62.531344 26.198212 -5.007481 -16.351691 56.505000 55.631573

hist(c(asdf))
## still skews to the right of zero
mean(asdf)
# [1] 8.729814

## pattern of equal values has disappeared for all but cluster=1
## no easily recognizable pattern
## possibly still faulty calculations in nMm functions
## histogram indicates, that mclust is still better, giving lower BIC
## what parameter counts does mclust use?
## 

ansnMmmcl <- fit.norMmix(x26, k=1:5, models=1:10, ini="mcl", trafo="clr1")

## lower cluster number for time reasons

asas <- ansnMmmcl$BIC + ansmcl2$BIC[1:5,]
#             EII           VII      EEI        VEI       EVI
# 1  1.455192e-11  1.455192e-11 6.214608   6.214608  6.214608
# 2 -1.395893e-01 -4.756254e-01 5.501241   3.078465 11.934082
# 3 -1.323177e+00 -5.898251e-01 4.763753   7.082149 15.931463
# 4 -3.618634e+00 -1.052905e+01 5.204012  -3.010698 19.164103
# 5 -3.878573e+00 -4.649895e+01 3.329679 -21.178059 25.452601
#         VVI      EEE        VEE       EVV       VVV
# 1  6.214608 6.214608   6.214608  6.214608  6.214608
# 2  9.631218 5.029163   2.663085 11.286669  8.834554
# 3 14.411384 5.341183   7.238553 13.556935 10.458178
# 4 10.827287 5.398480  -5.686936 18.098031  1.836237
# 5 12.682385 2.607073 -42.152267 28.450389 25.378072

hist(c(asas))
mean(asas)
# [1] 4.116118

## seems finding 1 cluster solution is very stable.
## always gives result between 6,7.
## still, mclust seems marginally better, and a _lot_ faster



####
#-------------------------------------------------------------------------------
####
# work on 2019-08-07
####

## speeding up llnorMmix

## setup

x <- rnorMmix(n=5000,MW26)

clus <- cluster::clara(x=x, 20, rngR=T, pamLike=T, samples=10)
index <- clus$clustering
tau <- matrix(0,5000,20)
tau[cbind(1:5000,index)] <- 1

nMm.temp <- mstep.nMm(x, tau)
# create par. vector out of m-step
initpar. <- nMm2par(obj=nMm.temp, trafo="clr1", model="VVV")


bench <- function() {
    system.time( {
        for (i in 1:100) {
            llnorMmix(initpar., x, 2, 20, trafo="clr1", model="VVV")
# [1] -23856.56
        }
    })
}

bench()
#    user  system elapsed 
#   0.564   0.000   0.562 

D. <- matrix(runif(40), 2, 20)
alpha <- runif(20)

## this should be faster
system.time({for (i in 1:1000) {invalpha <- 1/exp(D.+rep(alpha,each=2))}})
#    user  system elapsed 
#   0.004   0.000   0.006 


system.time({for (i in 1:1000){for (j in 1:20){invalpha <- 1/exp(alpha[j]+D.[,j])}}})
#    user  system elapsed 
#   0.020   0.000   0.018 

bench()
#    user  system elapsed 
#   0.564   0.000   0.563 
#    user  system elapsed 
#   0.568   0.000   0.570 
#    user  system elapsed 
#   0.744   0.000   0.746 

## first execution seems slower than subsequent exs??



## make substitution of assignment of D.
## D. <- matrix(....
## D. <- apply(D.,2, function(j) j-sum(j)/p)


bench()
#    user  system elapsed 
#   0.588   0.000   0.597 
#    user  system elapsed 
#   0.580   0.000   0.585 
#    user  system elapsed 
#   0.576   0.000   0.576 
#    user  system elapsed 
#   0.856   0.004   0.868 

## seems also slower

bench()
#    user  system elapsed 
#   0.572   0.000   0.574 
#    user  system elapsed 
#   0.568   0.000   0.569 
#    user  system elapsed 
#   0.732   0.000   0.735 

## as untangling code seems to be making things faster, I try to make 
## j-sum(j)/p its own function

bench()
#    user  system elapsed 
#   0.572   0.000   0.582 
#    user  system elapsed 
#   0.576   0.000   0.577 
#    user  system elapsed 
#   0.576   0.004   0.580 
#    user  system elapsed 
#   0.580   0.000   0.578 
#    user  system elapsed 
#   0.724   0.008   0.732 

## doesn't work
## changed back

bench()
#    user  system elapsed 
#   0.568   0.000   0.566 
#    user  system elapsed 
#   0.564   0.004   0.568 
#    user  system elapsed 
#   0.568   0.000   0.565 
#    user  system elapsed 
#    0.56    0.00    0.56 
#    user  system elapsed 
#   0.704   0.000   0.707 
#    user  system elapsed 
#    0.58    0.00    0.58 

## about the same

## trying preallocating invl
ans <- bench()
#    user  system elapsed 
#   0.584   0.000   0.585 
#    user  system elapsed 
#   0.580   0.004   0.593 
#    user  system elapsed 
#   0.596   0.000   0.597 
#    user  system elapsed 
#   0.728   0.000   0.728 

## no


## change bench function for better visibility of speed


bench <- function() {
    times <- vector("numeric", length(30:100))
    for (k in 30:100) {
        times[k-29] <- system.time( {
            for (i in 1:k) {
                llnorMmix(initpar., x, 2, 20, trafo="clr1", model="VVV")
                # [1] -23856.56
            }
        })[[3]]
        print(times[k-29])
    }
    return(times)
}


tre <- bench()
plot(tre, type="l")

## from now on tre is the mark to beat



## can take L. <- diag(1,p) out of the loop, needs to be done only once

rew1 <- bench()
points(rew1, type="l", col="green")
## about the same

## try again to predo invalpha

rew2 <- bench()
points(rew2, type="l", col="red")
## now we see it is about the same

## weird spikes


## try bytecomplile the function

library(compiler)

llnorcmp <- cmpfun(llnorMmix)


benchcmp <- function() {
    times <- vector("numeric", length(30:100))
    for (k in 30:100) {
        times[k-29] <- system.time( {
            for (i in 1:k) {
                llnorcmp(initpar., x, 2, 20, trafo="clr1", model="VVV")
                # [1] -23856.56
            }
        })[[3]]
        print(times[k-29])
    }
    return(times)
}


rew3 <- benchcmp()
points(rew3, type="l", col="blue")
## didn't help


## tried make function out of backsolve part

rew4 <- bench()
points(rew4, type="l", col="gray")

## does nothing


## calculations

tables <- cbind(30:100,tre)
df <- data.frame(tables)
colnames(df) <- c("iterations", "time")
lmfit <- lm(time~iterations, df)
lmfit
# 
# Call:
# lm(formula = time ~ iterations, data = df)
# 
# Coefficients:
# (Intercept)   iterations  
#   -0.006563     0.005994  
# 


## after rm(list=ls()) somehow runs faster, does R save functions into local 
## memory???

## turns out I need to re-source bench() to get it to recognize llnorMmix
## can redo all tests

## substitution of assignment

rew5 <- bench()
points(rew5, type="l", col="purple")

# maybe better? hard to tell, probably irrelevant


## try byte compile again

rew6 <- benchcmp()
points(rew6, type="l", col="orange")

## not significantly better




## again quit R console and restarted it, all of a sudden numbers are
## improved

## now R starts with --vanilla argument, maybe that helps.

## writing benchmark as it is now

tre <- bench()
tre
#  [1] 0.333 0.179 0.180 0.190 0.195 0.206 0.216 0.214 0.221
# [10] 0.226 0.231 0.235 0.245 0.246 0.257 0.263 0.263 0.269
# [19] 0.278 0.281 0.326 0.294 0.300 0.303 0.307 0.317 0.319
# [28] 0.321 0.324 0.333 0.339 0.347 0.360 0.358 0.366 0.366
# [37] 0.374 0.392 0.387 0.436 0.399 0.399 0.410 0.421 0.417
# [46] 0.422 0.428 0.440 0.437 0.450 0.450 0.466 0.461 0.464
# [55] 0.517 0.481 0.482 0.487 0.502 0.496 0.508 0.517 0.538
# [64] 0.521 0.533 0.533 0.576 0.546 0.549 0.548 0.564
plot(tre, type="l")

## this is a bit better than it was this morning



system.time(for (l in 1:100) alpha <- par.[f:f2] )
#    user  system elapsed 
#   0.004   0.000   0.002 
system.time(for (l in 1:100) D. <- apply(matrix(par.[f2.1:f22],p,k),2, function(j) j-sum(j)/p) )
#    user  system elapsed 
#   0.012   0.000   0.010 
system.time(for (l in 1:100) invalpha <- (1/exp(rep(alpha, each=p)+D.)))
#    user  system elapsed 
#   0.004   0.000   0.004 
system.time(for (l in 1:100) L.temp <- matrix(par.[f22.1:f221],p*(p-1)/2,k))
#    user  system elapsed 
#   0.004   0.000   0.004 
system.time(for (l in 1:100) L. <- diag(1,p))
#    user  system elapsed 
#   0.000   0.000   0.002 
for (i in 1:k) {
 system.time(for (l in 1:100)    L.[lower.tri(L., diag=FALSE)] <- L.temp[,i])
#    user  system elapsed 
#   0.000   0.000   0.003 
 system.time(for (l in 1:100)    rss <- colSums(invalpha[,i]*backsolve(L., (x-mu[,i]), upper.tri=FALSE)^2))
#    user  system elapsed 
#   0.008   0.000   0.007 
 system.time(for (l in 1:100)    invl <- invl+w[i]*exp(-0.5*(p*(alpha[i]+l2pi)+rss)))
#    user  system elapsed 
#   0.004   0.000   0.006 
}
system.time(for (l in 1:100) sum(log(invl)))
#    user  system elapsed 
#   0.004   0.000   0.002 

##



        alpha <- par.[f:f2]
        D. <- apply(matrix(par.[f2.1:f22],p,k),2, function(j) j-sum(j)/p)
        invalpha <- (1/exp(rep(alpha, each=p)+D.))
        L.temp <- matrix(par.[f22.1:f221],p*(p-1)/2,k)
        L. <- diag(1,p)
system.time(
        for (l in 1:100) {
            for (i in 1:k) {
            L.[lower.tri(L., diag=FALSE)] <- L.temp[,i]
            rss <- colSums(invalpha[,i]*backsolve(L., (x-mu[,i]), upper.tri=FALSE)^2)
            invl <- invl+w[i]*exp(-0.5*(p*(alpha[i]+l2pi)+rss))
        }})

# user  system elapsed
# 0.100   0.012   0.111
        sum(log(invl))



system.time( 

    for (l in 1:100) {

        alpha <- par.[f:f2]
        D. <- apply(matrix(par.[f2.1:f22],p,k),2, function(j) j-sum(j)/p)
        invalpha <- (1/exp(rep(alpha, each=p)+D.))
        L.temp <- matrix(par.[f22.1:f221],p*(p-1)/2,k)
        L. <- diag(1,p)
        for (i in 1:k) {
            L.[lower.tri(L., diag=FALSE)] <- L.temp[,i]
            rss <- colSums(invalpha[,i]*backsolve(L., (x-mu[,i]), upper.tri=FALSE)^2)
            invl <- invl+w[i]*exp(-0.5*(p*(alpha[i]+l2pi)+rss))
        }
        sum(log(invl))

    })


## try Rprof

Rprof(file="profile.out", line.profiling=TRUE)
        alpha <- par.[f:f2]
        D. <- apply(matrix(par.[f2.1:f22],p,k),2, function(j) j-sum(j)/p)
        invalpha <- (1/exp(rep(alpha, each=p)+D.))
        L.temp <- matrix(par.[f22.1:f221],p*(p-1)/2,k)
        L. <- diag(1,p)
        for (i in 1:k) {
            L.[lower.tri(L., diag=FALSE)] <- L.temp[,i]
            rss <- colSums(invalpha[,i]*backsolve(L., (x-mu[,i]), upper.tri=FALSE)^2)
            invl <- invl+w[i]*exp(-0.5*(p*(alpha[i]+l2pi)+rss))
        }
        sum(log(invl))

Rprof(NULL)

summaryRprof("profile.out", lines = "show")


Rprof(file="profile.out", line.profiling="TRUE")
eval(parse(file="~/ethz/BA/norMmix/R/llnorMmix.R", keep.source=TRUE))
Rprof(NULL)
summaryRprof("profile.out", lines= "show")


## no clue what to do.
## e-mail to Mächler??
## tomorrows problem.


####
#-------------------------------------------------------------------------------
####
# work on 2019-08-09
####
# implement smaller parameters by using fewer params in D.


## finally tests don't return errors
## this has to be an error


####
#-------------------------------------------------------------------------------
####
# work on 2019-08-10
####


## fixed unrelated issue and all of a sudden all bugs are gone???

## not gonna complain

## random test as sanity check

ansnmm <- fit.norMmix(x, k= 1:5, models=1:10,trafo="clr1", ini="cla")

ansmcl <- mclust::Mclust(x, G=1:5, modelNames=models)

diffbic <- ansnmm$BIC --ansmcl$BIC
# Bayesian Information Criterion (BIC): 
#             EII           VII        EEI        VEI        EVI
# 1  2.910383e-10  2.910383e-10   0.000000   0.000000   0.000000
# 2 -1.544405e+00 -4.805815e+00  -1.815754  -5.152510  -1.756925
# 3 -1.280272e+01 -7.160097e+00 -15.053828  -6.772972 -10.997506
# 4 -9.986395e+00 -5.756005e+00 -15.835385 -12.186757  -7.619998
# 5 -8.631051e+00 -1.526588e+01 -14.627508 -13.272501 -13.462888
#             VVI           EEE           VEE           EVV
# 1  5.239697e-04  2.182787e-11  2.182787e-11  2.182787e-11
# 2 -4.955601e+00 -1.918936e+00 -5.843762e+00 -4.231558e+00
# 3 -9.309414e+00 -4.802646e+00 -3.740256e+00 -5.803574e+00
# 4 -1.773906e+01 -8.903814e+00 -8.377495e+00 -1.135212e+01
# 5 -1.829755e+01 -9.935866e+00 -1.106693e+01 -2.856427e+01
#             VVV
# 1  2.182787e-11
# 2 -1.021222e+01
# 3 -5.392810e+00
# 4 -1.773004e+01
# 5 -2.383822e+01

## really good now

mean(diffbic)
# [1] -7.93045
hist(diffbic)

colSums(diffbic)
#       EII       VII       EEI       VEI       EVI       VVI 
# -32.96457 -32.98780 -47.33248 -37.38474 -33.83732 -50.30110 
#       EEE       VEE       EVV       VVV 
# -25.56126 -29.02844 -49.95153 -57.17328 
rowSums(diffbic)
#             1             2             3             4 
#  5.239703e-04 -4.223748e+01 -8.183583e+01 -1.154871e+02 
#             5 
# -1.569627e+02 
colMeans(diffbic)
#        EII        VII        EEI        VEI        EVI 
#  -6.592915  -6.597560  -9.466495  -7.476948  -6.767463 
#        VVI        EEE        VEE        EVV        VVV 
# -10.060220  -5.112252  -5.805688  -9.990306 -11.434656 
rowMeans(diffbic)
#             1             2             3             4 
#  5.239703e-05 -4.223748e+00 -8.183583e+00 -1.154871e+01 
#             5 
# -1.569627e+01 

## big improvement, seems strictly better in all but k=1 cases.
## no longer weird fluctuations in values, seems proportional to k.
## not as much proportional to models.

## TODO
## figure out what stat tests to use. maybe relative likelihood or something
## finish plot function
## larger cases p>2 k>3

## out of curiosity
hist(c(ansnmm$BIC[-1,], -ansmcl$BIC[-1,]))

plot(sort(ansmcl$BIC))
plot(sort(ansnmm$BIC), col="red")
points(sort(-ansmcl$BIC),col="blue")


####
#-------------------------------------------------------------------------------
####
## work on 2019-08-11
####

## fix up plot function
## uses plot functionality



## started plotnd.norMmix,  doesn't work yet, but close


####
#-------------------------------------------------------------------------------
####
## work on 2019-08-12
####

## today, work on larger cases and trying to finish plot fctn
## maybe see if mcl still gives errors

## wrote some cases with p=3,5 and plot function works now for p>2

x3 <- rnorMmix(2000, MW33)
rnmm <- fit.norMmix(x3, k=1:5, models=1:10, ini="cla", trafo="clr1")

rmcl <- mclust::Mclust(x3, G=1:5, modelNames=models)

diffbic <- rnmm$BIC --rmcl$BIC
# Bayesian Information Criterion (BIC): 
#             EII           VII         EEI        VEI        EVI
# 1  8.731149e-11  8.731149e-11   0.0000000   0.000000   0.000000
# 2 -2.766107e-01 -6.967244e-01  -0.8540788  -3.071501  -2.911382
# 3 -1.587003e+00 -1.729446e+00  -5.8811812  -9.260836  -4.150406
# 4 -5.631269e+00 -2.476995e+00 -10.7263792  -8.995189 -23.402054
# 5 -2.579224e+01 -5.149856e+00 -20.1837536 -29.090563  -6.702159
#             VVI           EEE           VEE           EVV
# 1  7.930794e-10 -1.455192e-11 -1.455192e-11 -1.455192e-11
# 2 -2.801682e+00 -1.199729e+00 -3.503617e+00 -1.805226e+00
# 3 -8.893355e+00 -7.843882e+00 -1.579852e+01 -1.752565e+01
# 4 -1.968783e+01 -9.987366e+00 -3.088096e+01 -2.535595e+01
# 5 -1.926292e+01 -1.415199e+01  1.598953e+01 -2.495358e+01
#             VVV
# 1 -1.455192e-11
# 2 -6.335436e+00
# 3 -1.309323e+01
# 4 -2.410534e+01
# 5 -4.285166e+01


tmp <- fit.norMmix(x3, k=1:10, model=1:10, ini="mcl", trafo="clr1", maxiter=6)
tmp

## "mcl" issue could not be reproduced here, possibly other mixture is ill 
## conditioned


####
#-------------------------------------------------------------------------------
####
## work on 2019-08-16
####

## establish connection to sfs servers

## maybe fix up rnorMmix a bit
## now rnorMmix supports cluster index and sampling

####
#-------------------------------------------------------------------------------
####
## work on 2019-08-17
####

## test saving objects

x <- rnorMmix(1000, MW26, index=TRUE, sampling=TRUE)


saveRDS(x, file="MW26e3")

x26 <- readRDS("MW26e3")

mods <- list(MW210,MW213,MW28,MW22,MW25,MW33,MW211,MW26,MW29,MW23,MW31,MW34,MW212,MW27,MW21,MW24,MW32,MW51)

rnorMmix(12, mods[[3]])

set.seed(2019)

asdf <- function(j) {
    set.seed(2019)
    x <- rnorMmix(100, j, index=T, sampling=T)
    saveRDS(x, file=paste0("mw.1e2.", j$dim, j$k, j$model,".RDS"))
}

lapply(mods, asdf)


xx <- rnorMmix(10000, MW34)

ans <- fit.norMmix(xx, k=1:10, models=1:10, trafo="clr1", ini="cla")

####
#-------------------------------------------------------------------------------
####
## work on 2019-08-17
####

## get readRDS writeRDS to work with norMmix

## list.files?

aa <- list.files("~/ethz/BA/Rscripts")

## regex?

aa <- list.files("~/ethz/BA/Rscripts", pattern="mw.*")

## yes, close to gettin it right



aa <- list.files("~/ethz/BA/Rscripts", pattern="mw.*")

for (i in aa) {
    x <- readRDS(paste0("~/ethz/BA/Rscripts/",i))
    print(x[1,])
}

## works!!

## now write file Auswertung-fit-mw.R


## see if merging worked

x <- rnorMmix(100, MW26)

ans <- fit.norMmix(x, k=1:10, models=1:10, trafo="clr1", ini="clara")

## important things work


####
##------------------------------------------------------------------------------
####
## work on 2019-08-20

## llnorMmix wasn't used for norMmixMLE until now and now its broken!!

x <- rnorMmix(100, MW26)

ans <- norMmixMLE(x,3,trafo="clr1", model="EII", ini="clara")

## non-finite finite-difference value

pp <- nMm2par(MW23, trafo="clr1", model="VVV")

llnorMmix(pp, x, 2, trafo="clr1", model="VVV")


####
##------------------------------------------------------------------------------
####
## work on 2019-08-21


## fixed error with llnorMmix
## problem was the  (if(weights>1) return inf)  statement
## commented that out and suddenly it doesnt crash everytime

## see if it is stable


x <- rnorMmix(100, MW26)

ans <- fit.norMmix(x, k=1:10, models=1:10, trafo="clr1", ini="clara")

aa <- BIC(ans)
bb <- logLik(ans)

mle <- norMmixMLE(x,k=1, model="EII", trafo="clr1", ini="clara")

## testing tests

asdf <- fit.norMmix(x, k=1:4, model=1:10, trafo="clr1", ini="clara",ll="nmm")


## larger numbers:

(MWdat <- Filter(function(.) is.norMmix(get(., "package:norMmix")),
                 ls("package:norMmix", pattern = "^MW[1-9]")))

ret <- list()

for (i in MWdat) {

    nm <- get(i, "package:norMmix")

    set.seed(2019); x <- rnorMmix(100, nm)

    aa <- fit.norMmix(x, k=1:4, model=1:10, trafo="clr1", ini="clara",
                      ll="nmm")

    bb <- fit.norMmix(x, k=1:4, model=1:10, trafo="clr1", ini="clara",
                      ll="mvt")

    a <- BIC(aa)[[1]]
    b <- BIC(bb)[[1]]

    ret[[i]] <- a-b

}


for (i in MWdat) {print(max(abs(c(ret[[i]]))))}
# [1] 9.71454e-09
# [1] 0.9334958
# [1] 0.2815176
# [1] 3.326458
# [1] 6.264971e-05
# [1] 9.390533e-11
# [1] 2.036131e-10
# [1] 0.6392211
# [1] 0.6392211
# [1] 1.947797e-09
# [1] 1.748394e-05
# [1] 1.21986e-09
# [1] 1.283524e-10
# [1] 4.729372e-11
# [1] 9.409709e-07
# [1] 8.458301e-11
# [1] 2.711181
# [1] NA
# NULL

for (i in MWdat) {
    print(ret[[i]])
    readline(prompt="enter to cont")
}

## seems good but some weird values, maybe llnmm gives some weird retvals.
## check manually

retnmm <- list()

for (i in MWdat) {

    nm <- get(i, "package:norMmix")

    set.seed(2019); x <- rnorMmix(100, nm)

    aa <- fit.norMmix(x, k=1:4, model=1:10, trafo="clr1", ini="clara",
                      ll="nmm")

    retnmm[[i]] <- BIC(aa)
}


for (i in MWdat) {
    print(retnmm[[i]])
    readline(prompt="enter to cont")
}

## ah, I see an issue. deleted value <- -value in MLE
## change in logLik


uio <- fit.norMmix(x, k=1:4, model=1:10, trafo="clr1", ini="clara", ll="nmm")

BIC(uio)


## done, put it into a commit

ij <- norMmixMLE(x, k=3, model="EII", trafo="clr1", ini="clara", ll="nmm")

plot(ij$norMmix)
ij$norMmix$mu
#          [,1]      [,2]     [,3]
# [1,] 10.95137  6.267234 1.122263
# [2,] 12.04606  6.907823 1.702433
# [3,] 12.91923  8.034070 3.433608
# [4,] 14.02475  8.738717 3.905640
# [5,] 14.98502 10.057415 4.543396
MW51$mu
#      [,1] [,2] [,3]
# [1,]    1    6   11
# [2,]    2    7   12
# [3,]    3    8   13
# [4,]    4    9   14
# [5,]    5   10   15

## reasonable

x <- rnorMmix(1000, MW51)
ij <- norMmixMLE(x, k=3, model="EII", trafo="clr1", ini="clara", ll="nmm")

ij$norMmix$mu
#          [,1]      [,2]     [,3]
# [1,] 10.97616 0.9919569 6.043642
# [2,] 11.97297 1.9991118 7.003464
# [3,] 12.98681 2.9527947 8.042657
# [4,] 13.99832 4.0824587 9.026949
# [5,] 15.08960 5.0045929 9.900400
MW51$mu
#      [,1] [,2] [,3]
# [1,]    1    6   11
# [2,]    2    7   12
# [3,]    3    8   13
# [4,]    4    9   14
# [5,]    5   10   15


x <- rnorMmix(10000, MW51)
ij <- norMmixMLE(x, k=3, model="EII", trafo="clr1", ini="clara", ll="nmm")

ij$norMmix$mu
#           [,1]     [,2]     [,3]
# [1,] 0.9916475 10.97603 5.997823
# [2,] 1.9929993 11.97723 6.977599
# [3,] 2.9983107 13.02145 8.021266
# [4,] 3.9778741 14.01703 8.974741
# [5,] 4.9951813 15.00860 9.998929


## expected improvement in accuracy


## time improvement?

t1 <- system.time(norMmixMLE(x, k=3, model="EII", trafo="clr1", ini="clara", ll="nmm"))
t2 <- system.time(norMmixMLE(x, k=3, model="EII", trafo="clr1", ini="clara", ll="mvt"))

t1
#    user  system elapsed 
#   3.764   0.128   3.895 
t2
#    user  system elapsed 
#   5.380   0.176   5.563 
t2/t1
#     user   system  elapsed 
# 1.429330 1.375000 1.428241 

## after removing x <- t(x)
t2/t1
#     user   system  elapsed 
# 1.483402 1.928571 1.496593 


## Auswertung-fit-various given to ada server

## write analysis tools for various data sets


####
##------------------------------------------------------------------------------
####
## work on 2019-08-22

## ada server done analysing data sets:

## results:
## smi: failed for k>=5, models=EII,VII,VVI
## loss, iris, crashed because is.numeric failed

data(SMI.12, package="copula")

smi <- SMI.12

ans <- norMmixMLE(smi,k=1, model="EII", trafo="clr1", ini="clara", ll="nmm")
#Error in optim(initpar., neglogl, method = method, control = control) :
#initial value in 'vmmin' is not finite
#Error in inherits(ok, "try-error") : object 'ok' not found

## no clue what that means

ans <- norMmixMLE(smi,k=7, model="EVI", trafo="clr1", ini="clara", ll="nmm")
## D!>=0

## need to remove the check for non degenerate case

iri <- iris[-5]
ans <- norMmixMLE(iri,k=1, model="EII", trafo="clr1", ini="clara", ll="nmm")
## is ok
ans <- norMmixMLE(iri,k=4, model="EVI", trafo="clr1", ini="clara", ll="nmm")

data(loss, package="copula")

ans <- norMmixMLE(loss,k=4, model="EVI", trafo="clr1", ini="clara", ll="mvt")

## issues with non-finite finite-difference value, with both nmm and mvt

ret <- fit.norMmix(iri, k=1:4, model=1:10, trafo="clr1", ini="clara", ll="nmm")


####
##------------------------------------------------------------------------------
####
## work on 2019-08-23

## see if I can't speed up llnorMmix

x <- rnorMmix(5000, MW29)

models <- c("EII","VII","EEI","VEI","EVI",
	    "VVI","EEE","VEE","EVV","VVV")

retnmm <- list()
retmvt <- list()

for (i in models) {
    pars <- nMm2par(MW29, trafo="clr1", model=i)
    retnmm[[i]] <- system.time(llnorMmix(pars, t(x), 2, trafo="clr1", model=i))
}

for (i in models) {
    pars <- nMm2par(MW29, trafo="clr1", model=i)
    retmvt[[i]] <- system.time(llmvtnorm(pars, x, 2, trafo="clr1", model=i))
}

## no noticeable difference, more loops



for (i in models) {
    pars <- nMm2par(MW29, trafo="clr1", model=i)
    retnmm[[i]] <- system.time(
        for (j in 1:100) {llnorMmix(pars, t(x), 2, trafo="clr1", model=i)})
}

for (i in models) {
    pars <- nMm2par(MW29, trafo="clr1", model=i)
    retmvt[[i]] <- system.time(
        for (j in 1:100) {llmvtnorm(pars, x, 2, trafo="clr1", model=i)})
}

## llnorMmix strictly better!!

ret <- list()
for (i in models) {
    ret[[i]] <- retmvt[[i]]/retnmm[[i]]
}
ret
# $EII
#     user   system  elapsed 
# 1.692308      NaN 1.784314 
# 
# $VII
#     user   system  elapsed 
# 1.769231      NaN 1.840000 
# 
# $EEI
#     user   system  elapsed 
# 1.571429      NaN 1.750000 
# 
# $VEI
#     user   system  elapsed 
# 1.923077      NaN 2.000000 
# 
# $EVI
#     user   system  elapsed 
# 1.785714      NaN 1.771930 
# 
# $VVI
#     user   system  elapsed 
# 1.866667      NaN 1.931034 
# 
# $EEE
#    user  system elapsed 
# 1.56250     NaN 1.59375 
# 
# $VEE
#     user   system  elapsed 
# 1.562500      NaN 1.584615 
# 
# $EVV
#     user   system  elapsed 
# 1.555556      NaN 1.513514 
# 
# $VVV
#     user   system  elapsed 
# 1.473684      NaN 1.473684 
# 

a <- lapply(ret, function(j) j[[1]] )
a <- unlist(a, use.names=FALSE)
a
#  [1] 1.692308 1.769231 1.571429 1.923077 1.785714 1.866667
#  [7] 1.562500 1.562500 1.555556 1.473684
mean(a)
# [1] 1.676266
sqrt(var(a))
# [1] 0.1529464


## look again into error in smi.12

data(SMI.12, package="copula")
smi <- SMI.12
ans <- norMmixMLE(smi,k=1, model="EII", trafo="clr1", ini="clara", ll="nmm")



####
##------------------------------------------------------------------------------
####
## work on 2019-08-23


## forcing positive definite matrix

d <- 1:4
L <- matrix(c(1,2,3,4,0,1,2,1,0,0,1,4,0,0,0,1),4,4)

A <- L%*% diag(d) %*% t(L)

x <- rnorMmix(1000, MW51)

ans <- norMmixMLE(x,3,model="VVI", trafo="clr1", ini="clara", ll="nmm")

## wrote forcePositive and inserted into MLE
ans <- norMmixMLE(x,3,model="VVI", trafo="clr1", ini="clara", ll="nmm")

## test with known 
data(SMI.12, package="copula")
smi <- SMI.12
ans <- norMmixMLE(smi,k=7, model="EVI", trafo="clr1", ini="clara", ll="nmm")

## doesn't work...

ds <- .Machine$double.eps
# [1] 2.220446e-16

## at p~20 numerical error is large enough to become non-negligible
log2(.Machine$double.eps)
# [1] -52

## sleep on this


## analyse fit-various

fv_dir <- normalizePath("~/ethz/BA/fit-various")

data(SMI.12, package="copula")
smi <- SMI.12

library(mclust)
models <- c("EII","VII","EEI","VEI","EVI",
	    "VVI","EEE","VEE","EVV","VVV")
mcl <- Mclust(smi, G=1:9, modelNames=models)
mcl$BIC
# Bayesian Information Criterion (BIC): 
#         EII       VII       EEI       VEI       EVI       VVI
# 1 -27040.40 -27040.40 -18221.49 -18221.49 -18221.49 -18221.49
# 2 -23616.88 -23210.41 -15963.08 -15855.59 -15804.29 -15706.51
# 3 -21993.92 -21779.22 -15093.34 -15026.02 -15020.48 -14930.06
# 4 -21145.54 -21020.72 -14449.88 -14612.14 -14333.75 -14457.47
# 5 -20451.94 -20132.82 -14183.27 -14169.27 -14158.41 -14112.15
# 6 -20374.68 -19867.57 -13902.63 -13920.47 -14259.14 -13974.76
# 7 -19936.92 -20074.01 -13591.04 -13663.24 -13963.23 -13818.66
# 8 -19962.27 -19592.58 -13502.57 -13323.99 -14020.86 -13614.06
# 9 -19479.90 -19429.46 -13168.27 -13033.79 -13504.23 -13359.90
#         EEE       VEE       EVV       VVV
# 1 -12246.64 -12246.64 -12246.64 -12246.64
# 2 -12113.03 -12108.64 -12025.43 -11805.19
# 3 -11996.94        NA        NA        NA
# 4 -11801.07        NA        NA        NA
# 5 -11640.26        NA        NA        NA
# 6 -11520.27        NA        NA        NA
# 7 -11494.28        NA        NA        NA
# 8 -11389.07        NA        NA        NA
# 9 -11322.15        NA        NA        NA
# 
# Top 3 models based on the BIC criterion: 
#     EEE,9     EEE,8     EEE,7 
# -11322.15 -11389.07 -11494.28 

## also has NAs for m=VEE,EVV,VVV k>2
## possibly because it is ill conditioned?

parcond(smi, k=3, model="VVV")
# [1] 0.2037572
## cutoff
parcond(smi, k=2, model="VVV")
# [1] 0.3058568

## how did MLE do?

nmm <- readRDS(file=file.path(fv_dir, "fit_smi_clr1_clara_nmm.rds"))

BIC(nmm$fit)
# [[1]]
#   EII VII      EEI      VEI      EVI VVI      EEE      VEE
# 1  NA  NA 18221.49 18221.49 18221.49  NA 12246.64 12246.64
# 2  NA  NA 15963.10 15893.15 15804.21  NA 12142.69 12132.12
# 3  NA  NA 15146.46 15085.00 14912.88  NA 12060.13 12031.97
# 4  NA  NA 14449.88 14421.11 14353.18  NA 11796.38 11843.84
# 5  NA  NA       NA       NA       NA  NA       NA       NA
# 6  NA  NA       NA       NA       NA  NA       NA       NA
# 7  NA  NA       NA       NA       NA  NA       NA       NA
# 8  NA  NA       NA       NA       NA  NA       NA       NA
# 9  NA  NA       NA       NA       NA  NA       NA       NA
#        EVV      VVV
# 1 12246.64 12246.64
# 2 11844.60 11698.89
# 3 11882.80 11872.70
# 4 12214.68 12151.10
# 5       NA       NA
# 6       NA       NA
# 7       NA       NA
# 8       NA       NA
# 9       NA       NA
# 
# $best
# [1] "2"   "VVV"
# 

BIC(nmm$fit)[[1]] + mcl$BIC
# Bayesian Information Criterion (BIC): 
#   EII VII          EEI        VEI           EVI VVI       EEE
# 1  NA  NA  0.000000000    0.00000    0.00000000  NA  0.000000
# 2  NA  NA  0.016141443   37.56204   -0.08607942  NA 29.654752
# 3  NA  NA 53.115187442   58.98080 -107.59887927  NA 63.191110
# 4  NA  NA  0.001328303 -191.03100   19.42585791  NA -4.696372
# 5  NA  NA           NA         NA            NA  NA        NA
# 6  NA  NA           NA         NA            NA  NA        NA
# 7  NA  NA           NA         NA            NA  NA        NA
# 8  NA  NA           NA         NA            NA  NA        NA
# 9  NA  NA           NA         NA            NA  NA        NA
#       VEE       EVV           VVV
# 1  0.0000    0.0000  1.818989e-12
# 2 23.4841 -180.8275 -1.063070e+02
# 3      NA        NA            NA
# 4      NA        NA            NA
# 5      NA        NA            NA
# 6      NA        NA            NA
# 7      NA        NA            NA
# 8      NA        NA            NA
# 9      NA        NA            NA


## varies wildly, not surprising as parcond is terrible


nmmm <- readRDS(file=file.path(fv_dir, "fit_smi_clr1_mclVVV_nmm.rds"))

BIC(nmmm$fit)[[1]] + mcl$BIC
# Bayesian Information Criterion (BIC): 
#   EII VII          EEI        VEI          EVI VVI       EEE
# 1  NA  NA 0.000000e+00    0.00000    0.0000000  NA   0.00000
# 2  NA  NA 3.896003e-04   37.56039   -0.1053436  NA  17.82603
# 3  NA  NA 7.330341e+01   58.98109 -107.5995117  NA 103.11319
# 4  NA  NA 3.985540e-04 -191.03112  100.5534188  NA 126.96456
# 5  NA  NA           NA         NA           NA  NA        NA
# 6  NA  NA           NA         NA           NA  NA        NA
# 7  NA  NA           NA         NA           NA  NA        NA
# 8  NA  NA           NA         NA           NA  NA        NA
# 9  NA  NA           NA         NA           NA  NA        NA
#        VEE      EVV          VVV
# 1  0.00000   0.0000 1.818989e-12
# 2 72.64549 147.8577 3.741246e+02
# 3       NA       NA           NA
# 4       NA       NA           NA
# 5       NA       NA           NA
# 6       NA       NA           NA
# 7       NA       NA           NA
# 8       NA       NA           NA
# 9       NA       NA           NA
# 
# Top 3 models based on the BIC criterion: 
#    VVV,2    EVV,2    EEE,4 
# 374.1246 147.8577 126.9646 




data(loss , package="copula")
mcll <- Mclust(loss, G=1:9, modelNames=models)

nmml <- readRDS(file=file.path(fv_dir, "fit_los_clr1_clara_nmm.rds"))

mcll$BIC
# Bayesian Information Criterion (BIC): 
#         EII       VII       EEI       VEI       EVI       VVI
# 1 -164785.1 -164785.1 -115645.4 -115645.4 -115645.4 -115645.4
# 2 -162295.2 -161948.7        NA        NA        NA        NA
# 3 -161913.9 -153697.6        NA        NA        NA        NA
# 4 -156442.1 -144833.1        NA        NA        NA        NA
# 5 -156478.7 -143439.6        NA        NA        NA        NA
# 6 -156337.3 -139234.4        NA        NA        NA        NA
# 7 -152566.7 -137717.3        NA        NA        NA        NA
# 8 -152603.3 -137276.6        NA        NA        NA        NA
# 9 -152639.9 -135576.2        NA        NA        NA        NA
#         EEE       VEE       EVV       VVV
# 1 -115267.2 -115267.2 -115267.2 -115267.2
# 2        NA        NA        NA        NA
# 3        NA        NA        NA        NA
# 4        NA        NA        NA        NA
# 5        NA        NA        NA        NA
# 6        NA        NA        NA        NA
# 7        NA        NA        NA        NA
# 8        NA        NA        NA        NA
# 9        NA        NA        NA        NA
# 
# Top 3 models based on the BIC criterion: 
#     EEE,1     EVV,1     VEE,1 
# -115267.2 -115267.2 -115267.2 

BIC(nmml$fit)
# [[1]]
#        EII      VII      EEI       VEI      EVI       VVI EEE
# 1       NA       NA 115645.4 115645.45 115645.4 115645.45  NA
# 2       NA       NA 115762.5  92460.41       NA  90259.66  NA
# 3       NA       NA       NA        NA       NA  79833.80  NA
# 4       NA       NA       NA        NA       NA  62367.79  NA
# 5 156541.7       NA       NA        NA       NA        NA  NA
# 6       NA 139884.9       NA        NA       NA        NA  NA
# 7       NA 139529.2       NA        NA       NA        NA  NA
# 8       NA 138531.6       NA        NA       NA        NA  NA
# 9       NA       NA       NA        NA       NA        NA  NA
#   VEE EVV VVV
# 1  NA  NA  NA
# 2  NA  NA  NA
# 3  NA  NA  NA
# 4  NA  NA  NA
# 5  NA  NA  NA
# 6  NA  NA  NA
# 7  NA  NA  NA
# 8  NA  NA  NA
# 9  NA  NA  NA
# 
# $best
# [1] "4"   "VVI"
# 

BIC(nmml$fit)[[1]] + mcll$BIC
# Bayesian Information Criterion (BIC): 
#        EII       VII EEI VEI EVI          VVI EEE VEE EVV VVV
# 1       NA        NA   0   0   0 9.538839e-06  NA  NA  NA  NA
# 2       NA        NA  NA  NA  NA           NA  NA  NA  NA  NA
# 3       NA        NA  NA  NA  NA           NA  NA  NA  NA  NA
# 4       NA        NA  NA  NA  NA           NA  NA  NA  NA  NA
# 5 63.04807        NA  NA  NA  NA           NA  NA  NA  NA  NA
# 6       NA  650.5279  NA  NA  NA           NA  NA  NA  NA  NA
# 7       NA 1811.8753  NA  NA  NA           NA  NA  NA  NA  NA
# 8       NA 1254.9862  NA  NA  NA           NA  NA  NA  NA  NA
# 9       NA        NA  NA  NA  NA           NA  NA  NA  NA  NA


irt <- iris[,-5]
mclit <- Mclust(irt, G=1:9, modelNames=models)

nmmit <- readRDS(file=file.path(fv_dir, "fit_irt_clr1_clara_nmm.rds"))

mclit$BIC
# Bayesian Information Criterion (BIC): 
#          EII        VII        EEI        VEI        EVI
# 1 -1804.0854 -1804.0854 -1522.1202 -1522.1202 -1522.1202
# 2 -1123.4117 -1012.2352 -1042.9679  -956.2823 -1007.3082
# 3  -878.7650  -853.8144  -813.0504  -779.1566  -797.8342
# 4  -893.6140  -812.6048  -827.4036  -748.4529  -837.5452
# 5  -782.6441  -742.6083  -741.9185  -688.3463  -766.8158
# 6  -715.7136  -705.7811  -693.7908  -676.1697  -774.0673
# 7  -731.8821  -698.5413  -713.1823  -680.7377  -813.5220
# 8  -725.0805  -701.4806  -691.4133  -679.4640  -740.4068
# 9  -694.5205  -700.0276  -696.2607  -702.0143  -767.8044
#          VVI       EEE       VEE       EVV       VVV
# 1 -1522.1202 -829.9782 -829.9782 -829.9782 -829.9782
# 2  -857.5515 -688.0972 -656.3270 -658.3306 -574.0178
# 3  -744.6382 -632.9647 -605.3982 -656.0359 -580.8396
# 4  -751.0198 -646.0258 -604.8371 -725.2925 -630.6000
# 5  -711.4502 -604.8131        NA        NA -676.6061
# 6  -707.2901 -609.8543 -609.5584        NA -754.7938
# 7  -766.6500 -632.4947        NA -809.8276 -806.9277
# 8  -764.1969 -639.2640 -654.8237 -831.7520 -830.6373
# 9  -755.8290 -653.0878        NA -882.4391 -883.6931
# 
# Top 3 models based on the BIC criterion: 
#     VVV,2     VVV,3     EEE,5 
# -574.0178 -580.8396 -604.8131 

BIC(nmmit$fit)
# [[1]]
#         EII       VII       EEI       VEI       EVI       VVI
# 1 1804.0854 1804.0854 1522.1202 1522.1202 1522.1202 1522.1202
# 2 1123.4113 1012.2352 1042.9679  956.2823 1007.3082  857.5515
# 3  878.7639  853.8090  813.0425  779.1502  797.8329  743.9974
# 4  784.2954  783.8168  735.4786  716.5224  732.4552  714.9112
# 5  734.3842  746.9900  694.3888  703.0489  695.6716  700.9021
# 6  715.6928  705.7721  693.7673  675.5589  722.1392  696.8944
# 7  712.0852  712.5276  671.6553  662.5188  704.1570  699.2258
# 8  686.0821  692.0590  661.0738  661.7881  703.6519  709.8228
# 9  694.5147  688.0760  678.5744  671.4087  737.3026  733.3708
#        EEE      VEE      EVV      VVV
# 1 829.9782 829.9782 829.9782 829.9782
# 2 688.0972 656.3270 658.3306 574.0178
# 3 632.9633 605.3968 621.5184 580.8389
# 4 591.4057 611.9207 629.9138 616.8043
# 5 600.5329 602.5561 670.4908 669.9706
# 6 621.8101 618.3594 706.4132 433.2243
# 7 617.5893 617.4666 758.0385 393.4651
# 8 622.4162 626.2729 795.3359 787.2027
# 9 638.2022 640.1674 859.4296 857.2945
# 
# $best
# [1] "7"   "VVV"
# 

BIC(nmmit$fit)[[1]] + mclit$BIC
# Bayesian Information Criterion (BIC): 
#             EII           VII           EEI           VEI
# 1  2.273737e-13  2.273737e-13  6.821210e-13  6.821210e-13
# 2 -4.448742e-04  6.431594e-08 -2.055176e-05  6.435812e-09
# 3 -1.104799e-03 -5.417018e-03 -7.953653e-03 -6.399837e-03
# 4 -1.093186e+02 -2.878799e+01 -9.192505e+01 -3.193045e+01
# 5 -4.825986e+01  4.381677e+00 -4.752970e+01  1.470255e+01
# 6 -2.076121e-02 -8.947392e-03 -2.354640e-02 -6.108759e-01
# 7 -1.979692e+01  1.398625e+01 -4.152700e+01 -1.821894e+01
# 8 -3.899832e+01 -9.421606e+00 -3.033954e+01 -1.767594e+01
# 9 -5.803030e-03 -1.195158e+01 -1.768631e+01 -3.060569e+01
#             EVI           VVI           EEE           VEE
# 1  6.821210e-13  2.989964e-10  0.000000e+00  0.000000e+00
# 2  2.010631e-08  5.820621e-07  3.864592e-07 -6.269241e-07
# 3 -1.257164e-03 -6.407446e-01 -1.399013e-03 -1.404603e-03
# 4 -1.050899e+02 -3.610858e+01 -5.462006e+01  7.083578e+00
# 5 -7.114416e+01 -1.054812e+01 -4.280213e+00            NA
# 6 -5.192812e+01 -1.039569e+01  1.195577e+01  8.800940e+00
# 7 -1.093650e+02 -6.742418e+01 -1.490546e+01            NA
# 8 -3.675494e+01 -5.437412e+01 -1.684782e+01 -2.855085e+01
# 9 -3.050181e+01 -2.245817e+01 -1.488556e+01            NA
#             EVV           VVV
# 1  0.000000e+00  0.000000e+00
# 2  7.579956e-09  1.122521e-08
# 3 -3.451749e+01 -7.226431e-04
# 4 -9.537873e+01 -1.379562e+01
# 5            NA -6.635417e+00
# 6            NA -3.215696e+02
# 7 -5.178908e+01 -4.134627e+02
# 8 -3.641605e+01 -4.343452e+01
# 9 -2.300949e+01 -2.639865e+01
# 
# Top 3 models based on the BIC criterion: 
#    VEI,5    VII,7    EEE,6 
# 14.70255 13.98625 11.95577 


## results all over the place.., no idea what to do




## try mahalanobis distance for matching

mu <- t(MW32$mu)
set.seed(2019);x <- rnorMmix(500, MW32)
nm <- norMmixMLE(x, k=3, model="VEI",trafo="clr1", ini="clara")

ret <- matrix(0, 3, 5)

for (i in 1:5) {
    ret[,i] <- mahalanobis(t(nm$norMmix$mu), center=mu[i,], cov=MW32$Sigma[,,i])
}

ret
#              [,1]      [,2]      [,3]      [,4]      [,5]
# [1,]   0.02723181 14.322284 37.087301 61.969810 87.699315
# [2,] 336.73335067 86.521204 21.117156  1.915131  1.193917
# [3,]  56.86232857  2.751678  2.714795 16.196353 35.085288


####
##------------------------------------------------------------------------------
####
## work on 2019-08-25


## continuation of yesterday
mumw <- t(MW32$mu)
set.seed(2019);x <- rnorMmix(500, MW32)
muf <- norMmixMLE(x, k=3, model="VEI", trafo="clr1", ini="clara")
ret <- matrix(0, 5, 3)
for (i in 1:3) {
    ret[,i] <- mahalanobis(mumw, center=muf$norMmix$mu[,i], cov=muf$norMmix$Sigma[,,i])
}
ret
#              [,1]       [,2]      [,3]
# [1,]   0.02910772 49.8109160 13.581450
# [2,]  31.15541733 25.5992264  1.311877
# [3,] 121.05523549  9.3737762  1.951380
# [4,] 269.72856218  1.1345655 15.499959
# [5,] 477.17539740  0.8815942 41.957614



####
##------------------------------------------------------------------------------
####
## work on 2019-08-27


## reproduce norMmixMLE up untill mstep


data(SMI.12, package="copula")
smi <- SMI.12

options(error = recover)
ans <- norMmixMLE(smi, k=1, model="EII", trafo="clr1", ini="mclVVV", ll="nmm")
## here issue vmmin??

ans <- norMmixMLE(smi, k=7, model="EVI", trafo="clr1", ini="clara", ll="nmm")
## here issue D<0


####
##------------------------------------------------------------------------------
####
## work on 2019-08-30


## fix vmmin issue. put mean in first few cases
ans <- norMmixMLE(smi, k=1, model="EII", trafo="clr1", ini="mclVVV", ll="nmm")

## in theory fixed for first 3 cases
ans <- norMmixMLE(smi, k=2, model="EII", trafo="clr1", ini="mclVVV", ll="nmm")
## no
ans <- norMmixMLE(smi, k=1, model="VII", trafo="clr1", ini="mclVVV", ll="nmm")
## no
ans <- norMmixMLE(smi, k=1, model="EEI", trafo="clr1", ini="mclVVV", ll="nmm")
## converged
ans <- norMmixMLE(smi, k=2, model="VII", trafo="clr1", ini="mclVVV", ll="nmm")
## no


ans <- fit.norMmix(smi, k=1:3, model=1:2, trafo="clr1", ini="clara", ll="nmm")
BIC(ans)
# [[1]]
#   EII VII
# 1  NA  NA
# 2  NA  NA
# 3  NA  NA
# 
# $best
# character(0)
# 


## try again with force positive
ans <- norMmixMLE(smi, k=1, model="EII", trafo="clr1", ini="mclVVV", ll="nmm")


ans <- fit.norMmix(smi, k=1:3, model=1:2, trafo="clr1", ini="clara", ll="nmm")
BIC(ans)
# [[1]]
#   EII VII
# 1  NA  NA
# 2  NA  NA
# 3  NA  NA
# 
# $best
# character(0)
# 

## no
ans <- norMmixMLE(smi, k=7, model="VVI", trafo="clr1", ini="mclVVV", ll="nmm")
## something wrong here too


## error handling
## fit.norMmix now saves proper error

## list with dim

x9 <- rnorMmix(100, MW29)
ans <- fit.norMmix(x9, k=1:3, model=1:2, trafo="clr1", ini="clara", ll="nmm")
rr <- ans$nMm
dim(rr) <- c(3,2)

## can't acces sublists after assigning dims


####
##------------------------------------------------------------------------------
####
## work on 2019-08-31


## try again with numerical stability
## make epsilon accessible variable
ans <- norMmixMLE(smi, k=1, model="VII", trafo="clr1", ini="mclVVV", ll="nmm", epsilon=1e-1)
## first converges with epsilon at 1e-1

ans <- fit.norMmix(smi, k=1:3, model=1:2, trafo="clr1", ini="clara", ll="nmm", epsilon=1e-1)
BIC(ans)
# [[1]]
#        EII      VII
# 1 27040.40 27040.40
# 2 23616.70 23210.41
# 3 21952.26 21779.22
# 
# $best
# [1] "VII"
# 

## this fixes it
