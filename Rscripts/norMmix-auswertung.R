## 

devtools::load_all("~/ethz/BA/norMmix")
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
