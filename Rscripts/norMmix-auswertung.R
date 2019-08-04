## 


devtools::load_all("~/ethz/BA/norMmix")
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
