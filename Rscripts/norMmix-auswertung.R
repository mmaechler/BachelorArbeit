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

