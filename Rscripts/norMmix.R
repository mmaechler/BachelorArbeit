#### the extra m stands for multivariate





## norMmix constructor

tole <- 1000*.Machine$double.eps

library(abind)
library(matrixcalc)
library(MASS)
library(scatterplot3d)

Shilf <- function(L, p){ # help function for norMmix
	Shilfarray <- array(0, c(p,p,length(L)))
	for (i in 1:length(L)){
		Shilfarray[,,i] <- L[i]*diag(p)
	}
	return(Shilfarray)
}








hilf <- function(x, k){ # help function evals to true if x is sym pos def
	for (i in 1:k){
		if ( !(isSymmetric(x[,,i], tol=tole)&& is.positive.semi.definite(x[,,i],tol=tole)) )
			stop("Sigma is not sym pos def")
	}
	return(TRUE)
}







norMmix <- function(
		    mu,
		    sigma = abind( rep(list(diag(p)), k), along = 3),
		    weight = rep(1/k, k),
		    name = NULL
#		    type = c("")
		    )
{
	## Purpose: constructor for 'norMmix' (multivariate normix)
	## makes sure values are sane.
	## --------------------------------------------------------
	## Arguments:
	##	mu: matrix with columns as vector means of dist.
	##	sigma: 	option 0: default, generates all identity
	##			covariance mat.
	##		option 1: scalar, generates EII dist
	##		option 2: vector of length k, generates VII
	##			distribution
	##		option 3: array of dimension p x p x k. 
	##			covariance matrices of distributions
	##	weight: vector of length k, sums to 1
	##	name: name attribute
	##	type: type of distribution VVV, IVV etc.
	## --------------------------------------------------------
	## Value: returns objext of class 'norMmix'
	## --------------------------------------------------------
	## Author: nicolas trutmann, Date:2019-06-19

	# p dimension 
	# k number of components

	# ispect mu 
	if (!is.numeric(mu)) stop("'mu' must be numeric")
	if (is.matrix(mu) == FALSE){
		p <- 1
		k <- length(mu)
	} else {
		k <- ncol(mu)
		p <- nrow(mu)
	}


	#ispect sigma
	if (!is.numeric(sigma)) stop("'sigma' must be numeric")
	if (is.vector(sigma) && length(sigma)==1)
		Sigma <- abind( rep(list(sigma*diag(p)), k), along = 3)
	else if (is.vector(sigma) && length(sigma)==k)
		Sigma <- Shilf(L=sigma, p=p) 
	else if (is.array(sigma) && dim(sigma) == c(p,p,k))
		Sigma <- sigma
	else stop("'sigma not among recognized formats, see source code for help'")
	if (!hilf(Sigma,k)) stop("error with sym pos def Sigma")

	#inspect weight
	if (!is.numeric(weight)) stop("'weight' must be numeric")
	if (! (is.vector(weight) && length(weight)==k) )
		stop("weight is not of correct dimension")
	if (!(all(weight >=0) && (sum(weight) - 1 < 1000*.Machine$double.eps)))
		stop("weight doesn't sum to 1 or isn't positive")

	structure( name = name,
		  class = "norMmix",
		  .Data = list(mu = mu , Sigma = Sigma, weight = weight,
		  		k=k, dim=p)
	)

}






is.norMmix <- function(obj){
	inherits(obj, "norMmix")
}






mean.norMmix <- function(obj){
	if (!is.norMmix(obj)) stop("object is not norMmix")
	k <- obj$k
	mu <- obj$mu
	w <- obj$weight

	me <- rep(0, obj$dim)

	for (i in 1:k) {
		me <- me + w[i]*mu[,i]
	}
	return(me)
}






### rnorMmix

rnorMmix <- function(
		     obj,
		     n = 511
	  )
{

	## Purpose: generates random values distributed by NMD
	## -------------------------------------------------------------------
	## Arguments:
	## 	obj: of type norMmix
	##	n: number of values desired
       	## -------------------------------------------------------------------
	## Value: matrix, columns are vectors
	## -------------------------------------------------------------------
	## Author: nicolas trutmann, Date:2019-06-21

	if(!inherits(obj, "norMmix")) stop("'obj' must be of type norMmix")

	mu <- obj$mu
	Sigma <- obj$Sigma
	weight <- obj$weight
	p <- obj$dim

	#need case n=1 distiction here

	nj <- rmultinom(n=1, size=n, prob=weight)
#	a <- matrix(unlist(lapply( seq(along=nj), function(j)mvrnorm(n=nj[j], mu=mu[,j], Sigma=Sigma[,,j]) )), ncol=p, byrow=TRUE)
	## this approach doesnt work matrices arent concatenated properly
	a <- do.call( rbind, lapply( seq(along=nj), function(j)mvrnorm(n=nj[j], mu=mu[,j], Sigma=Sigma[,,j]) ))
}
