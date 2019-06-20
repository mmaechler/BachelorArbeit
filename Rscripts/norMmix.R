### the extra m stands for multivariate


## norMmix constructor

norMmix <- function(
		    mu = ,
		    Sigma = ,
		    weight = NULL,
		    name = NULL,
		    type = c(""), 
		    )
{
	## Purpose: constructor for 'norMmix' (multivariate normix)
	## makes sure values are sane.
	## --------------------------------------------------------
	## Arguments:
	##	mu: matrix with columns as vector means of dist.
	##	Sigma: array of dimension p x p x k. covariance mat
	##		rices of distributions
	##	weight: vector of length k, sums to 1
	##	name: name attribute
	##	type: type of distribution VVV, IVV etc.
	## --------------------------------------------------------
	## Value: returns objext of class 'norMmix'
	## --------------------------------------------------------
	## Author: nicolas trutmann, Date:2019-06-19

	# p dimension 
	# k number of components

	if(!is.numeric(mu)) stop("'mu' must be numeric")
	k <- ncol(mu)
	p <- nrow(mu)



	structure( name = name, class = norMmix, .Data = list(mu, Sigma, weight) )

}
