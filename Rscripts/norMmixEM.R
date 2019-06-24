### EM-algorithm for norMmix datasets

## have not figured out how to do the details here,
## for now only scetch of programm !!not functional!!

library(mixtools)


# M step

mstep.eMm <- function( y, par ){
	
	if (!is.matrix(y)) stop("in mstep y is not matrix")

	#  calculating T_i

	T1[,i] <- rowSums(tau[i,])

	T2[,i] <- rowSums(tau[i,]*y[??])  # sum_j tau_ij*y_j

	T3[,i] <- rowSums(tau[i,]* y[??]%*%y[??])


	#mu

	mu[i,] <- T2[,i]/T1[,i]

	#sigma

	Sigma[,,i] <- ( T3[,i]-T1[,i]^-1 *(T2[,i]%*%T2[,i]) ) / T1[,i]

	# weight

	weight <- T1[,i]/nrow(y)

	##return params

}


estep.nMm <- function(y, par){

	phi <- matrix(0,k,n)
	
	for (i in 1:k){
	phi[i,] <- dmvnorm(y, mu=mu[i,], sigma=Sigma[,,i]) # a k x n vec
	}

	T5 <- rep(0,n) 
	T5 <- colSums( phi %*% diag(weight))^-1

	tau <- (weight %*% phi) %*% T5
	

	## return params
}

	
