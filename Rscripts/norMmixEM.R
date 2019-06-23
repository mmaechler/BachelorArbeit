### EM-algorithm for norMmix datasets

## have not figured out how to do the details here,
## for now only scetch of programm !!not functional!!


# M step

mstep.eMm <- function( x, par ){
	
	if (!is.matrix(x)) stop("in mstep x is not matrix")

	#  calculating T_i

	T1[,i] <- rowSums(tau[i,])

	T2[,i] <- rowSums(tau[i,]*y[??])  # sum_j tau_ij*y_j

	T3[,i] <- rowSums(tau[i,]* y[??]%*%y[??])


	#mu

	mu[i,] <- T2[,i]/T1[,i]

	#sigma

	Sigma[,,i] <- ( T3[,i]-T1[,i]^-1 *(T2[,i]%*%T2[,i]) ) / T1[,i]

	# weight

	weight <- T1[,i]/nrow(x)

	##return params
}


estep.nMm <- function(x, par){

	phi[i,j] <- #probability of y_j with param mu_i Sigma_i

	T5[,i] <- rowSums( weight[,i]*phi[j,i] )

	tau[,i] <- weight[,i]*phi[i,j]/T5[,j]
	

	## return params
}

	
