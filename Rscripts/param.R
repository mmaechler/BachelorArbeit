#### functions handling parameter manipulation. par2nM nM2par etc


source(file="Cholesky.R")

## a way to make upper/lower triangular matrix into vector without lower/upper half

# suppose l is a matrix:
#l.index <- upper.tri(l, diag=TRUE)
#x <- l[l.index==TRUE]

#reverse
#l[upper.tri(l, diag=TRUE)==TRUE] <- x

## works


#wrapper function for lower.tri to vector

ld. <- function(mat){
	x <- mat[lower.tri(mat, diag=FALSE)]
	x
}

dl. <- function(d,x,p){
	mat <- matrix(0,p,p)
	mat[lower.tri(mat,diag=FALSE)==TRUE] <- x
	mat.temp <- mat + diag(rep(1,p))  
	mat <- mat.temp %*% diag(d) %*% t(mat.temp)
	mat
}



nMm2par <- function(obj, 
		    trafo=c("clr1", "logit"),
		    model=c("EII","VII","EEI","VEI","EVI","VVI","EEE",
		    "EVE","VEE","VVE","EEV","VEV","EVV","VVV")
		    ){

	#transferring values of obj to handier variables
	sig <- obj$Sigma
	mu <- obj$mu
	w <- obj$w
	k <- obj$k
	p <- obj$dim

	trafo <- match.arg(trafo)
	model <- match.arg(model)


	#output vector of parameter values
	c(
	  k, #number of clusters
	  p, #dimension
	  model, #model name
	  w <- switch(trafo, #weights either logit or centered log ratio
		    "logit" = {
			    w <- w/sum(w) #necessary?
			    wlog <- log(w/(1+w))[-1L]
		    },

		    "clr1" = {
			    lw <- log(w)
			    lp <- lw - mean(lw)
			    lp
		    	     },

		    stop("invalid argument trafo, ",trafo)
		    ),
	  mu, #means
	  Sigma <- switch(model, #model dependent covariance values
			  "EII" = sig[1,1,1],

			  "VII" = sig[1,1,],

			  "EEI" = {D. <- diag(sig[,,1])
			  	   alpha <- (prod(D.)^(1/p))
			  	   c(alpha, D)
			  	  },

			  "VEI" = {alpha <- apply(sig, 3, function(j) prod(diag(j))^(1/p) )
			  	   D. <- diag(sig[,,1])/alpha[1]
			  	   c(alpha, D.)
				  },
			  
			  "EVI" = {alpha <- prod(diag(sig[,,1]))^(1/p)
			  	   D. <- apply(sig,3,diag)/alpha
			  	   c(alpha,D.)
			  	  },

			  "VVI" = {alpha <- apply(sig,3, function(j) det(j)^(1/p)) 
			  	   D.temp <- apply(sig,3,diag)
			 	   D. <- D.temp %*% diag(1/alpha) # this is fastest?? https://stackoverflow.com/questions/20596433/how-to-divide-each-row-of-a-matrix-by-elements-of-a-vector-in-r
				  },

			  "EEE" = {A.temp <- ldl(sig)
			  	   L <- A.temp$L
			  	   D.temp <- A.temp$Diag
				   D. <- D.temp/(lambda <- prod(D.temp)^(1/p))
				   c( alpha, D., ld.(L) )
				  },

			  "EVE" = {A.temp <- ldl(sig[,,1])
			  	   D.temp <- A.temp$Diag

			  	   L <- A.temp$L
				   alpha <- prod(D.temp)^(1/p)
				   D. <- c( apply(sig,3, function(j) ldl(j)$Diag/alpha) )

				   c( alpha, D., ld.(L) )
			  	  },

			  "VEE" = {alpha <- apply(sig,3, function(j) prod(ldl(j)$Diag)^(1/p))
			  	   L <- ldl(sig[,,1])$L
			  	   D. <- ldl(sig[,,1])/alpha[1]

				   c(alpha, D., ld.(L))
				  },

			  "VVE" = {alpha <- apply(sig,3, function(j) prod(ldl(j)$Diag)^(1/p))
			  	   L <- ldl(sig[,,1]$L)
			  	   D.temp <- apply(sig,3, function(j) ldl(j)$Diag)
				   D. <- D.temp %*% diag(1/alpha)

				   c(alpha, D., ld.(L))
				  },

			  "EEV" = ,

			  "VEV" = ,

			  #"EVV" = ,same as VVV

			  #"VVV" = {S <- matrix(p*(p+1L),k)
#			  for (i in 1:k){
#					S[,i] <- ldl(sig[,,i])$L
#				}
#			  },

			  stop("invalid argument in 'model'")
			  )
	)
}




par2nMm <- function(p){

	##
	##

	k <- p[1L]
	p <- p[2L]
	model <- p[3L]

	w <- switch(trafo,
		    "plogit" = ,

		    "clr1" = 
		    )

	mu <- matrix(0, c(p,k))

	for (i in 1:k){
		mu[,i] <- p[(3L+i*p):(2L+(i+1L)*p)]
	}


	Sigma <- switch()


}



