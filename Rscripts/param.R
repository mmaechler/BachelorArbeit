#### functions handling parameter manipulation. par2nM nM2par etc

### for info on models see Celeux and Govaert 1995


source(file="Cholesky.R")
source(file="weight.R")

# temp source file
source(file="norMmix.R")
source(file="zmarrwandnMm.R")

## map lower.tri to vec
ld. <- function(mat){
	x <- mat[lower.tri(mat, diag=FALSE)]
	x
}

## map vec to lower.tri
dl. <- function(d,x,p){
	mat <- matrix(0,p,p)
	mat[lower.tri(mat,diag=FALSE)==TRUE] <- x
	mat.temp <- mat + diag(rep(1,p))  
	mat <- mat.temp %*% diag(d) %*% t(mat.temp)
	mat
}



nMm2par <- function(obj, 
		    trafo=c("clr1", "logit"),
		    model=c("EII","VII","EEI","VEI",
			    "EVI","VVI","EVV","VVV")
		    ){
	#-------------------------------------------------
	#
	#
	#
	#
	#-------------------------------------------------

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
	  w <- switch(trafo, #weights either logit or centered log ratio
		    "logit" = qlogis(w[-1L]),

		    "clr1" = clr1(w),

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

			  "EVV" = {alpha <- prod( ldl(sig[,,1])$Diag )^(1/p)
			  	   D. <- apply(sig,3, function(j) ldl(j)$Diag/alpha)
			  	   L. <- apply(sig,3, function(j) ld.( ldl(j)$L ))
				   c(alpha, D., L.)},

			  "VVV" = {alpha <- apply(sig,3, function(j) prod(ldl(j)$Diag )^(1/p))
			  	   D. <- apply(sig,3, function(j) ldl(j)$Diag/alpha)
			  	   L. <- apply(sig,3, function(j) ld.( ldl(j)$L ))
			  	   c(alpha, D., L.)},

			  stop("invalid argument in 'model'")
			  )
	)
}




par2nMm <- function(par., p, k,
		    trafo=c("clr1", "logit"),
		    model=c("EII","VII","EEI","VEI",
			    "EVI","VVI","EVV","VVV"),
		    MLE=FALSE) {

	trafo <- match.arg(trafo)
	model <- match.arg(model)

	f <- k+p*k # start of Sigma/alpha
	f1 <- f + k - 1L # end of alpha if variable

	w.temp <- par.[1:(k-1)]
	w <- switch(trafo,
	"plogit" = {w.temp <- plogis(w.temp)
		    if ((sp. <- sum(w.temp)) > 1) {
			    stop("weights add to greater than 1")
		    }
		    w <- c((1-sp.), w.temp)
		   },

	"clr1" = clr1inv(w.temp), 
		    )

	mu <- matrix(par[k:(k+p*k-1)], p, k)



	Sigma <- switch(model,
	
	# diagonal cases
	
	"EII" = {lambda <- par.[f]

		 sig <- array( rep(diag(lambda, p),k), c(p,p,k) )},

	"VII" = {lambda <- par.[f:f1]

		 sig <- array(unlist(lapply( lambda, function(j) diag(j,p) )), c(p,p,k)) },

	"EEI" = {lambda <- par.[f]

		 D.temp <- par.[f+1:f+p]
		 
		 sig <- array( rep(diag(lambda*D.temp),k), c(p,p,k) )},

	"VEI" = {lambda <- par.[f:f1]
			
		 D.temp <- par.[f1+1:f1+p]
		
		 D. <- tcrossprod(D.temp,lambda)

		 sig <- array( apply(D.,2, diag), c(p,p,k)) },

	"EVI" = {lambda <- par.[f]
			
		 D.temp <- matrix(par.[f1+1:f1+p*k],p,k)

		 D. <- D.temp*lambda

		 sig <- array( apply(D.,2, diag), c(p,p,k)) },

	"VVI" = {lambda <- par.[f:f1]
			
		 D.temp <- matrix(par.[f1+1:f1+p*k],p,k)
	
		 D. <- D.temp %*% diag(lambda)

		 sig <- array( apply(D.,2, diag), c(p,p,k)) },

	# variable cases

	"EVV" = {lambda <- par.[f]

		 D.temp <- matrix(par.[f1+1:f1+p*k],p,k)

		 f2 <- f1 + p*k +1 # start of L_k
		 f3 <- (p*(p-1)/2)
		 f4 <- (p*(p-1)/2)*k
		 f5 <- f2 + f4 -1

		 L.temp <- matrix(par.[f2:f5],f3,k)

		 sig <- array(0, c(p,p,k))
		 for (i in 1:k) {
			 sig[,,i] <- dl.(D.temp[,i],L.temp[,i],p)
		 }
		},
		 

	"VVV" = {lambda <- par.[f:f1]

		 D.temp <- matrix(par.[f1+1:f1+p*k],p,k)

		 D. <- D.temp %*% diag(lambda)

		 f2 <- f1 + p*k +1 # start of L_k
		 f3 <- (p*(p-1)/2)
		 f4 <- (p*(p-1)/2)*k
		 f5 <- f2 + f4 -1

		 L.temp <- matrix(par.[f2:f5],f3,k)

		 sig <- array(0, c(p,p,k))
		 for (i in 1:k) {
			 sig[,,i] <- dl.(D.temp[,i],L.temp[,i],p)
		 }
		},
	stop("error in Sigma switch statement")
	
	)


	list( w, mu, Sigma )

}



