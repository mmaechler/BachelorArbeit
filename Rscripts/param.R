#### functions handling parameter manipulation. par2nM nM2par etc



## a way to make upper/lower triangular matrix into vector without lower/upper half

# suppose l is a matrix:
#l.index <- upper.tri(l, diag=TRUE)
#x <- l[l.index==TRUE]

#reverse
#l[upper.tri(l, diag=TRUE)==TRUE] <- x

## works

udtv <- function(mat){
	upper.tri(chol(mat), diag=TRUE)
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

			  "EEI" = diag(sig[,,1]),

			  "VEI" = c( sig[1,1,]/min(sig[1,1,]), sig[,,1]*min(sig[1,1,])/sig[1,1,1] ),
			  
			  "EVI" = ,

			  "VVI" = ,

			  "EEE" = udtv(sig[,,1]),

			  "EVE" = ,

			  "VEE" = ,

			  "VVE" = ,

			  "EEV" = ,

			  "VEV" = ,

			  #"EVV" = ,same as VVV

			  "VVV" = {S <- matrix(p*(p+1L),k)
			  	for (i in 1:k){
					S[,i] <- udtv(sig[,,i])
				}
			  }

			  stop("invalid argument in 'model'")
			  )
	)
}




par2nMm <- function(p,  ??){

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



