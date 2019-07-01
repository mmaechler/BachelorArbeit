### testing accuracy of clr trafo for different logs


#### conclusion: should use natural log, smallest error.


clr1 <- function(w, base=exp(1)){

	stopifnot(is.numeric(w) && all(w >= 0))
	# calculate clr1
	ln <- log(w, base=base)
	res <- ln[-1L] - mean(ln)
	res
}


clr1inv <- function(p, base=exp(1)){

	stopifnot(is.numeric(p))
	# calc weights
	p1 <- exp(c(-sum(p), p))
	f <- sum(p1)
	pp <- p1/f
	pp
}

clrtest <- function(w, base=exp(1)){
	w - clr1inv(clr1(w,base=base), base=base)
}


test <- 1:7000
test <- test/sum(test)
test

diffe <- clrtest(test)
difff <- clrtest(test,base=10)
diffg <- clrtest(test,base=3)

plot(diffe) # interresting error pattern
plot(difff) # a lot worse than exp(1)
plot(diffg) # same



## smaller sample -> bigger values


test <- 1:50
test <- test/sum(test)

diffh <- clrtest(test)

plot(diffh)


test <- c(1,50,999)
test <- test/sum(test)

diffi <- clrtest(test)

plot(diffi)

## error seems stable in +- 1e-17, smaller than machine precision




