context("test results of llnorMmix function for calculable cases")


test_that("EII test", {

		  tr <- "clr1"
		  m <- "EII"

		  par. <- c(0,0,0)
		  x <- rbind(c(0,0), c(1,1))
		  p <- 2
		  k <- 1
		  
		  ret <- llnorMmix(par., x, p, k, trafo=tr, model=m)
		  expect_equal(ret, (log(1/(2*pi))+log(1/(2*pi*exp(1)))))

		  par. <- c(0,0,log(2))
		  x <- rbind(c(0,0), c(1,1))
		  p <- 2
		  k <- 1
		  
		  ret <- llnorMmix(par., x, p, k, trafo=tr, model=m)
		  expect_equal(ret, 2*(-log(2*pi)- 0.5*log(4)) -2)



})


test_that("VII test", {

		  tr <- "clr1"
		  m <- "VII"

		  par. <- c(0,0,0,0,log(2))
		  x <- rbind(c(0,0), c(1,1))
		  p <- 2
		  k <- 2
		  
		  ret <- llnorMmix(par., x, p, k, trafo=tr, model=m)
#		  expect_equal(ret,



})
