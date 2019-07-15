context("test input/output of par2nMm")



test_that("io test using MWnm", {
		  
		  n2p <- function(obj) { nMm2par(obj, trafo="clr1", model=obj$model) }
		  tr <- "clr1"

		  m <- MW2nm1
		  a <- n2p(m)
		  b <- par2nMm(a,m$dim,m$k,trafo=tr,m$model)

		  expect_equal(m$weight,b$weight)
		  expect_equal(m$mu,b$mu)
		  expect_equal(m$Sigma[,,1],b$Sigma[,,1])
		  expect_equal(m$k,b$k)
		  expect_equal(m$dim,b$dim)

		  m <- MW2nm2
		  a <- n2p(m)
		  b <- par2nMm(a,m$dim,m$k,trafo=tr,m$model)

		  expect_equal(m$weight,b$weight)
		  expect_equal(m$mu,b$mu)
		  expect_equal(m$Sigma[,,1],b$Sigma[,,1])
		  expect_equal(m$k,b$k)
		  expect_equal(m$dim,b$dim)







})
