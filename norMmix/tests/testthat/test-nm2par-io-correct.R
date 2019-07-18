context("test io correctness of nMm2par")


test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("test length", {
		  tra <- "clr1"
		  t1 <- nMm2par(MW2nm1, trafo=tra, model=MW2nm1$model)
		  expect_equal(length(t1),3)
		  t2 <- nMm2par(MW2nm2, trafo=tra, model=MW2nm2$model)
		  expect_equal(length(t2),11)
		  t4 <- nMm2par(MW2nm4, trafo=tra, model=MW2nm4$model)
		  expect_equal(length(t4),7)
		  t6 <- nMm2par(MW26, trafo=tra, model=MW26$model)
		  expect_equal(length(t6),8)
		  t7 <- nMm2par(MW27, trafo=tra, model=MW27$model)
		  expect_equal(length(t7),9)
		  t8 <- nMm2par(MW28, trafo=tra, model=MW28$model)
		  expect_equal(length(t8),15)
		  t9 <- nMm2par(MW29, trafo=tra, model=MW29$model)
		  expect_equal(length(t9),11)
		  t10 <- nMm2par(MW210, trafo=tra, model=MW210$model)
		  expect_equal(length(t10),9)
		  t11 <- nMm2par(MW211, trafo=tra, model=MW211$model)
		  expect_equal(length(t11),10)
		  t12 <- nMm2par(MW212, trafo=tra, model=MW212$model)
		  expect_equal(length(t12),12)
		  t13 <- nMm2par(MW213, trafo=tra, model=MW213$model)
		  expect_equal(length(t13),13)
})




test_that("test if n2p(p2n()) == id ", {

		  tr <- "clr1"
		  m <- MW211
		  l <- par211
		  k <- nMm2par(par2nMm(l,m$dim,m$k,trafo=tr,m$model),
			       trafo=tr,model=m$model)
		  expect_equal(k,l)

		  m <- MW212
		  l <- par212
		  k <- nMm2par(par2nMm(l,m$dim,m$k,trafo=tr,m$model),
			       trafo=tr,model=m$model)
		  expect_equal(k,l)

		  m <- MW213
		  l <- par213
		  k <- nMm2par(par2nMm(l,m$dim,m$k,trafo=tr,m$model),
			       trafo=tr,model=m$model)
		  expect_equal(k,l)

		  m <- MW26
		  l <- par26
		  k <- nMm2par(par2nMm(l,m$dim,m$k,trafo=tr,m$model),
			       trafo=tr,model=m$model)
		  expect_equal(k,l)

		  m <- MW27
		  l <- par27
		  k <- nMm2par(par2nMm(l,m$dim,m$k,trafo=tr,m$model),
			       trafo=tr,model=m$model)
		  expect_equal(k,l)

		  m <- MW28
		  l <- par28
		  k <- nMm2par(par2nMm(l,m$dim,m$k,trafo=tr,m$model),
			       trafo=tr,model=m$model)
		  expect_equal(k,l)

		  m <- MW29
		  l <- par29
		  k <- nMm2par(par2nMm(l,m$dim,m$k,trafo=tr,m$model),
			       trafo=tr,model=m$model)
		  expect_equal(k,l)

		  m <- MW2nm1
		  l <- par2nm1
		  k <- nMm2par(par2nMm(l,m$dim,m$k,trafo=tr,m$model),
			       trafo=tr,model=m$model)
		  expect_equal(k,l)

		  m <- MW2nm1.2
		  l <- par2nm1.2
		  k <- nMm2par(par2nMm(l,m$dim,m$k,trafo=tr,m$model),
			       trafo=tr,model=m$model)
		  expect_equal(k,l)

		  m <- MW2nm2
		  l <- par2nm2
		  k <- nMm2par(par2nMm(l,m$dim,m$k,trafo=tr,m$model),
			       trafo=tr,model=m$model)
		  expect_equal(k,l)

		  m <- MW2nm4
		  l <- par2nm4
		  k <- nMm2par(par2nMm(l,m$dim,m$k,trafo=tr,m$model),
			       trafo=tr,model=m$model)
		  expect_equal(k,l)

		  m <- MW2nm5
		  l <- par2nm5
		  k <- nMm2par(par2nMm(l,m$dim,m$k,trafo=tr,m$model),
			       trafo=tr,model=m$model)
		  expect_equal(k,l)

})




test_that("test wrong inputs in mu", {
		  
		  tr <- "clr1"
		  o1 <- MW210
		  o1$mu <- cbind( c(0,0),c(0,0),c(0,0) )
		  expect_error(nMm2par(o1, trafo=tr,model=o1$model))

		  o1$mu <- cbind( c(0,0,0),c(0,0,0) )
		  expect_error(nMm2par(o1, trafo=tr,model=o1$model))

		  o1$mu <- cbind( c(0,0),c(0,"hello") )
		  expect_error(nMm2par(o1, trafo=tr,model=o1$model))

		  o1$mu <- cbind( c(0,0),c(0,NA) )
		  expect_error(nMm2par(o1, trafo=tr,model=o1$model))
})


test_that("test wrong inputs in w", {
		  tr <- "clr1"
		  o2 <- MW28

		  o2$weight <- c(1,1)
		  expect_error(nMm2par(o2, trafo=tr,model=o1$model))
		  
		  o2$weight <- c(1,1)/2
		  expect_error(nMm2par(o2, trafo=tr,model=o1$model))
		  
		  o2$weight <- c(1,1,1,1)/4
		  expect_error(nMm2par(o2, trafo=tr,model=o1$model))



})




