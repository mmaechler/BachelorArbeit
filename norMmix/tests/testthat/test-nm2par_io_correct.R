context("test io correctness of nMm2par")


test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("test length", {
		  tra <- "clr1"
		  t1 <- nMm2par(MW2nm1, trafo=tra, model="EII")
		  expect_equal(length(t1),3)
		  t2 <- nMm2par(MW2nm2, trafo=tra, model=MW2nm2$model)
		  expect_equal(length(t2),11)
		  t4 <- nMm2par(MW2nm4, trafo=tra, model=MW2nm4$model)
		  expect_equal(length(t4),7)
		  t6 <- nMm2par(MW26, trafo=tra, model=MW26$model)
		  expect_equal(length(t6),8)
		  t7 <- nMm2par(MW27, trafo=tra, model=MW27$model)
		  expect_equal(length(t7),9)
})
