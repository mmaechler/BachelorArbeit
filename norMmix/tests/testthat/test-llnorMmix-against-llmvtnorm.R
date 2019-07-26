context("test llnorMmix against llmvtnorm to check values of llnorMmix")





test_that("EII llnorMmix against llmvtnorm", {
		  tr <- "clr1"
		  mo <- "EII"

		  par. <- par2nm1
		  x <- rnorMmix(obj=MW2nm1)
		  p <- 2
		  k <- 1

		  retnMm <- llnorMmix(par., x, p, k, trafo=tr, model=mo)
		  retmvt <- llmvtnorm(par., x, p, k, trafo=tr, model=mo)

		  expect_equal(retnMm,retmvt)
})
			


test_that("VII llnorMmix against llmvtnorm", {
		  tr <- "clr1"
		  mo <- "VII"

		  set.seed(2019)
		  par. <- par2nm2
		  x <- rnorMmix(obj=MW2nm2)
		  p <- 2
		  k <- 3

		  retnMm <- llnorMmix(par., x, p, k, trafo=tr, model=mo)
		  retmvt <- llmvtnorm(par., x, p, k, trafo=tr, model=mo)

		  expect_equal(retnMm,retmvt)
})
			



test_that("EEI llnorMmix against llmvtnorm", {
		  tr <- "clr1"
		  mo <- "EEI"

		  par. <- par26
		  x <- rnorMmix(obj=MW26)
		  p <- 2
		  k <- 2

		  retnMm <- llnorMmix(par., x, p, k, trafo=tr, model=mo)
		  retmvt <- llmvtnorm(par., x, p, k, trafo=tr, model=mo)

		  expect_equal(retnMm,retmvt)
})


test_that("VEI llnorMmix against llmvtnorm", {
		  tr <- "clr1"
		  mo <- "VEI"

		  par. <- par27
		  x <- rnorMmix(obj=MW27)
		  p <- 2
		  k <- 2

		  retnMm <- llnorMmix(par., x, p, k, trafo=tr, model=mo)
		  retmvt <- llmvtnorm(par., x, p, k, trafo=tr, model=mo)

		  expect_equal(retnMm,retmvt)
})


test_that("EVI llnorMmix against llmvtnorm", {
		  tr <- "clr1"
		  mo <- "EVI"

		  par. <- par28
		  x <- rnorMmix(obj=MW28)
		  p <- 2
		  k <- 3

		  retnMm <- llnorMmix(par., x, p, k, trafo=tr, model=mo)
		  retmvt <- llmvtnorm(par., x, p, k, trafo=tr, model=mo)

		  expect_equal(retnMm,retmvt)
})


test_that("VVI llnorMmix against llmvtnorm", {
		  tr <- "clr1"
		  mo <- "VVI"

		  par. <- par29
		  x <- rnorMmix(obj=MW29)
		  p <- 2
		  k <- 2

		  retnMm <- llnorMmix(par., x, p, k, trafo=tr, model=mo)
		  retmvt <- llmvtnorm(par., x, p, k, trafo=tr, model=mo)

		  expect_equal(retnMm,retmvt)
})


test_that("EEE llnorMmix against llmvtnorm", {
		  tr <- "clr1"
		  mo <- "EEE"

		  par. <- par210
		  x <- rnorMmix(obj=MW210)
		  p <- 2
		  k <- 2

		  retnMm <- llnorMmix(par., x, p, k, trafo=tr, model=mo)
		  retmvt <- llmvtnorm(par., x, p, k, trafo=tr, model=mo)

		  expect_equal(retnMm,retmvt)
})


test_that("VEE llnorMmix against llmvtnorm", {
		  tr <- "clr1"
		  mo <- "VEE"

		  par. <- par211
		  x <- rnorMmix(obj=MW211)
		  p <- 2
		  k <- 2

		  retnMm <- llnorMmix(par., x, p, k, trafo=tr, model=mo)
		  retmvt <- llmvtnorm(par., x, p, k, trafo=tr, model=mo)

		  expect_equal(retnMm,retmvt)
})


test_that("EVV llnorMmix against llmvtnorm", {
		  tr <- "clr1"
		  mo <- "EVV"

		  par. <- par212
		  x <- rnorMmix(obj=MW212)
		  p <- 2
		  k <- 2

		  retnMm <- llnorMmix(par., x, p, k, trafo=tr, model=mo)
		  retmvt <- llmvtnorm(par., x, p, k, trafo=tr, model=mo)

		  expect_equal(retnMm,retmvt)
})


test_that("VVV llnorMmix against llmvtnorm", {
		  tr <- "clr1"
		  mo <- "VVV"

		  par. <- par213
		  x <- rnorMmix(obj=MW213)
		  p <- 2
		  k <- 2

		  retnMm <- llnorMmix(par., x, p, k, trafo=tr, model=mo)
		  retmvt <- llmvtnorm(par., x, p, k, trafo=tr, model=mo)

		  expect_equal(retnMm,retmvt)
})
