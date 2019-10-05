## script tests a) behaviour under smple size change b) behaviour under seed change

testsize <- FALSE
testseed <- TRUE

GH_BA_dir <- normalizePath("~/BachelorArbeit")
save_dir  <- normalizePath("~/4NT/smallseed")
##-------------------------------------------------
stopifnot(dir.exists(GH_BA_dir),
          dir.exists(save_dir))

## MM: Dies ist *NICHT* die Art ein Paket zu laden .. 
##     aber "good enough for now":
devtools::load_all(file.path(GH_BA_dir, "norMmix"))
size <- 5:15*100
seedx <- 1:10
seedmle <- 1:100

if (testsize) {
    for (seed in seedx) {
        for (s in size) {

            nm <- MW214
            set.seed(2019+seed); x <- rnorMmix(s,nm)
            st <- system.time(
                r <- tryCatch(fitnMm(x, k=1:8, models=1:10,
                                          trafo="clr1", maxit=1e4),
                              error = identity)
                )

            bics <- BIC(r)

            sFilepdf <- sprintf("bic_fit_MW214_size=%d_seed=%d.pdf",
                                s, seed)
            pdf(file=file.path(save_dir,sFilepdf))
            plot(r, name=sFilepdf)
            dev.off()

            cat("result of fit: "); str(r, max=1)
            sFile <- sprintf("fit_MW214_size=%d_seed=%d.rds",
                             s, seed)
            cat("--> saving to file:", sFile, "\n")
            saveRDS(list(st = st, fit = r, bics),
                    file=file.path(save_dir, sFile))
        }
    }
}


if (testseed) {
    for (seed in seedx) {
        for (s in seedmle) {

            nm <- MW214
            set.seed(2019+seed); x <- rnorMmix(500,nm) # set size fixed to 500
            st <- system.time(
                r <- tryCatch({set.seed(2019+s); fitnMm(x, k=1:7, models=1:10,
                                                 trafo="clr1", maxit=1e4)},
                              error = identity)
                )

            bics <- BIC(r)

            sFilepdf <- sprintf("bic_fit_MW214_mleseed=%d_seed=%d.pdf",
                                s, seed)
            pdf(file=file.path(save_dir,sFilepdf))
            plot(r, name=sFilepdf)
            dev.off()

            cat("result of fit: "); str(r, max=1)
            sFile <- sprintf("fit_MW214_mleseed=%d_seed=%d.rds",
                             s, seed)
            cat("--> saving to file:", sFile, "\n")
            ## save both time measurement *and* fit 
            saveRDS(list(st = st, fit = r),
                    file=file.path(save_dir, sFile))
        }
    }
}
