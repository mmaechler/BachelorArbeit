## nm2

GH_BA_dir <- normalizePath("~/ethz/BA")
save_dir  <- normalizePath("~/ethz/BA/Rscripts")
##-------------------------------------------------
stopifnot(dir.exists(GH_BA_dir),
          dir.exists(save_dir))

## MM: Dies ist *NICHT* die Art ein Paket zu laden .. 
##     aber "good enough for now":
devtools::load_all(file.path(GH_BA_dir, "norMmix"))

##set.seed(2019); x <- rnorMmix(300, MW27)
##
##ans <- fit.norMmix(x)
##
##pdf(file=file.path(save_dir, "temp.pdf"))
##plot(ans, name="test")
##dev.off()
##
##plot(ans, name="test", plotbest=TRUE)



(rdat <- list.files(file.path(GH_BA_dir, "norMmix/data"), pattern="mw.*"))

(MWdat <- Filter(function(.) is.norMmix(get(., "package:norMmix")),
                 ls("package:norMmix", pattern = "^MW[1-9]")))

trafos <- c("clr1", "logit")
inits <- c("clara", "mclVVV")
lls <- c("nmm", "mvt")
size <- c(400,500)
seeds <- 1:50

for (dat in MWdat) {
    for (trafo in trafos) {
        for (ini in inits) {
            for (ll in lls) {
                for (s in size) {
                    for (seed in seeds) {

                        nm <- get(dat ,"package:norMmix")
                        set.seed(2019+seed); x <- rnorMmix(s,nm)
                        st <- system.time(
                            r <- tryCatch(fit.norMmix(x, k=1:7, models=1:10,
                                                      trafo=trafo, ini=ini, maxit=1e4),
                                          error = identity)
                            )

                        bics <- BIC(r)

                        sFilepdf <- sprintf("bic_fit_%s_n=%d_%s_%s_seed=%d.pdf",
                                            dat, s, trafo, ini, seed)
                        pdf(file=file.path(save_dir,sFilepdf))
                        plot(r, name=sFilepdf)
                        dev.off()


                        #sFilepdf2 <- sprintf("best_fit_%s_n=%d_%s_%s_seed=%d.pdf",
                        #                     dat, s, trafo, ini, seed)
                        #pdf(file=file.path(save_dir,sFilepdf2))
                        #plot(r, name=sFilepdf, plotbest=TRUE)
                        #dev.off()

                        cat("result of fit: "); str(r, max=1)
                        sFile <- sprintf("fit_%s_n=%d_%s_%s_seed=%d.rds",
                                         dat, s, trafo, ini, seed)
                        cat("--> saving to file:", sFile, "\n")
                        ## save both time measurement *and* fit 
                        saveRDS(list(st = st, fit = r),
                                file=file.path(save_dir, sFile))
                    }
                }
            }
        }
    }
}
