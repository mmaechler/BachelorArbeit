## test effect of size on clustering accuracy

GH_BA_dir <- normalizePath("~/BachelorArbeit")
save_dir  <- normalizePath("~/4NT/largesize")
stopifnot(dir.exists(GH_BA_dir),
          dir.exists(save_dir))
devtools::load_all(file.path(GH_BA_dir, "norMmix"))
##--------------------------------------------------

####
# MW214 is a very difficult cluster. We would like to see improved
# accuracy for larger sample sizes.
# Since this will take a long time, we will vary nothing but initial
# clustering.
####

inits <- c("clara", "mclVVV")
sizes <- c(500, 1000, 2000, 5000, 10000, 20000, 50000)
seeds <- 1:50
filelist <- vector(mode="character", length=700)
i <- 0 # counter for filelist

for (seed in seeds) {
    for (size in sizes) {
        for (ini in inits) {
            set.seed(2019+seed); x <- rnorMmix(size, MW214)
            st <- system.time(
                r <- tryCatch(fitnMm(x, k=1:7, models=1:10,
                                          ini=ini, maxit=1e4),
                              error = identity)
                )
            models <- c("EII","VII","EEI","VEI","EVI",
                        "VVI","EEE","VEE","EVV","VVV")
            mr <- mclust::Mclust(x, G=1:7, modelNames=models)
            sFile <- sprintf("fit_MW214_n=%d_%s_seed=%d.rds",
                             size, ini, seed)
            cat("--> saving to file:", sFile, "\n")
            saveRDS(list(st = st, fit = r, mclust=mr),
                    file=file.path(save_dir, sFile))
            i <- i+1
            filelist[i] <- sFile
        }
    }
}

saveRDS(filelist, file=file.path(save_dir, "filelist"))
