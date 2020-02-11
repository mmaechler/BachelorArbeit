## difference of init on clustering

GH_BA_dir <- normalizePath("~/BachelorArbeit")
save_dir  <- normalizePath("~/4NT/smallinit")
stopifnot(dir.exists(GH_BA_dir),
          dir.exists(save_dir))
devtools::load_all(file.path(GH_BA_dir, "norMmix"))
library(mclust)
##--------------------------------------------------

####
# Intuition here is that clara creates a 'patchwork' of clusters
# whereas agglomerative clustering might find 'pockets' of clusters
# inside larger clusters.
# see https://scikit-learn.org/stable/modules/clustering.html#clustering
# for nice graphic showing agglomerative hierarchical clustering
# behaviour.
# On the other hand, ahc might make mistakes with closely parallel or
# homogenous patterns
####

####
# Setup:
# in ahc's favour: MW24, MW214
# against ahc: MW210, MW29,
# difficult for both: MW28
# easy for both: MW213
#
# Since these are all low dimensional vary over larger sizes
# 500, 1000, 2000, 5000
####

inits <- c("clara", "mclVVV")
cases <- c("MW24", "MW214", "MW210", "MW29", "MW28", "MW213")
sizes <- c(500, 1000, 2000, 5000)
seeds <- 1:50
filelist <- vector(mode="character", length=2400)
i <- 0

for (size in sizes) {
    for (case in cases) {
        nm <- get(case ,"package:norMmix")

        for (ini in inits) {
            for (seed in seeds) {
                set.seed(2019+seed); x <- rnorMmix(size, nm)
                st <- system.time(
                    r <- tryCatch(fitnMm(x, k=1:7, models=1:10,
                                              ini=ini, maxit=1e4),
                                  error = identity)
                    )

                models <- c("EII","VII","EEI","VEI","EVI",
                            "VVI","EEE","VEE","EVV","VVV")
                mr <- mclust::Mclust(x, G=1:7, modelNames=models)

                sFile <- sprintf("fit_%s_n=%d_%s_seed=%d.rds",
                                 case, size, ini, seed)
                cat("--> saving to file:", sFile, "\n")
                saveRDS(list(st = st, fit = r, mclust=mr),
                        file=file.path(save_dir, sFile))
                i <- i+1
                filelist[i] <- sFile
            }
        }
    }
}

saveRDS(filelist, file=file.path(save_dir, "filelist"))
