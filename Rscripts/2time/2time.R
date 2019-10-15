## Intent: analyse time as function of p,k,n

nmmdir <- normalizePath("~/BachelorArbeit/norMmix.Rcheck/")
savdir <- normalizePath("~/BachelorArbeit/Rscripts/2time")
stopifnot(dir.exists(nmmdir), dir.exists(savdir))
library(norMmix, lib.loc=nmmdir)
library(mclust)

## at n=500,p=2 can do about 250xfitnMm(x,1:10) in 24h
seeds <- 1:10
sizes <- c(500, 1000, 2000)
nmm <- list(MW214, MW34, MW51)
## => about 100 cases

# for naming purposes
nmnames <- c("MW214", "MW34", "MW51")
sizenames <- c("500", "1000", "2000")
files <- vector(mode="character")

for (nm in 1:3) {
    for (size in sizes) {
    set.seed(2019); x <- rnorMmix(size, nmm[[nm]])
        for (seed in seeds) {
            set.seed(2019+seed)
            r <- tryCatch(fitnMm(x, k=1:8,
                                 optREPORT=1e4, maxit=1e4),
                          error = identity)
            filename <- sprintf("%s_size=%0.4d_seed=%0.2d.rds",
                                nmnames[nm], size, seed)
            files <- append(files, filename)
            cat("===> saving to file:", filename, "\n")
            saveRDS(list(fit=r), file=file.path(savdir, filename))
        }
    }
}

fillis <- list()
for (i in seq_along(sizes)) {
    for (j in seq_along(nmnames)) {
        # for lack of AND matching, OR match everything else and invert
        ret <- grep(paste(sizenames[-i], nmnames[-j], sep="|"), 
                    files, value=TRUE, invert=TRUE)
        fillis[[paste0(sizenames[i], nmnames[j])]] <- ret
    }
}

epfl(fillis, savdir)
