## Intent: analyse SMI.12 dataset

nmmdir <- normalizePath("~/ethz/BA/norMmix.Rcheck/")
savdir <-  normalizePath("~/ethz/BA/Rscripts/2time")
stopifnot(dir.exists(nmmdir), dir.exists(savdir))
library(norMmix, lib.loc=nmmdir)
library(mclust)

## at n=500,p=2 can do about 250xfitnMm(x,1:10) in 24h
seeds <- 1:30
inits <- c("clara", "mclVVV")
data(SMI.12, package="copula")
smi <- SMI.12
## => about 100 cases

# for naming purposes
files <- vector(mode="character")


for (ini in inits) {
    for (seed in seeds) {
        set.seed(2019+seed)
        r <- tryCatch(fitnMm(smi, k=1:8, ini=ini
                             optREPORT=1e4, maxit=1e4),
                      error = identity)
        filename <- sprintf("smi_ini=%s_seed=%0.2d.rds",
                            ini, seed)
        files <- append(files, filename)
        cat("===> saving to file:", filename, "\n")
        saveRDS(list(fit=r), file=file.path(savdir, filename))
    }
}

fillis <- list()
for (i in inits) {
    # for lack of AND matching, OR match everything else and invert
    ret <- grep(i, 
                files, value=TRUE)
    fillis[[i]] <- ret
}

epfl(fillis, savdir)
