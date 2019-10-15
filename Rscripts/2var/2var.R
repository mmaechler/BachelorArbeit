## Intent: analyse various datasets, not normal mixtures

nmmdir <- normalizePath("~/BachelorArbeit/norMmix.Rcheck/")
savdir <- normalizePath("~/BachelorArbeit/Rscripts/2var")
stopifnot(dir.exists(nmmdir), dir.exists(savdir))
library(norMmix, lib.loc=nmmdir)
library(mclust)

## at n=500,p=2 can do about 250xfitnMm(x,1:10) in 24h
seeds <- 1:25

irt <- data.matrix(iris[,-5]) # irt= IRis Truncated
iri <- data.matrix(iris)
data(loss, package="copula")
los <- data.matrix(loss)
dats <- list(irt, iri, los)

# for naming purposes
datnames <- c("tr_iris", "iris", "loss")
files <- vector(mode="character")

for (i in seq_along(dats)) {
    for (seed in seeds) {
        set.seed(2019+seed)
        r <- tryCatch(fitnMm(dats[[i]], k=1:8,
                             optREPORT=1e4, maxit=1e4),
                      error = identity)
        filename <- sprintf("%s_seed=%0.2d.rds",
                            datnames[i], seed)
        files <- append(files, filename)
        cat("===> saving to file:", filename, "\n")
        saveRDS(list(fit=r), file=file.path(savdir, filename))
    }
}

fillis <- list()
for (i in datnames) {
        ret <- grep(i, files, value=TRUE)
        fillis[[i]] <- ret
}

epfl(fillis, savdir)
