library(norMmix, lib.loc="~/ethz/BA/norMmix.Rcheck/")
# change this dir to whereever the simulations are saved
mainsav <- normalizePath("~/ethz/BA/Rscripts/")

savdir <- file.path(mainsav, "2time")
filelist <- list.files(savdir, pattern=".rds")
filelist <- grep("mcl.rds", filelist, invert=TRUE, value=TRUE)

size05 <- grep("size=0500", filelist, value=TRUE)
size10 <- grep("size=1000", filelist, value=TRUE)
size20 <- grep("size=2000", filelist, value=TRUE)
    
size05mw34 <- grep("MW34", size05, value=TRUE)
size05mw51 <- grep("MW51", size05, value=TRUE)
size05mw214 <- grep("MW214", size05, value=TRUE)

size10mw34 <- grep("MW34", size10, value=TRUE)
size10mw51 <- grep("MW51", size10, value=TRUE)
size10mw214 <- grep("MW214", size10, value=TRUE)

size20mw34 <- grep("MW34", size20, value=TRUE)
size20mw51 <- grep("MW51", size20, value=TRUE)
size20mw214 <- grep("MW214", size20, value=TRUE)


f <- readRDS(file=file.path(savdir,size05mw34[1]))

## apply clara

tmp <- function(o, name) {
    v <- readRDS(file.path(savdir,o))
    x <- v$fit$x
    size <- v$fit$n
    r <- tryCatch(fitnMm(x, k=1:8, ini="mclVVV",
                         optREPORT=1e4, maxit=1e4),
                  error = identity)
    filename <- sprintf("%s_size=%0.4d_mclVVV.rds",
                        name, size)
    cat("===> saving to file:", filename, "\n")
    saveRDS(list(fit=r), file=file.path(savdir, filename))
}

tmp(size05mw214[3], "MW214")
tmp(size10mw214[3], "MW214")
tmp(size20mw214[3], "MW214")

tmp(size05mw34[3], "MW34")
tmp(size10mw34[3], "MW34")
tmp(size20mw34[3], "MW34")

tmp(size05mw51[3], "MW51")
tmp(size10mw51[3], "MW51")
tmp(size20mw51[3], "MW51")

