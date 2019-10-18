## 2time had mistake in it that produced incorrect plots

nmmdir <- normalizePath("~/ethz/BA/norMmix.Rcheck/")
savdir <- normalizePath("~/ethz/BA/Rscripts/2time")
stopifnot(dir.exists(nmmdir), dir.exists(savdir))
library(norMmix, lib.loc=nmmdir)
library(mclust)


# for naming purposes
nmnames <- c("MW214", "MW34", "MW51")
sizenames <- c("500", "1000", "2000")
files <- list.files(savdir, pattern=".rds")

fillis <- list()
for (i in seq_along(sizenames)) {
    for (j in seq_along(nmnames)) {
        # for lack of AND matching, OR match everything else and invert
        ret <- grep(paste(sizenames[-i], nmnames[-j], sep="|", collapse="|"), 
                    files, value=TRUE, invert=TRUE)
        fillis[[paste0(nmnames[j], sizenames[i])]] <- ret
    }
}

epfl(fillis, savdir)
