## dong a small run with smi data


GH_BA_dir <- normalizePath("~/BachelorArbeit")
save_dir  <- normalizePath("~/")
sourcedir <- normalizePath("~/BachelorArbeit/Rscripts")
stopifnot(dir.exists(GH_BA_dir),
          dir.exists(save_dir))
devtools::load_all(file.path(GH_BA_dir, "norMmix"))
source(file.path(sourcedir,"adafuncs.R"))
##--------------------------------------------------


## data sets:
data(SMI.12, package="copula")
smi <- SMI.12

## variables:
inits <- c("clara", "mclVVV")
seeds <- 1:50
# inits to study behaviour under different inits
# and seed to get statistical results

## 5/parcond(smi, k=7, model="VVV") = 57.30
## 5/parcond(smi, k=6, model="VVV") = 49.11
## this should inform how many simulations
## need to be done to get signal to noise high
## If we want s/n=50, need k<=6


for (ini in inits) {
    files <- character(length(seeds))

    for (s in seeds) {
        set.seed(2019+s)
        st <- system.time(
            r <- tryCatch(fit.norMmix(smi, k=1:8, models=1:10, ini=ini,
                                      optREPORT=20, maxit=1e4),
                          error = identity)
            )

        bics <- BIC(r)

        sFilepdf <- sprintf("bic_smi_%s_seed=%d.pdf",
                            ini, s)
        pdf(file=file.path(save_dir,sFilepdf))
        plot(r, name=sFilepdf)
        dev.off()

        cat("result of fit: "); str(r, max=1)
        sFile <- sprintf("fit_smi_%s_seed=%d.rds",
                         ini, s)
        cat("--> saving to file:", sFile, "\n")
        files[s] <- sFile

        summa <- list(systemtime=st, bic=bics)
        saveRDS(list(s = summa, fit = r),
                file=file.path(save_dir, sFile))
    }

    saveRDS(files, file=file.path(save_dir,sprintf("filenames_%s.rds", ini)))
    g <- massbic(files, savedir)$v

    pdf(file=file.path(save_dir, sprintf("summary_%s.pdf", ini)),20,16)
    op <- par(mfrow=c(4,5))
    for (i in 1:10) {
        matplot(g[,i,], lty=1, col=rainbow(10)[i], main=models[i], type="l")
    }
    for (i in 1:10) {
        boxplot(t(g[,i,]), lty=1, col=rainbow(10)[i], main=models[i], type="l")
    }
    par(op)
    dev.off()
} 
