GH_BA_dir <- normalizePath("~/ethz/BA/")
save_dir  <- normalizePath("~/")
##--------------------------------------------------
stopifnot(dir.exists(GH_BA_dir),
          dir.exists(save_dir))


devtools::load_all(file.path(GH_BA_dir, "norMmix"))


## data sets:
data(SMI.12, package="copula")
smi <- SMI.12

# options to vary over
dd <- "smi.12"
trafo <- "clr1"
init <- "clara"
lls <- c("nmm", "mvt")


seed.n <- 1

cat(">--> init = ", dQuote(init), "\n")
for (ll in lls) {
    cat(">--> loglikelihood = ", dQuote(ll), "\n")
    set.seed(seed.n)
    st <- system.time(
        r <- tryCatch(
            fit.norMmix(smi, k=1:9, models=1:10,
                        trafo=trafo, ini=init, ll=ll,
                        maxit=300, epsilon=0.1),
            error = identity)
    )
    cat("result of fit: "); str(r, max.level=1)
    sFile <- sprintf("fit_%s_%s_%s_%s.rds",
                     dd, trafo, ini, ll)
    cat("--saving to file:", sFile, "\n")
    saveRDS(list(fit=r, st=st),
            file=file.path(save_dir, sFile))
    seed.n <- seed.n+1
}
cat("\n---end-of-dataset---\n")
cat("\n***--------------------------------------------***\n")



