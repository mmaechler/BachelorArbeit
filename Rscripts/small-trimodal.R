## try MW215

nmmdir <- normalizePath("")
savdir <- normalizePath("")
stopifnot(dir.exists(GH_BA_dir),
          dir.exists(save_dir))
library(mclust)
devtools::load_all(file.path(nmdir, "norMmix"))

seeds <- 1:50
inits <- c("clara", "mclVVV")
set.seed(2019); x <- rnorMmix(1000, MW215)
## correct: 3, VEE
files_clara <- vector()
files_mclVVV <- vector()

for (ini in inits) {
    for (seed in seeds) {
        r <- fitnMm(x, k=1:7, ini=ini, maxit=1e4)
        cat("result of fit: "); str(r, max=1)
        sFile <- sprintf("fit_MW215_n=1000_%s_seed=%d.rds",
                         ini, seed)
        cat("--> saving to file:", sFile, "\n")
        ## save both time measurement *and* fit 
        saveRDS(r,file=file.path(savdir, sFile))
        switch(ini,
               "clara" = append(files_clara, sFile),
               "mclVVV" = append(files_mclVVV, sFile)
               )
    }
}

bicsnmm_clara <- massbic(files_clara, savdir)
bicsnmm_mclVVV <- massbic(files_mclVVV, savdir)
bicsmcl_clara <- massbicm(files_clara, savdir)
bicsmcl_mclVVV <- massbicm(files_mclVVV, savdir)

saveRDS(bicsmcl_clara, file=file.path(savdir, "bicsmcl_clara.rds"))
saveRDS(bicsmcl_mclVVV, file=file.path(savdir, "bicsmcl_mclVVV.rds"))

pdf(file=file.path(savdir, "nmm_clara.pdf"))
massplot(bicsnmm_clara, main="nmm_clara, 50 seeds"); dev.off()
pdf(file=file.path(savdir, "mcl_clara.pdf"))
massplot(bicsmcl_clara, main="mcl_clara, 50 seeds"); dev.off()
pdf(file=file.path(savdir, "nmm_mclVV.pdf"))
massplot(bicsnmm_mclVVV, main="nmm_mclVV, 50 seeds"); dev.off()
pdf(file=file.path(savdir, "mcl_mclVV.pdf"))
massplot(bicsmcl_mclVVV, main="nmm_mclVV, 50 seeds"); dev.off()

pdf(file=file.path(savdir, "comp_nmm_clara.pdf"))
massplot(bicsnmm_clara, bicsmcl_clara); dev.off()
pdf(file=file.path(savdir, "comp_nmm_mclVV.pdf"))
massplot(bicsnmm_mclVVV, bicsmcl_mclVVV); dev.off()

cat("done. \n")
