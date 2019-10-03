savdir <- normalizePath("~/BScThesis/4MM/small-tri")
nmmdir <- normalizePath("~/BachelorArbeit")
library(mclust)
devtools::load_all(file.path(nmmdir, "norMmix"))
source("~/BachelorArbeit/Rscripts/adafuncs.R")
source("~/BachelorArbeit/Rscripts/tempmassbic.R")

seeds <- 1:50
inits <- c("clara", "mclVVV")
files_clara <- vector()
files_mclVVV <- vector()

for (ini in inits) {
    for (seed in seeds) {
        sFile <- sprintf("fit_MW215_n=1000_%s_seed=%d.rds",
                         ini, seed)
        switch(ini,
               "clara" = files_clara <- append(files_clara, sFile),
               "mclVVV" = files_mclVVV <- append(files_mclVVV, sFile)
               )
    }
}

print(files_clara)

bicsnmm_clara <- tmassbic(files_clara, savdir)
bicsnmm_mclVVV <- tmassbic(files_mclVVV, savdir)
bicsmcl_clara <- tmassbicm(files_clara, savdir)
bicsmcl_mclVVV <- tmassbicm(files_mclVVV, savdir)

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
