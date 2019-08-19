GH_BA_dir <- normalizePath("~/BachelorArbeit")
save_dir  <- normalizePath("~/BScThesis/4MM")
##-------------------------------------------------
stopifnot(dir.exists(GH_BA_dir),
          dir.exists(save_dir))

## MM: Dies ist *NICHT* die Art ein Paket zu laden .. 
##     aber "good enough for now":
devtools::load_all(file.path(GH_BA_dir, "norMmix"))

## Es gibt gar kein '/data' subdirectory !!
(rdat <- list.files(file.path(GH_BA_dir, "norMmix/data"), pattern="mw.*"))

## Die Datensätze sind *Teil* der Paket Objekte :
(MWdat <- Filter(function(.) is.norMmix(get(., "package:norMmix")),
                 ls("package:norMmix", pattern = "^MW[1-9]")))

trafos <- c("clr1", "logit")
inits <- c("clara", "mclVV") # <- neue Namen

seed.n <- 1

for (nm in MWdat) {
    cat("work on norMmix model ",nm,"\n")
    nMm <- get(nm ,"package:norMmix")
    for(n in c(100, 200, 500, 1000)) {
        cat(" seed = ",seed.n,"\n")
        x <- rnorMmix(n, nMm)
        cat("dataset x:"); str(x)
        
        for(trafo in trafos) {
            cat(">--> trafo = ", dQuote(trafo),"\n")
            for(ini in inits) {
                set.seed(seed.n) # for reproducibility
                st <- system.time(
                    r <- tryCatch(fit.norMmix(x, k=1:7, models=1:10,
                                              trafo=trafo, ini=ini, maxiter=300),
                                  error = identity)
                )
                cat("result of fit: "); str(r, max=1)
                sFile <- sprintf("fit_%s_n=%d_%s_%s.rds",
                                 nm, n, trafo, ini)
                cat("--> saving to file:", sFile, "\n")
                ## save both time measurement *and* fit 
                saveRDS(list(st = st, fit = r),
                        file=file.path(save_dir, sFile))
                seed.n <- seed.n + 1
            }
        }
        cat("\n---end-of-dataset---------------------\n")
    }
    cat("\n***============================================***\n")
}








### Nomenklatur:

## Daten in norMmix/data/  :

## xx.yy.zz.RDS

# x: Ursprungsdatensatz z.B: MW-- norMmix Objekte
# y: sample size: z.B 1e4
# z: "korrekte" cluster und model


## output von diesem Script :

## aa.bb.cc.d

# a: trafo
# b: ini
# c: range d.h. k=1:7 models=1:10 -> k17m110
# d: objekt auf dem operiert wird z.B. mw.1e3.22EEE.RDS
