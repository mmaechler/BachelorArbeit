#### plots of simulations, as shown in meeting of 2019-09-26


## correct dirs for local system

# dir of norMmix package
nmmdir <- "~/ethz/BA/norMmix"
#dir with simulation data
rdsdir <- "~/ethz/BA/Rscripts"
# save dir
savdir <- "~/ethz/BA/Rscripts/sim_plots"
stopifnot(dir.exists(nmmdir),
          dir.exists(rdsdir),
          dir.exists(savdir))


## load tools
devtools::load_all(nmmdir)
source(file.path(rdsdir, "adafuncs.R"))


## sim names
ssiz <- "smallsize"
ssee <- "smallseed"
nmm2 <- "nm2"
# smi2 <- "fit-smi2" too small for consideration
ssmi <- "smallsmi"
sini <- "smallinit"

## sim .rds files
sssiz <- list.files(file.path(rdsdir,ssiz), pattern="*rds")
sssee <- list.files(file.path(rdsdir,ssee), pattern="*rds")
snmm2 <- list.files(file.path(rdsdir,nmm2), pattern="*rds")
    snmm2 <- snmm2[grep("clr1", snmm2)] # prescreening logit files
sssmi <- list.files(file.path(rdsdir,ssmi), pattern="*rds")
ssini <- list.files(file.path(rdsdir,sini), pattern="*rds")

## subdivide files
# sssiz
sssiz_0500 <- sssiz[grep("=500", sssiz)]
sssiz_0600 <- sssiz[grep("=600", sssiz)]
sssiz_0700 <- sssiz[grep("=700", sssiz)]
sssiz_0800 <- sssiz[grep("=800", sssiz)]
sssiz_0900 <- sssiz[grep("=900", sssiz)]
sssiz_1000 <- sssiz[grep("=1000", sssiz)]
sssiz_1100 <- sssiz[grep("=1100", sssiz)]
sssiz_1200 <- sssiz[grep("=1200", sssiz)]
sssiz_1300 <- sssiz[grep("=1300", sssiz)]
sssiz_1400 <- sssiz[grep("=1400", sssiz)]
sssiz_1500 <- sssiz[grep("=1500", sssiz)]

# sssee
sssee_seed01 <- sssee[grep("1.rds", sssee)]
sssee_seed02 <- sssee[grep("2.rds", sssee)]
sssee_seed03 <- sssee[grep("3.rds", sssee)]
sssee_seed04 <- sssee[grep("4.rds", sssee)]
sssee_seed05 <- sssee[grep("5.rds", sssee)]
sssee_seed06 <- sssee[grep("6.rds", sssee)]
sssee_seed07 <- sssee[grep("7.rds", sssee)]
sssee_seed08 <- sssee[grep("8.rds", sssee)]
sssee_seed09 <- sssee[grep("9.rds", sssee)]
sssee_seed10 <- sssee[grep("0.rds", sssee)]

# snmm2
snmm2_MW210 <- snmm2[grep("MW210_", snmm2)]
    snmm2_MW210_400 <- snmm2_MW210[grep("=400", snmm2_MW210)]
        snmm2_MW210_400_clara <- snmm2_MW210_400[grep("clara", snmm2_MW210_400)]
        snmm2_MW210_400_mclVVV <- snmm2_MW210_400[grep("mclV", snmm2_MW210_400)]
    snmm2_MW210_500 <- snmm2_MW210[grep("=500", snmm2_MW210)]
        snmm2_MW210_500_clara <- snmm2_MW210_500[grep("clara", snmm2_MW210_500)]
        snmm2_MW210_500_mclVVV <- snmm2_MW210_500[grep("mclV", snmm2_MW210_500)]
snmm2_MW21 <- snmm2[grep("MW21_", snmm2)]
    snmm2_MW21_400 <- snmm2_MW21[grep("=400", snmm2_MW21)]
        snmm2_MW21_400_clara <- snmm2_MW21_400[grep("clara", snmm2_MW21_400)]
        snmm2_MW21_400_mclVVV <- snmm2_MW21_400[grep("mclV", snmm2_MW21_400)]
    snmm2_MW21_500 <- snmm2_MW21[grep("=500", snmm2_MW21)]
        snmm2_MW21_500_clara <- snmm2_MW21_500[grep("clara", snmm2_MW21_500)]
        snmm2_MW21_500_mclVVV <- snmm2_MW21_500[grep("mclV", snmm2_MW21_500)]

# sssmi no variables to separate

# ssini
ssini_MW24 <- ssini[grep("MW24", ssini)]
    ssini_MW24_1000 <- ssini_MW24[grep("1000", ssini_MW24)]
        ssini_MW24_1000_clara <- ssini_MW24_1000[grep("clara", ssini_MW24_1000)]
        ssini_MW24_1000_mclVVV <- ssini_MW24_1000[grep("mclV", ssini_MW24_1000)]
    ssini_MW24_500 <- ssini_MW24[grep("500", ssini_MW24)]
        ssini_MW24_500_clara <- ssini_MW24_500[grep("clara", ssini_MW24_500)]
        ssini_MW24_500_mclVVV <- ssini_MW24_500[grep("mclV", ssini_MW24_500)]
ssini_MW28 <- ssini[grep("MW28", ssini)]
    ssini_MW28_500 <- ssini_MW28[grep("500", ssini_MW28)]
        ssini_MW28_500_clara <- ssini_MW28_500[grep("clara", ssini_MW28_500)]
        ssini_MW28_500_mclVVV <- ssini_MW28_500[grep("mclV", ssini_MW28_500)]
    ssini_MW28_1000 <- ssini_MW28[grep("1000", ssini_MW28)]
        ssini_MW28_1000_clara <- ssini_MW28_1000[grep("clara", ssini_MW28_1000)]
        ssini_MW28_1000_mclVVV <- ssini_MW28_1000[grep("mclV", ssini_MW28_1000)]
ssini_MW29 <- ssini[grep("MW29", ssini)]
    ssini_MW29_500 <- ssini_MW29[grep("500", ssini_MW29)]
        ssini_MW29_500_clara <- ssini_MW29_500[grep("clara", ssini_MW29_500)]
        ssini_MW29_500_mclVVV <- ssini_MW29_500[grep("mclV", ssini_MW29_500)]
    ssini_MW29_1000 <- ssini_MW29[grep("1000", ssini_MW29)]
        ssini_MW29_1000_clara <- ssini_MW29_1000[grep("clara", ssini_MW29_1000)]
        ssini_MW29_1000_mclVVV <- ssini_MW29_1000[grep("mclV", ssini_MW29_1000)]
ssini_MW210 <- ssini[grep("MW210", ssini)]
    ssini_MW210_500 <- ssini_MW210[grep("500", ssini_MW210)]
        ssini_MW210_500_clara <- ssini_MW210_500[grep("clara", ssini_MW210_500)]
        ssini_MW210_500_mclVVV <- ssini_MW210_500[grep("mclV", ssini_MW210_500)]
    ssini_MW210_1000 <- ssini_MW210[grep("1000", ssini_MW210)]
        ssini_MW210_1000_clara <- ssini_MW210_1000[grep("clara", ssini_MW210_1000)]
        ssini_MW210_1000_mclVVV <- ssini_MW210_1000[grep("mclV", ssini_MW210_1000)]
ssini_MW213 <- ssini[grep("MW213", ssini)]
    ssini_MW213_500 <- ssini_MW213[grep("500", ssini_MW213)]
        ssini_MW213_500_clara <- ssini_MW213_500[grep("clara", ssini_MW213_500)]
        ssini_MW213_500_mclVVV <- ssini_MW213_500[grep("mclV", ssini_MW213_500)]
    ssini_MW213_1000 <- ssini_MW213[grep("1000", ssini_MW213)]
        ssini_MW213_1000_clara <- ssini_MW213_1000[grep("clara", ssini_MW213_1000)]
        ssini_MW213_1000_mclVVV <- ssini_MW213_1000[grep("mclV", ssini_MW213_1000)]
ssini_MW214 <- ssini[grep("MW214", ssini)]
    ssini_MW214_500 <- ssini_MW214[grep("500", ssini_MW214)]
        ssini_MW214_500_clara <- ssini_MW214_500[grep("clara", ssini_MW214_500)]
        ssini_MW214_500_mclVVV <- ssini_MW214_500[grep("mclV", ssini_MW214_500)]
    ssini_MW214_1000 <- ssini_MW214[grep("1000", ssini_MW214)]
        ssini_MW214_1000_clara <- ssini_MW214_1000[grep("clara", ssini_MW214_1000)]
        ssini_MW214_1000_mclVVV <- ssini_MW214_1000[grep("mclV", ssini_MW214_1000)]

#### get bic arrays from files        

        
cat("massbic from here \n")        
        


nmsssiz <-  massbic(sssiz, file.path(rdsdir, ssiz))
mmsssiz <-  massbicm(sssiz, file.path(rdsdir, ssiz))
nmsssee <-  massbic(sssee, file.path(rdsdir, ssee))
mmsssee <-  massbicm(sssee, file.path(rdsdir, ssee))
nmsnmm2 <-  massbic(snmm2, file.path(rdsdir, nmm2))
mmsnmm2 <-  massbicm(snmm2, file.path(rdsdir, nmm2))
nmsssmi <-  massbic(sssmi, file.path(rdsdir, ssmi))
mmsssmi <-  massbicm(sssmi, file.path(rdsdir, ssmi))
nmssini <-  massbic(ssini, file.path(rdsdir, sini))
mmssini <-  massbicm(ssini, file.path(rdsdir, sini))

cat("full bics done \n")



nmsssiz_0500 <- massbic(sssiz_0500, file.path(rdsdir, ssiz ))
mmsssiz_0500 <- massbicm(sssiz_0500, file.path(rdsdir, ssiz ))
nmsssiz_0600 <- massbic(sssiz_0600, file.path(rdsdir, ssiz ))
mmsssiz_0600 <- massbicm(sssiz_0600, file.path(rdsdir, ssiz ))
nmsssiz_0700 <- massbic(sssiz_0700, file.path(rdsdir, ssiz ))
mmsssiz_0700 <- massbicm(sssiz_0700, file.path(rdsdir, ssiz ))
nmsssiz_0800 <- massbic(sssiz_0800, file.path(rdsdir, ssiz ))
mmsssiz_0800 <- massbicm(sssiz_0800, file.path(rdsdir, ssiz ))
nmsssiz_0900 <- massbic(sssiz_0900, file.path(rdsdir, ssiz ))
mmsssiz_0900 <- massbicm(sssiz_0900, file.path(rdsdir, ssiz ))
nmsssiz_1000 <- massbic(sssiz_1000, file.path(rdsdir, ssiz ))
mmsssiz_1000 <- massbicm(sssiz_1000, file.path(rdsdir, ssiz ))
nmsssiz_1100 <- massbic(sssiz_1100, file.path(rdsdir, ssiz ))
mmsssiz_1100 <- massbicm(sssiz_1100, file.path(rdsdir, ssiz ))
nmsssiz_1200 <- massbic(sssiz_1200, file.path(rdsdir, ssiz ))
mmsssiz_1200 <- massbicm(sssiz_1200, file.path(rdsdir, ssiz ))
nmsssiz_1300 <- massbic(sssiz_1300, file.path(rdsdir, ssiz ))
mmsssiz_1300 <- massbicm(sssiz_1300, file.path(rdsdir, ssiz ))
nmsssiz_1400 <- massbic(sssiz_1400, file.path(rdsdir, ssiz ))
mmsssiz_1400 <- massbicm(sssiz_1400, file.path(rdsdir, ssiz ))
nmsssiz_1500 <- massbic(sssiz_1500, file.path(rdsdir, ssiz ))
mmsssiz_1500 <- massbicm(sssiz_1500, file.path(rdsdir, ssiz ))

cat("ssize done \n")

nmsssee_seed01 <- massbic(sssee_seed01, file.path(rdsdir, ssee))
mmsssee_seed01 <- massbicm(sssee_seed01, file.path(rdsdir, ssee))
nmsssee_seed02 <- massbic(sssee_seed02, file.path(rdsdir, ssee))
mmsssee_seed02 <- massbicm(sssee_seed02, file.path(rdsdir, ssee))
nmsssee_seed03 <- massbic(sssee_seed03, file.path(rdsdir, ssee))
mmsssee_seed03 <- massbicm(sssee_seed03, file.path(rdsdir, ssee))
nmsssee_seed04 <- massbic(sssee_seed04, file.path(rdsdir, ssee))
mmsssee_seed04 <- massbicm(sssee_seed04, file.path(rdsdir, ssee))
nmsssee_seed05 <- massbic(sssee_seed05, file.path(rdsdir, ssee))
mmsssee_seed05 <- massbicm(sssee_seed05, file.path(rdsdir, ssee))
nmsssee_seed06 <- massbic(sssee_seed06, file.path(rdsdir, ssee))
mmsssee_seed06 <- massbicm(sssee_seed06, file.path(rdsdir, ssee))
nmsssee_seed07 <- massbic(sssee_seed07, file.path(rdsdir, ssee))
mmsssee_seed07 <- massbicm(sssee_seed07, file.path(rdsdir, ssee))
nmsssee_seed08 <- massbic(sssee_seed08, file.path(rdsdir, ssee))
mmsssee_seed08 <- massbicm(sssee_seed08, file.path(rdsdir, ssee))
nmsssee_seed09 <- massbic(sssee_seed09, file.path(rdsdir, ssee))
mmsssee_seed09 <- massbicm(sssee_seed09, file.path(rdsdir, ssee))
nmsssee_seed10 <- massbic(sssee_seed10, file.path(rdsdir, ssee))
mmsssee_seed10 <- massbicm(sssee_seed10, file.path(rdsdir, ssee))

cat("sseed done \n")

nmsssmi <- massbic(sssmi, file.path(rdsdir, ssmi))
mmsssmi <- massbicm(sssmi, file.path(rdsdir, ssmi))

cat("ssmi done \n")

nmsnmm2_MW210 <-  massbic(snmm2_MW210, file.path(rdsdir, nmm2))
mmsnmm2_MW210 <-  massbicm(snmm2_MW210, file.path(rdsdir, nmm2))
    nmsnmm2_MW210_400 <- massbic(snmm2_MW210_400, file.path(rdsdir, nmm2))
    mmsnmm2_MW210_400 <- massbicm(snmm2_MW210_400, file.path(rdsdir, nmm2))
        nmsnmm2_MW210_400_clara <- massbic(snmm2_MW210_400_clara, file.path(rdsdir, nmm2))
        mmsnmm2_MW210_400_clara <- massbicm(snmm2_MW210_400_clara, file.path(rdsdir, nmm2))
        nmsnmm2_MW210_400_mclVVV <- massbic(snmm2_MW210_400_mclVVV, file.path(rdsdir, nmm2))
        mmsnmm2_MW210_400_mclVVV <- massbicm(snmm2_MW210_400_mclVVV, file.path(rdsdir, nmm2))
    nmsnmm2_MW210_500 <- massbic(snmm2_MW210_500, file.path(rdsdir, nmm2))
    mmsnmm2_MW210_500 <- massbicm(snmm2_MW210_500, file.path(rdsdir, nmm2))
        nmsnmm2_MW210_500_clara <- massbic(snmm2_MW210_500_clara, file.path(rdsdir, nmm2))
        mmsnmm2_MW210_500_clara <- massbicm(snmm2_MW210_500_clara, file.path(rdsdir, nmm2))
        nmsnmm2_MW210_500_mclVVV <- massbic(snmm2_MW210_500_mclVVV, file.path(rdsdir, nmm2))
        mmsnmm2_MW210_500_mclVVV <- massbicm(snmm2_MW210_500_mclVVV, file.path(rdsdir, nmm2))
nmsnmm2_MW21 <- massbic(snmm2_MW21, file.path(rdsdir, nmm2))
mmsnmm2_MW21 <- massbicm(snmm2_MW21, file.path(rdsdir, nmm2))
    nmsnmm2_MW21_400 <- massbic(snmm2_MW21_400, file.path(rdsdir, nmm2))
    mmsnmm2_MW21_400 <- massbicm(snmm2_MW21_400, file.path(rdsdir, nmm2))
        nmsnmm2_MW21_400_clara <- massbic(snmm2_MW21_400_clara, file.path(rdsdir, nmm2))
        mmsnmm2_MW21_400_clara <- massbicm(snmm2_MW21_400_clara, file.path(rdsdir, nmm2))
        nmsnmm2_MW21_400_mclVVV <- massbic(snmm2_MW21_400_mclVVV, file.path(rdsdir, nmm2))
        mmsnmm2_MW21_400_mclVVV <- massbicm(snmm2_MW21_400_mclVVV, file.path(rdsdir, nmm2))
    nmsnmm2_MW21_500 <- massbic(snmm2_MW21_500, file.path(rdsdir, nmm2))
    mmsnmm2_MW21_500 <- massbicm(snmm2_MW21_500, file.path(rdsdir, nmm2))
        nmsnmm2_MW21_500_clara <- massbic(snmm2_MW21_500_clara, file.path(rdsdir, nmm2))
        mmsnmm2_MW21_500_clara <- massbicm(snmm2_MW21_500_clara, file.path(rdsdir, nmm2))
        nmsnmm2_MW21_500_mclVVV <- massbic(snmm2_MW21_500_mclVVV, file.path(rdsdir, nmm2))
        mmsnmm2_MW21_500_mclVVV <- massbicm(snmm2_MW21_500_mclVVV, file.path(rdsdir, nmm2))

cat("nm2 done \n")


nmssini_MW24 <- massbic(ssini_MW24, file.path(rdsdir, sini))
mmssini_MW24 <- massbicm(ssini_MW24, file.path(rdsdir, sini))
    nmssini_MW24_1000 <- massbic(ssini_MW24_1000, file.path(rdsdir, sini))
    mmssini_MW24_1000 <- massbicm(ssini_MW24_1000, file.path(rdsdir, sini))
        nmssini_MW24_1000_clara <- massbic(ssini_MW24_1000_clara, file.path(rdsdir, sini))
        mmssini_MW24_1000_clara <- massbicm(ssini_MW24_1000_clara, file.path(rdsdir, sini))
        nmssini_MW24_1000_mclVVV <- massbic(ssini_MW24_1000_mclVVV, file.path(rdsdir, sini))
        mmssini_MW24_1000_mclVVV <- massbicm(ssini_MW24_1000_mclVVV, file.path(rdsdir, sini))
    nmssini_MW24_500 <- massbic(ssini_MW24_500, file.path(rdsdir, sini))
    mmssini_MW24_500 <- massbicm(ssini_MW24_500, file.path(rdsdir, sini))
        nmssini_MW24_500_clara <- massbic(ssini_MW24_500_clara, file.path(rdsdir, sini))
        mmssini_MW24_500_clara <- massbicm(ssini_MW24_500_clara, file.path(rdsdir, sini))
        nmssini_MW24_500_mclVVV <- massbic(ssini_MW24_500_mclVVV, file.path(rdsdir, sini))
        mmssini_MW24_500_mclVVV <- massbicm(ssini_MW24_500_mclVVV, file.path(rdsdir, sini))
nmssini_MW28 <- massbic(ssini_MW28, file.path(rdsdir, sini))
mmssini_MW28 <- massbicm(ssini_MW28, file.path(rdsdir, sini))
    nmssini_MW28_500 <- massbic(ssini_MW28_500, file.path(rdsdir, sini))
    mmssini_MW28_500 <- massbicm(ssini_MW28_500, file.path(rdsdir, sini))
        nmssini_MW28_500_clara <- massbic(ssini_MW28_500_clara, file.path(rdsdir, sini))
        mmssini_MW28_500_clara <- massbicm(ssini_MW28_500_clara, file.path(rdsdir, sini))
        nmssini_MW28_500_mclVVV <- massbic(ssini_MW28_500_mclVVV, file.path(rdsdir, sini))
        mmssini_MW28_500_mclVVV <- massbicm(ssini_MW28_500_mclVVV, file.path(rdsdir, sini))
    nmssini_MW28_1000 <- massbic(ssini_MW28_1000, file.path(rdsdir, sini))
    mmssini_MW28_1000 <- massbicm(ssini_MW28_1000, file.path(rdsdir, sini))
        nmssini_MW28_1000_clara <- massbic(ssini_MW28_1000_clara, file.path(rdsdir, sini))
        mmssini_MW28_1000_clara <- massbicm(ssini_MW28_1000_clara, file.path(rdsdir, sini))
        nmssini_MW28_1000_mclVVV <- massbic(ssini_MW28_1000_mclVVV, file.path(rdsdir, sini))
        mmssini_MW28_1000_mclVVV <- massbicm(ssini_MW28_1000_mclVVV, file.path(rdsdir, sini))
nmssini_MW29 <- massbic(ssini_MW29, file.path(rdsdir, sini))
mmssini_MW29 <- massbicm(ssini_MW29, file.path(rdsdir, sini))
    nmssini_MW29_500 <- massbic(ssini_MW29_500, file.path(rdsdir, sini))
    mmssini_MW29_500 <- massbicm(ssini_MW29_500, file.path(rdsdir, sini))
        nmssini_MW29_500_clara <- massbic(ssini_MW29_500_clara, file.path(rdsdir, sini))
        mmssini_MW29_500_clara <- massbicm(ssini_MW29_500_clara, file.path(rdsdir, sini))
        nmssini_MW29_500_mclVVV <- massbic(ssini_MW29_500_mclVVV, file.path(rdsdir, sini))
        mmssini_MW29_500_mclVVV <- massbicm(ssini_MW29_500_mclVVV, file.path(rdsdir, sini))
    nmssini_MW29_1000 <- massbic(ssini_MW29_1000, file.path(rdsdir, sini))
    mmssini_MW29_1000 <- massbicm(ssini_MW29_1000, file.path(rdsdir, sini))
        nmssini_MW29_1000_clara <- massbic(ssini_MW29_1000_clara, file.path(rdsdir, sini))
        mmssini_MW29_1000_clara <- massbicm(ssini_MW29_1000_clara, file.path(rdsdir, sini))
        nmssini_MW29_1000_mclVVV <- massbic(ssini_MW29_1000_mclVVV, file.path(rdsdir, sini))
        mmssini_MW29_1000_mclVVV <- massbicm(ssini_MW29_1000_mclVVV, file.path(rdsdir, sini))
nmssini_MW210 <- massbic(ssini_MW210, file.path(rdsdir, sini))
mmssini_MW210 <- massbicm(ssini_MW210, file.path(rdsdir, sini))
    nmssini_MW210_500 <- massbic(ssini_MW210_500, file.path(rdsdir, sini))
    mmssini_MW210_500 <- massbicm(ssini_MW210_500, file.path(rdsdir, sini))
        nmssini_MW210_500_clara <- massbic(ssini_MW210_500_clara, file.path(rdsdir, sini))
        mmssini_MW210_500_clara <- massbicm(ssini_MW210_500_clara, file.path(rdsdir, sini))
        nmssini_MW210_500_mclVVV <- massbic(ssini_MW210_500_mclVVV, file.path(rdsdir, sini))
        mmssini_MW210_500_mclVVV <- massbicm(ssini_MW210_500_mclVVV, file.path(rdsdir, sini))
    nmssini_MW210_1000 <- massbic(ssini_MW210_1000, file.path(rdsdir, sini))
    mmssini_MW210_1000 <- massbicm(ssini_MW210_1000, file.path(rdsdir, sini))
        nmssini_MW210_1000_clara <- massbic(ssini_MW210_1000_clara, file.path(rdsdir, sini))
        mmssini_MW210_1000_clara <- massbicm(ssini_MW210_1000_clara, file.path(rdsdir, sini))
        nmssini_MW210_1000_mclVVV <- massbic(ssini_MW210_1000_mclVVV, file.path(rdsdir, sini))
        mmssini_MW210_1000_mclVVV <- massbicm(ssini_MW210_1000_mclVVV, file.path(rdsdir, sini))
nmssini_MW213 <- massbic(ssini_MW213, file.path(rdsdir, sini))
mmssini_MW213 <- massbicm(ssini_MW213, file.path(rdsdir, sini))
    nmssini_MW213_500 <- massbic(ssini_MW213_500, file.path(rdsdir, sini))
    mmssini_MW213_500 <- massbicm(ssini_MW213_500, file.path(rdsdir, sini))
        nmssini_MW213_500_clara <- massbic(ssini_MW213_500_clara, file.path(rdsdir, sini))
        mmssini_MW213_500_clara <- massbicm(ssini_MW213_500_clara, file.path(rdsdir, sini))
        nmssini_MW213_500_mclVVV <- massbic(ssini_MW213_500_mclVVV, file.path(rdsdir, sini))
        mmssini_MW213_500_mclVVV <- massbicm(ssini_MW213_500_mclVVV, file.path(rdsdir, sini))
    nmssini_MW213_1000 <- massbic(ssini_MW213_1000, file.path(rdsdir, sini))
    mmssini_MW213_1000 <- massbicm(ssini_MW213_1000, file.path(rdsdir, sini))
        nmssini_MW213_1000_clara <- massbic(ssini_MW213_1000_clara, file.path(rdsdir, sini))
        mmssini_MW213_1000_clara <- massbicm(ssini_MW213_1000_clara, file.path(rdsdir, sini))
        nmssini_MW213_1000_mclVVV <- massbic(ssini_MW213_1000_mclVVV, file.path(rdsdir, sini))
        mmssini_MW213_1000_mclVVV <- massbicm(ssini_MW213_1000_mclVVV, file.path(rdsdir, sini))
nmssini_MW214 <- massbic(ssini_MW214, file.path(rdsdir, sini))
mmssini_MW214 <- massbicm(ssini_MW214, file.path(rdsdir, sini))
    nmssini_MW214_500 <- massbic(ssini_MW214_500, file.path(rdsdir, sini))
    mmssini_MW214_500 <- massbicm(ssini_MW214_500, file.path(rdsdir, sini))
        nmssini_MW214_500_clara <- massbic(ssini_MW214_500_clara, file.path(rdsdir, sini))
        mmssini_MW214_500_clara <- massbicm(ssini_MW214_500_clara, file.path(rdsdir, sini))
        nmssini_MW214_500_mclVVV <- massbic(ssini_MW214_500_mclVVV, file.path(rdsdir, sini))
        mmssini_MW214_500_mclVVV <- massbicm(ssini_MW214_500_mclVVV, file.path(rdsdir, sini))
    nmssini_MW214_1000 <- massbic(ssini_MW214_1000, file.path(rdsdir, sini))
    mmssini_MW214_1000 <- massbicm(ssini_MW214_1000, file.path(rdsdir, sini))
        nmssini_MW214_1000_clara <- massbic(ssini_MW214_1000_clara, file.path(rdsdir, sini))
        mmssini_MW214_1000_clara <- massbicm(ssini_MW214_1000_clara, file.path(rdsdir, sini))
        nmssini_MW214_1000_mclVVV <- massbic(ssini_MW214_1000_mclVVV, file.path(rdsdir, sini))
        mmssini_MW214_1000_mclVVV <- massbicm(ssini_MW214_1000_mclVVV, file.path(rdsdir, sini))

cat("sini done \n")

cat("starting plots \n")

pdf(file = file.path(savdir, "small_size.pdf"))
massplot(nmsssiz); dev.off()
pdf(file = file.path(savdir, "mcl_small_size.pdf"))
massplot(mmsssiz); dev.off()
pdf(file = file.path(savdir, "small_seed.pdf"))
massplot(nmsssee); dev.off()
pdf(file = file.path(savdir, "mcl_small_seed.pdf"))
massplot(mmsssee); dev.off()
pdf(file = file.path(savdir, "nm2.pdf"))
massplot(nmsnmm2); dev.off()
pdf(file = file.path(savdir, "mcl_nm2.pdf"))
massplot(mmsnmm2); dev.off()
pdf(file = file.path(savdir, "small_smi.pdf"))
massplot(nmsssmi); dev.off()
pdf(file = file.path(savdir, "mcl_small_smi.pdf"))
massplot(mmsssmi); dev.off()
pdf(file = file.path(savdir, "small_ini.pdf"))
massplot(nmssini); dev.off()
pdf(file = file.path(savdir, "mcl_small_ini.pdf"))
massplot(mmssini); dev.off()

pdf(file = file.path(savdir, "small_size_0500.pdf"))
massplot(nmsssiz_0500); dev.off()
pdf(file = file.path(savdir, "mcl_small_size_0500.pdf"))
massplot(mmsssiz_0500); dev.off()
pdf(file = file.path(savdir, "small_size_0600.pdf"))
massplot(nmsssiz_0600); dev.off()
pdf(file = file.path(savdir, "mcl_small_size_0600.pdf"))
massplot(mmsssiz_0600); dev.off()
pdf(file = file.path(savdir, "small_size_0700.pdf"))
massplot(nmsssiz_0700); dev.off()
pdf(file = file.path(savdir, "mcl_small_size_0700.pdf"))
massplot(mmsssiz_0700); dev.off()
pdf(file = file.path(savdir, "small_size_0800.pdf"))
massplot(nmsssiz_0800); dev.off()
pdf(file = file.path(savdir, "mcl_small_size_0800.pdf"))
massplot(mmsssiz_0800); dev.off()
pdf(file = file.path(savdir, "small_size_0900.pdf"))
massplot(nmsssiz_0900); dev.off()
pdf(file = file.path(savdir, "mcl_small_size_0900.pdf"))
massplot(mmsssiz_0900); dev.off()
pdf(file = file.path(savdir, "small_size_1000.pdf"))
massplot(nmsssiz_1000); dev.off()
pdf(file = file.path(savdir, "mcl_small_size_1000.pdf"))
massplot(mmsssiz_1000); dev.off()
pdf(file = file.path(savdir, "small_size_1100.pdf"))
massplot(nmsssiz_1100); dev.off()
pdf(file = file.path(savdir, "mcl_small_size_1100.pdf"))
massplot(mmsssiz_1100); dev.off()
pdf(file = file.path(savdir, "small_size_1200.pdf"))
massplot(nmsssiz_1200); dev.off()
pdf(file = file.path(savdir, "mcl_small_size_1200.pdf"))
massplot(mmsssiz_1200); dev.off()
pdf(file = file.path(savdir, "small_size_1300.pdf"))
massplot(nmsssiz_1300); dev.off()
pdf(file = file.path(savdir, "mcl_small_size_1300.pdf"))
massplot(mmsssiz_1300); dev.off()
pdf(file = file.path(savdir, "small_size_1400.pdf"))
massplot(nmsssiz_1400); dev.off()
pdf(file = file.path(savdir, "mcl_small_size_1400.pdf"))
massplot(mmsssiz_1400); dev.off()
pdf(file = file.path(savdir, "small_size_1500.pdf"))
massplot(nmsssiz_1500); dev.off()
pdf(file = file.path(savdir, "mcl_small_size_1500.pdf"))
massplot(mmsssiz_1500); dev.off()


pdf(file = file.path(savdir, "small_seed_01.pdf"))
massplot(nmsssee_seed01); dev.off()
pdf(file = file.path(savdir, "mcl_small_seed_01.pdf"))
massplot(mmsssee_seed01); dev.off()
pdf(file = file.path(savdir, "small_seed_02.pdf"))
massplot(nmsssee_seed02); dev.off()
pdf(file = file.path(savdir, "mcl_small_seed_02.pdf"))
massplot(mmsssee_seed02); dev.off()
pdf(file = file.path(savdir, "small_seed_03.pdf"))
massplot(nmsssee_seed03); dev.off()
pdf(file = file.path(savdir, "mcl_small_seed_03.pdf"))
massplot(mmsssee_seed03); dev.off()
pdf(file = file.path(savdir, "small_seed_04.pdf"))
massplot(nmsssee_seed04); dev.off()
pdf(file = file.path(savdir, "mcl_small_seed_04.pdf"))
massplot(mmsssee_seed04); dev.off()
pdf(file = file.path(savdir, "small_seed_05.pdf"))
massplot(nmsssee_seed05); dev.off()
pdf(file = file.path(savdir, "mcl_small_seed_05.pdf"))
massplot(mmsssee_seed05); dev.off()
pdf(file = file.path(savdir, "small_seed_06.pdf"))
massplot(nmsssee_seed06); dev.off()
pdf(file = file.path(savdir, "mcl_small_seed_06.pdf"))
massplot(mmsssee_seed06); dev.off()
pdf(file = file.path(savdir, "small_seed_07.pdf"))
massplot(nmsssee_seed07); dev.off()
pdf(file = file.path(savdir, "mcl_small_seed_07.pdf"))
massplot(mmsssee_seed07); dev.off()
pdf(file = file.path(savdir, "small_seed_08.pdf"))
massplot(nmsssee_seed08); dev.off()
pdf(file = file.path(savdir, "mcl_small_seed_08.pdf"))
massplot(mmsssee_seed08); dev.off()
pdf(file = file.path(savdir, "small_seed_09.pdf"))
massplot(nmsssee_seed09); dev.off()
pdf(file = file.path(savdir, "mcl_small_seed_09.pdf"))
massplot(mmsssee_seed09); dev.off()
pdf(file = file.path(savdir, "small_seed_10.pdf"))
massplot(nmsssee_seed10); dev.off()
pdf(file = file.path(savdir, "mcl_small_seed_10.pdf"))
massplot(mmsssee_seed10); dev.off()


pdf(file = file.path(savdir, "small_smi.pdf"))
massplot(nmsssmi); dev.off()
pdf(file = file.path(savdir, "mcl_small_smi.pdf"))
massplot(mmsssmi); dev.off()


pdf(file = file.path(savdir, "nm2_MW210.pdf"))
massplot(nmsnmm2_MW210); dev.off()
pdf(file = file.path(savdir, "mcl_nm2_MW210.pdf"))
massplot(mmsnmm2_MW210); dev.off()
pdf(file = file.path(savdir, "nm2_MW210_400.pdf"))
massplot(nmsnmm2_MW210_400); dev.off()
pdf(file = file.path(savdir, "mcl_nm2_MW210_400.pdf"))
massplot(mmsnmm2_MW210_400); dev.off()
pdf(file = file.path(savdir, "nm2_MW210_400_clara.pdf"))
massplot(nmsnmm2_MW210_400_clara ); dev.off()
pdf(file = file.path(savdir, "mcl_nm2_MW210_400_clara.pdf"))
massplot(mmsnmm2_MW210_400_clara ); dev.off()
pdf(file = file.path(savdir, "nm2_MW210_400_mclVVV.pdf"))
massplot(nmsnmm2_MW210_400_mclVVV); dev.off()
pdf(file = file.path(savdir, "mcl_nm2_MW210_400mclVVV.pdf"))
massplot(mmsnmm2_MW210_400_mclVVV); dev.off()
pdf(file = file.path(savdir, "nm2_MW210_500.pdf"))
massplot(nmsnmm2_MW210_500); dev.off()
pdf(file = file.path(savdir, "mcl_nm2_MW210_500.pdf"))
massplot(mmsnmm2_MW210_500); dev.off()
pdf(file = file.path(savdir, "nm2_MW210_500_clara.pdf"))
massplot(nmsnmm2_MW210_500_clara ); dev.off()
pdf(file = file.path(savdir, "mcl_nm2_MW210_500_clara.pdf"))
massplot(mmsnmm2_MW210_500_clara ); dev.off()
pdf(file = file.path(savdir, "nm2_MW210_500_mclVVV.pdf"))
massplot(nmsnmm2_MW210_500_mclVVV); dev.off()
pdf(file = file.path(savdir, "mcl_nm2_MW210_500_mclVVV.pdf"))
massplot(mmsnmm2_MW210_500_mclVVV); dev.off()
pdf(file = file.path(savdir, "nm2_MW21.pdf"))
massplot(nmsnmm2_MW21); dev.off()
pdf(file = file.path(savdir, "mcl_nm2_MW21.pdf"))
massplot(mmsnmm2_MW21); dev.off()
pdf(file = file.path(savdir, "nm2_MW21_400.pdf"))
massplot(nmsnmm2_MW21_400); dev.off()
pdf(file = file.path(savdir, "mcl_nm2_MW_400.pdf"))
massplot(mmsnmm2_MW21_400); dev.off()
pdf(file = file.path(savdir, "nm2_MW_400_clara.pdf"))
massplot(nmsnmm2_MW21_400_clara); dev.off()
pdf(file = file.path(savdir, "mcl_nm2_MW_400_clara.pdf"))
massplot(mmsnmm2_MW21_400_clara); dev.off()
pdf(file = file.path(savdir, "nm2_MW21_400_mclVVV.pdf"))
massplot(nmsnmm2_MW21_400_mclVVV ); dev.off()
pdf(file = file.path(savdir, "mcl_nm2_MW21_400_mclVVV.pdf"))
massplot(mmsnmm2_MW21_400_mclVVV ); dev.off()
pdf(file = file.path(savdir, "nm2_MW21_500.pdf"))
massplot(nmsnmm2_MW21_500); dev.off()
pdf(file = file.path(savdir, "mcl_nm2_MW21_500.pdf"))
massplot(mmsnmm2_MW21_500); dev.off()
pdf(file = file.path(savdir, "nm2_MW21_500_clara.pdf"))
massplot(nmsnmm2_MW21_500_clara); dev.off()
pdf(file = file.path(savdir, "mcl_nm2_MW21_500_clara.pdf"))
massplot(mmsnmm2_MW21_500_clara); dev.off()
pdf(file = file.path(savdir, "nm2_MW21_500_mclVVV.pdf"))
massplot(nmsnmm2_MW21_500_mclVVV ); dev.off()
pdf(file = file.path(savdir, "mcl_nm2_MW21_500_mclVVV.pdf"))
massplot(mmsnmm2_MW21_500_mclVVV ); dev.off()


pdf(file = file.path(savdir, "small_init_MW24.pdf"))
massplot(nmssini_MW24); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW24.pdf"))
massplot(mmssini_MW24); dev.off()
pdf(file = file.path(savdir, "small_init_MW24_1000.pdf"))
massplot(nmssini_MW24_1000); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW24_1000.pdf"))
massplot(mmssini_MW24_1000); dev.off()
pdf(file = file.path(savdir, "small_init_MW24_1000_clara.pdf"))
massplot(nmssini_MW24_1000_clara); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW24_1000_clara.pdf"))
massplot(mmssini_MW24_1000_clara); dev.off()
pdf(file = file.path(savdir, "small_init_MW24_1000_mclVVV.pdf"))
massplot(nmssini_MW24_1000_mclVVV ); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW24_1000_mclVVV.pdf"))
massplot(mmssini_MW24_1000_mclVVV ); dev.off()
pdf(file = file.path(savdir, "small_init_MW24_500.pdf"))
massplot(nmssini_MW24_500); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW24_500.pdf"))
massplot(mmssini_MW24_500); dev.off()
pdf(file = file.path(savdir, "small_init_MW24_500_clara.pdf"))
massplot(nmssini_MW24_500_clara); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW24_500_clara.pdf"))
massplot(mmssini_MW24_500_clara); dev.off()
pdf(file = file.path(savdir, "small_init_MW24_500_mclVVV.pdf"))
massplot(nmssini_MW24_500_mclVVV); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW24_500_mclVVV.pdf"))
massplot(mmssini_MW24_500_mclVVV); dev.off()
pdf(file = file.path(savdir, "small_init_MW28.pdf"))
massplot(nmssini_MW28); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW28.pdf"))
massplot(mmssini_MW28); dev.off()
pdf(file = file.path(savdir, "small_init_MW28_500.pdf"))
massplot(nmssini_MW28_500); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW28_500.pdf"))
massplot(mmssini_MW28_500); dev.off()
pdf(file = file.path(savdir, "small_init_MW28_500_clara.pdf"))
massplot(nmssini_MW28_500_clara); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW28_500_clara.pdf"))
massplot(mmssini_MW28_500_clara); dev.off()
pdf(file = file.path(savdir, "small_init_MW28_500_mclVVV.pdf"))
massplot(nmssini_MW28_500_mclVVV); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW28_500_mclVVV.pdf"))
massplot(mmssini_MW28_500_mclVVV); dev.off()
pdf(file = file.path(savdir, "small_init_MW28_1000.pdf"))
massplot(nmssini_MW28_1000); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW28_1000.pdf"))
massplot(mmssini_MW28_1000); dev.off()
pdf(file = file.path(savdir, "small_init_MW28_1000_clara.pdf"))
massplot(nmssini_MW28_1000_clara); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW28_1000_clara.pdf"))
massplot(mmssini_MW28_1000_clara); dev.off()
pdf(file = file.path(savdir, "small_init_MW28_1000_mclVVV.pdf"))
massplot(nmssini_MW28_1000_mclVVV ); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW28_1000_mclVVV.pdf"))
massplot(mmssini_MW28_1000_mclVVV ); dev.off()
pdf(file = file.path(savdir, "small_init_MW29.pdf"))
massplot(nmssini_MW29); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW29.pdf"))
massplot(mmssini_MW29); dev.off()
pdf(file = file.path(savdir, "small_init_MW29_500.pdf"))
massplot(nmssini_MW29_500); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW29_500.pdf"))
massplot(mmssini_MW29_500); dev.off()
pdf(file = file.path(savdir, "small_init_MW29_500_clara.pdf"))
massplot(nmssini_MW29_500_clara); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW29_500_clara.pdf"))
massplot(mmssini_MW29_500_clara); dev.off()
pdf(file = file.path(savdir, "small_init_MW29_500_mclVVV.pdf"))
massplot(nmssini_MW29_500_mclVVV); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW29_500_mclVVV.pdf"))
massplot(mmssini_MW29_500_mclVVV); dev.off()
pdf(file = file.path(savdir, "small_init_MW29_1000.pdf"))
massplot(nmssini_MW29_1000); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW29_1000.pdf"))
massplot(mmssini_MW29_1000); dev.off()
pdf(file = file.path(savdir, "small_init_MW29_1000_clara.pdf"))
massplot(nmssini_MW29_1000_clara); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW29_1000_clara.pdf"))
massplot(mmssini_MW29_1000_clara); dev.off()
pdf(file = file.path(savdir, "small_init_MW29_1000_mclVVV.pdf"))
massplot(nmssini_MW29_1000_mclVVV ); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW29_1000_mclVVV.pdf"))
massplot(mmssini_MW29_1000_mclVVV ); dev.off()
pdf(file = file.path(savdir, "small_init_MW210.pdf"))
massplot(nmssini_MW210); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW210.pdf"))
massplot(mmssini_MW210); dev.off()
pdf(file = file.path(savdir, "small_init_MW210_500.pdf"))
massplot(nmssini_MW210_500); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW210_500.pdf"))
massplot(mmssini_MW210_500); dev.off()
pdf(file = file.path(savdir, "small_init_MW210_500_clara.pdf"))
massplot(nmssini_MW210_500_clara); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW210_500_clara.pdf"))
massplot(mmssini_MW210_500_clara); dev.off()
pdf(file = file.path(savdir, "small_init_MW210_500_mclVVV.pdf"))
massplot(nmssini_MW210_500_mclVVV ); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW210_500_mclVVV.pdf"))
massplot(mmssini_MW210_500_mclVVV ); dev.off()
pdf(file = file.path(savdir, "small_init_MW210_1000.pdf"))
massplot(nmssini_MW210_1000); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW210_1000.pdf"))
massplot(mmssini_MW210_1000); dev.off()
pdf(file = file.path(savdir, "small_init_MW210_1000_clara.pdf"))
massplot(nmssini_MW210_1000_clara ); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW210_1000_clara.pdf"))
massplot(mmssini_MW210_1000_clara ); dev.off()
pdf(file = file.path(savdir, "small_init_MW210_1000_mclVVV.pdf"))
massplot(nmssini_MW210_1000_mclVVV); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW210_1000_mclVVV.pdf"))
massplot(mmssini_MW210_1000_mclVVV); dev.off()
pdf(file = file.path(savdir, "small_init_MW213.pdf"))
massplot(nmssini_MW213); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW213.pdf"))
massplot(mmssini_MW213); dev.off()
pdf(file = file.path(savdir, "small_init_MW213_500.pdf"))
massplot(nmssini_MW213_500); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW213_500.pdf"))
massplot(mmssini_MW213_500); dev.off()
pdf(file = file.path(savdir, "small_init_MW213_500_clara.pdf"))
massplot(nmssini_MW213_500_clara); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW213_500_clara.pdf"))
massplot(mmssini_MW213_500_clara); dev.off()
pdf(file = file.path(savdir, "small_init_MW213_500_mclVVV.pdf"))
massplot(nmssini_MW213_500_mclVVV ); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW213_500_mclVVV.pdf"))
massplot(mmssini_MW213_500_mclVVV ); dev.off()
pdf(file = file.path(savdir, "small_init_MW213_1000.pdf"))
massplot(nmssini_MW213_1000); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW213_1000.pdf"))
massplot(mmssini_MW213_1000); dev.off()
pdf(file = file.path(savdir, "small_init_MW213_1000_clara.pdf"))
massplot(nmssini_MW213_1000_clara ); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW213_1000_clara.pdf"))
massplot(mmssini_MW213_1000_clara ); dev.off()
pdf(file = file.path(savdir, "small_init_MW213_1000_mclVVV.pdf"))
massplot(nmssini_MW213_1000_mclVVV); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW213_1000_mclVVV.pdf"))
massplot(mmssini_MW213_1000_mclVVV); dev.off()
pdf(file = file.path(savdir, "small_init_MW214.pdf"))
massplot(nmssini_MW214); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW214.pdf"))
massplot(mmssini_MW214); dev.off()
pdf(file = file.path(savdir, "small_init_MW214_500.pdf"))
massplot(nmssini_MW214_500); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW214_500.pdf"))
massplot(mmssini_MW214_500); dev.off()
pdf(file = file.path(savdir, "small_init_MW214_500_clara.pdf"))
massplot(nmssini_MW214_500_clara); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW214_500_clara.pdf"))
massplot(mmssini_MW214_500_clara); dev.off()
pdf(file = file.path(savdir, "small_init_MW214_500_mclVVV.pdf"))
massplot(nmssini_MW214_500_mclVVV ); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW214_500_mclVVV.pdf"))
massplot(mmssini_MW214_500_mclVVV ); dev.off()
pdf(file = file.path(savdir, "small_init_MW214_1000.pdf"))
massplot(nmssini_MW214_1000); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW214_1000.pdf"))
massplot(mmssini_MW214_1000); dev.off()
pdf(file = file.path(savdir, "small_init_MW214_1000_clara.pdf"))
massplot(nmssini_MW214_1000_clara ); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW214_1000_clara.pdf"))
massplot(mmssini_MW214_1000_clara ); dev.off()
pdf(file = file.path(savdir, "small_init_MW214_1000_mclVVV.pdf"))
massplot(nmssini_MW214_1000_mclVVV); dev.off()
pdf(file = file.path(savdir, "mcl_small_init_MW214_1000_mclVVV.pdf"))
massplot(mmssini_MW214_1000_mclVVV); dev.off()


cat("plots done \n")
