devtools::load_all("~/ethz/BA/norMmix")


rdat <- list.files("~/ethz/BA/Rscripts", pattern="mw.*")


for (i in rdat) {
    print(paste("work on data set:",i))
    x <- readRDS(paste0("~/ethz/BA/Rscripts/", i))
    ans1 <- fit.norMmix(x, k=1:3, models=1:2, trafo="clr1", ini="cla", maxiter=300)
    saveRDS(ans1, file=paste0("~/ethz/BA/Rscripts/","clr1.cla.k17m110.",i))

    ans2 <- fit.norMmix(x, k=1:7, models=1:10, trafo="logit", ini="cla", maxiter=300)
    saveRDS(ans2, file=paste0("~/ethz/BA/Rscripts/","logit.cla.k17m110.",i))

    ans3 <- fit.norMmix(x, k=1:7, models=1:10, trafo="clr1", ini="mcl", maxiter=300)
    saveRDS(ans3, file=paste0("~/ethz/BA/Rscripts/","clr1.mcl.k17m110.",i))

    ans2 <- fit.norMmix(x, k=1:7, models=1:10, trafo="logit", ini="mcl", maxiter=300)
    saveRDS(ans2, file=paste0("~/ethz/BA/Rscripts/","logit.mcl.k17m110.",i))

}

