devtools::load_all("BachelorArbeit/norMmix")

rdat <- list.files("~/BachelorArbeit/norMmix/data", pattern="mw.*")


for (i in rdat) {
    print(paste("work on data set:",i))
    x <- readRDS(paste0("~/BachelorArbeit/norMmix/data/", i))
    ans1 <- fit.norMmix(x, k=1:3, models=1:2, trafo="clr1", ini="cla", maxiter=300)
    saveRDS(ans1, file=paste0("~/BScThesis/4MM/","clr1.cla.k17m110.",i))

    ans2 <- fit.norMmix(x, k=1:7, models=1:10, trafo="logit", ini="cla", maxiter=300)
    saveRDS(ans2, file=paste0("~/BScThesis/4MM/","logit.cla.k17m110.",i))

    ans3 <- fit.norMmix(x, k=1:7, models=1:10, trafo="clr1", ini="mcl", maxiter=300)
    saveRDS(ans3, file=paste0("~/BScThesis/4MM/","clr1.mcl.k17m110.",i))

    ans2 <- fit.norMmix(x, k=1:7, models=1:10, trafo="logit", ini="mcl", maxiter=300)
    saveRDS(ans2, file=paste0("~/BScThesis/4MM/","logit.mcl.k17m110.",i))
}

