####
## Auswertungsfunktionen für Ergebnisse von ada-simulationen

adabic <- function(string, na.rm=FALSE) {
    arr <- array(0, c(7,10,50))
    val <- list()

    for (i in 1:50) {
        dir1 <- "~/ethz/BA/Rscripts/nm2"
        nm <- readRDS(file.path(dir1,paste0(string,as.character(i),".rds")))
        arr[,,i] <- BIC(nm$fit)[[1]]
        val[[i]] <- BIC(nm$fit)[[2]]
    }

    mus <- apply(arr,c(1,2), function(j) mean(j, na.rm=na.rm))
    sds <- apply(arr,c(1,2), sd)

    list(mu=mus, sd=sds, best=val)
}


adabest <- function(string) {
    ret <- rle(sort(unlist(lapply(adabic(string1)$best, function(j)paste(j[[1]],j[[2]])))))
    cbind(ret$values, ret$lengths)
}


adall <- function(string, model) {
    rettrue <- numeric(50)
    retfit <- numeric(50)
    truenm <- get(model, "package:norMmix")

    for (i in 1:50) {
        dir1 <- "~/ethz/BA/Rscripts/nm2"
        nm <- readRDS(file.path(dir1,paste0(string,as.character(i),".rds")))
        val <- BIC(nm$fit)[[2]]
        fitnm <- nm$fit$nMm[val[1],val[2]][[1]]$norMmix
        retfit[i] <- tryCatch(sllnorMmix(nm$fit$x, fitnm), error = function(e) NA)
        rettrue[i] <- sllnorMmix(nm$fit$x, truenm)
    }

    list(retfit,rettrue,retfit/rettrue)
}



massbic <- function(string, DIR) {
    nm1 <- readRDS(file=file.path(DIR,string[1]))
    aa <- dim(nm1$fit$nMm)
    val <- array(0, c(aa[1],aa[2],length(string)))

    for (i in 1:length(string)) {
        nm <- readRDS(file=file.path(DIR,string[i]))
        val[,,i] <- BIC(nm$fit)[[1]]
    }

    mus <- apply(val,3, mean)
    sds <- apply(val,3, sd)
    list(mean=mus, sd=sds, v=val)
}



massplot <- function(f) {
    op <- par(mfrow=c(4,5))
    for (i in 1:10) {
        matplot(f[,i,], lty=1, col=rainbow(10)[i], main=models[i], type="l")
    }
    for (i in 1:10) {
        boxplot(t(f[,i,]), lty=1, col=rainbow(10)[i], main=models[i], type="l")
    }
    par(op)
}
