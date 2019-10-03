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
    cl <- nm1$fit$k
    mo <- nm1$fit$models

    val <- array(0, lengths(list(cl, mo, string)))

    for (i in 1:length(string)) {
        nm <- readRDS(file=file.path(DIR,string[i]))
        val[,,i] <- BIC(nm$fit)[[1]]
    }
    dimnames(val) <- list(clusters=cl, models=mo, simulation=string)

    val
}


massbicm <- function(string, DIR) {
    nm <- readRDS(file.path(DIR, string[1]))
    cl <- nm$fit$k
    mo <- nm$fit$models
    valm <- array(0, lengths(list(cl,mo,string)))
    for (i in 1:length(string)) {
        nm <- readRDS(file.path(DIR, string[i]))
        x <- nm$fit$x
        valm[,,i] <- Mclust(x, G=cl, modelNames=mo)$BIC
    }
    dimnames(valm) <- list(clusters=cl, models=mo, files=string)
    -valm
}


massbest <- function(f) {
    cl <- dimnames(f)[[1]]
    mo <- dimnames(f)[[2]]
    l <- dim(f)[3]
    be <- matrix(0, 0, 2)
    for (i in 1:l) {
        be <- rbind(be, which(f[,,i]==min(f[,,i]), arr.ind=TRUE))
    }
    best <- cbind(cl[be[,1]], mo[be[,2]])
    best
}


sortbest <- function(m) {
    r <- apply(m, 1, function(j) paste(j[1], j[2]))
    ss <- rle(sort(r))
    ix <- sort(ss$lengths, decreasing=TRUE, index.return=TRUE)$ix
    ret <- list(values=ss$values[ix], reps=ss$lengths[ix])
    ret
}


#massplot <- function(f) {
#    ran <- extendrange(f)
#    size <- dim(f)[3]
#    if (size>999) adj <- 0.05
#    else if (size) adj <- 0.1
#    else if (size) adj <- 0.3
#    else adj <- 0.9
#    models <- c("EII","VII","EEI","VEI","EVI",
#                "VVI","EEE","VEE","EVV","VVV")
#    op <- par(mfrow=c(4,5))
#    for (i in 1:10) {
#        matplot(f[,i,], lty=1, col=adjustcolor(rainbow(10)[i],adj), 
#                main=models[i], type="l", ylim=ran)
#    }
#    for (i in 1:10) {
#        boxplot(t(f[,i,]), lty=1, col=adjustcolor(rainbow(10)[i],0.4), 
#                main=models[i], type="l", ylim=ran)
#    }
#    par(op)
#}

massplot <- function(f, main="unnamed", p) {
    ran <- extendrange(f)
    size <- dim(f)[3]
    cl <- as.numeric(dimnames(f)$clusters)
    adj <- exp(-0.002*size) # set so at n=1000 alpha value is ~0.1
    models <- mods()
    ## FIXME: should not assume all models present.
    op <- sfsmisc::mult.fig(mfrow=c(4,5), main=main)
    for (i in 1:10) {
        matplot(f[,i,], lty=1, col=adjustcolor(rainbow(10)[i],adj), 
                main=models[i], type="l", ylim=ran)
        if (!missing(p)) {
            axis(3, at=seq_along(cl), 
                 label=parlen(cl,p,models[i]))
        }
    }
    for (i in 1:10) {
        boxplot(t(f[,i,]), lty=1, col=adjustcolor(rainbow(10)[i],0.4), 
                main=models[i], type="l", ylim=ran)
    }
    par(op$old.par)
}


compplot <- function(f, g) {
    ylim <- extendrange(c(f,g))
    adj <- 0.4
    op <- par(mfrow=c(2,5))
    for (i in 1:10) {
        matplot(f[,i,], lty=1, col=adjustcolor(rainbow(20)[i],adj), 
                main=models[i], type="l", ylim=ylim)
        matplot(g[,i,], lty=1, col=adjustcolor(rainbow(20)[i+10],adj), 
                main=models[i], type="l", ylim=ylim, add=TRUE)
    }
    par(op)
}


optrcount <- function(obj) {
    nm <- obj$fit
    cl <- nm$k
    mo <- nm$models
    valfn <- matrix(0,length(cl),length(mo))
    valgr <- matrix(0,length(cl),length(mo))
    dimnames(valfn) <- list(clusters=cl, models=mo)
    dimnames(valgr) <- list(clusters=cl, models=mo)

    for (i in cl) {
        for (j in mo) {
            valfn[i,j] <- nm$nMm[i,j][[1]]$optr$counts[1]
            valgr[i,j] <- nm$nMm[i,j][[1]]$optr$counts[2]
        }
    }

    list(fn=valfn, gr=valgr)
}

masscount <- function(string, DIR) {
    nm1 <- readRDS(file=file.path(DIR,string[1]))
    cl <- nm1$fit$k
    mo <- nm1$fit$models

    valfn <- array(0, lengths(list(cl, mo, string)))
    valgr <- array(0, lengths(list(cl, mo, string)))

    for (i in 1:length(string)) {
        nm <- readRDS(file=file.path(DIR,string[i]))
        valfn[,,i] <- optrcount(nm)$fn
        valgr[,,i] <- optrcount(nm)$fn
    }
    dimnames(valfn) <- list(clusters=cl, models=mo, simulation=string)
    dimnames(valgr) <- list(clusters=cl, models=mo, simulation=string)

    list(fn=valfn, gr=valgr)
}
