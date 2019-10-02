tmassbic <- function(string, DIR) {

    nm1 <- readRDS(file=file.path(DIR,string[1]))
    cl <- nm1$k
    mo <- nm1$models

    val <- array(0, lengths(list(cl, mo, string)))

    for (i in 1:length(string)) {
        nm <- readRDS(file=file.path(DIR,string[i]))
        val[,,i] <- BIC(nm)[[1]]
    }
    dimnames(val) <- list(clusters=cl, models=mo, simulation=string)

    val
}


tmassbicm <- function(string, DIR) {
    nm <- readRDS(file.path(DIR, string[1]))
    cl <- nm$k
    mo <- nm$models
    valm <- array(0, lengths(list(cl,mo,string)))
    for (i in 1:length(string)) {
        nm <- readRDS(file.path(DIR, string[i]))
        x <- nm$x
        valm[,,i] <- Mclust(x, G=cl, modelNames=models)$BIC
    }
    dimnames(valm) <- list(clusters=cl, models=mo, files=string)
    -valm
}
