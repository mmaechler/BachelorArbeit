pkgname <- "norMmix"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('norMmix')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("MarronWand")
### * MarronWand

flush(stderr()); flush(stdout())

### Name: MarronWand
### Title: Marron-Wand-like Specific Multivariate Normal Mixture 'norMmix'
###   Objects
### Aliases: MarronWand MW21 MW22 MW23 MW24 MW25 MW26 MW27 MW28 MW29 MW210
###   MW211 MW212 MW213 MW214 MW31 MW32 MW33 MW34 MW51
### Keywords: datasets distribution

### ** Examples

MW210
plot(MW214)



cleanEx()
nameEx("compplot")
### * compplot

flush(stderr()); flush(stdout())

### Name: compplot
### Title: composition plot
### Aliases: compplot

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (f, g, main = "unnamed") 
{
    ylim <- extendrange(c(f, g))
    adj <- 0.4
    op <- sfsmisc::mult.fig(mfrow = c(2, 5), main = main, mar = 0.1 + 
        c(2, 4, 4, 1))
    models <- dimnames(f)$models
    for (i in 1:10) {
        matplot(f[, i, ], lty = 1, col = adjustcolor(rainbow(20)[i], 
            adj), main = models[i], type = "l", ylim = ylim)
        matplot(g[, i, ], lty = 1, col = adjustcolor(rainbow(20)[i + 
            10], adj), main = models[i], type = "l", ylim = ylim, 
            add = TRUE)
    }
    par(op$old.par)
  }



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("epfl")
### * epfl

flush(stderr()); flush(stdout())

### Name: epfl
### Title: evaluate and plot from file list
### Aliases: epfl

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (files, savdir, subt = 11) 
{
    stopifnot(is.list(files), dir.exists(savdir))
    setwd(savdir)
    for (fi in files) {
        if (length(fi) == 0) {
            next
        }
        main <- substring(fi[1], 1, nchar(fi[1]) - subt)
        f <- massbic(fi, savdir)
        g <- massbicm(fi, savdir)
        pdf(file = paste0(main, ".pdf"))
        massplot(f, main = main)
        dev.off()
        pdf(file = paste0(main, "_mcl.pdf"))
        massplot(g, main = paste0(main, "_mcl"))
        dev.off()
        pdf(file = paste0(main, "_comp.pdf"))
        compplot(f, g, main = paste0(main, "_comp"))
        dev.off()
    }
  }



cleanEx()
nameEx("fitnMm")
### * fitnMm

flush(stderr()); flush(stdout())

### Name: fitnMm
### Title: Fit Several Normal Mixture Models to a Dataset
### Aliases: fitnMm

### ** Examples

x <- rnorMmix(500, MW21)
fitnMm(x, 1:3) ## will fit all models with 1:3 clusters



cleanEx()
nameEx("ldl")
### * ldl

flush(stderr()); flush(stdout())

### Name: ldl
### Title: LDL' Cholesky Decomposition
### Aliases: ldl

### ** Examples

(L <- rbind(c(1,0,0), c(3,1,0), c(-4,5,1)))
D <- c(4,1,9)
FF <- L %*% diag(D) %*% t(L)
FF
LL <- ldl(FF)
stopifnot(all.equal(L, LL$L),
          all.equal(D, LL$D))

## rank deficient :
FF0 <- L %*% diag(c(4,0,9)) %*% t(L)
((L0 <- ldl(FF0))) #  !! now fixed with the  if(Di == 0) test
## With the "trick", it works:
stopifnot(all.equal(FF0,
                    L0$L %*% diag(L0$D) %*% t(L0$L)))
## [hint: the LDL' is no longer unique when the matrix is singular]

system.time(for(i in 1:10000) ldl(FF) ) # ~ 0.2 sec

(L <- rbind(c( 1, 0, 0, 0),
            c( 3, 1, 0, 0),
            c(-4, 5, 1, 0),
            c(-2,20,-7, 1)))
D <- c(4,1, 9, 0.5)
F4 <- L %*% diag(D) %*% t(L)
F4
L4 <- ldl(F4)
stopifnot(all.equal(L, L4$L),
          all.equal(D, L4$D))

system.time(for(i in 1:10000) ldl(F4) )


## rank deficient :
F4.0 <- L %*% diag(c(4,1,9,0)) %*% t(L)
((L0 <- ldl(F4.0)))
stopifnot(all.equal(F4.0,
                    L0$L %*% diag(L0$D) %*% t(L0$L)))

F4_0 <- L %*% diag(c(4,1,0,9)) %*% t(L)
((L0 <- ldl(F4_0)))
stopifnot(all.equal(F4_0,
                    L0$L %*% diag(L0$D) %*% t(L0$L)))


## Large
mkLDL <- function(n, rF = function(n) sample.int(n), rFD = function(n) 1+ abs(rF(n))) {
    L <- diag(nrow=n)
    L[lower.tri(L)] <- rF(n*(n-1)/2)
    list(L = L, D = rFD(n))
}

(LD <- mkLDL(17))

chkLDL <- function(n, ..., verbose=FALSE, tol = 1e-14) {
    LD <- mkLDL(n, ...)
    if(verbose) cat(sprintf("n=%3d ", n))
    n <- length(D <- LD$D)
    L <- LD$L
    M <- L %*% diag(D) %*% t(L)
    r <- ldl(M)
    stopifnot(exprs = {
        all.equal(M,
                  r$L %*% diag(r$D) %*% t(r$L), tol=tol)
        all.equal(L, r$L, tol=tol)
        all.equal(D, r$D, tol=tol)
    })
    if(verbose) cat("[ok]\n")
    invisible(list(LD = LD, M = M, ldl = r))
}

(chkLDL(7))

N <- 99 ## test  N  random cases
set.seed(101)
for(i in 1:N) {
    cat(sprintf("i=%3d, ",i))
    chkLDL(rpois(1, lambda = 20), verbose=TRUE)
}

system.time(chkLDL( 500)) # 0.62

try( ## this almost never "works":
system.time(chkLDL( 500, rF = rnorm, rFD = function(n) 10 + runif(n))) # 0.64
)

if(interactive())
   system.time(chkLDL( 600)) # 1.09
## .. then it grows quickly for (on nb-mm4)
## for n = 1000  it typically *fails*: The matrix M  is typically very ill conditioned
## does not depend much on the RNG ?

"==> much better conditioned L and hence M : "
set.seed(120)
L <- as(Matrix::tril(toeplitz(exp(-(0:999)/50))), "matrix")
dimnames(L) <- NULL
D <- 10 + runif(nrow(L))
M <- L %*% diag(D) %*% t(L)
rcond(L) # 0.010006 !
rcond(M) # 9.4956e-5
if(FALSE) # ~ 4-5 sec
   system.time(r <- ldl(M))



cleanEx()
nameEx("massbic")
### * massbic

flush(stderr()); flush(stdout())

### Name: massbic
### Title: extract BIC from .rds files
### Aliases: massbic

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (string, DIR) 
{
    nm1 <- readRDS(file = file.path(DIR, string[1]))
    cl <- nm1$fit$k
    mo <- nm1$fit$models
    val <- array(0, lengths(list(cl, mo, string)))
    dims <- vector(mode = "integer", length = length(string))
    for (i in 1:length(string)) {
        nm <- readRDS(file = file.path(DIR, string[i]))
        val[, , i] <- BIC(nm$fit)[[1]]
        dims[i] <- nm$fit$p
    }
    dimnames(val) <- list(clusters = cl, models = mo, simulation = string)
    attr(val, "dims") <- dims
    val
  }



cleanEx()
nameEx("massbicm")
### * massbicm

flush(stderr()); flush(stdout())

### Name: massbicm
### Title: do mclust along .rds files from fitnMm
### Aliases: massbicm

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (string, DIR) 
{
    nm <- readRDS(file.path(DIR, string[1]))
    cl <- nm$fit$k
    mo <- nm$fit$models
    valm <- array(0, lengths(list(cl, mo, string)))
    dims <- vector(mode = "integer", length = length(string))
    for (i in 1:length(string)) {
        nm <- readRDS(file.path(DIR, string[i]))
        x <- nm$fit$x
        valm[, , i] <- mclust::Mclust(x, G = cl, modelNames = mo)$BIC
        dims[i] <- nm$fit$p
    }
    dimnames(valm) <- list(clusters = cl, models = mo, files = string)
    attr(valm, "dims") <- dims
    -valm
  }



cleanEx()
nameEx("massplot")
### * massplot

flush(stderr()); flush(stdout())

### Name: massplot
### Title: plot from massbic
### Aliases: massplot

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (f, main = "unnamed") 
{
    ran <- extendrange(f)
    size <- dim(f)[3]
    cl <- as.numeric(dimnames(f)$clusters)
    p <- attr(f, "dims")
    adj <- exp(-0.002 * size)
    models <- mods()
    op <- sfsmisc::mult.fig(mfrow = c(4, 5), main = main, mar = 0.1 + 
        c(2, 4, 4, 1))
    for (i in 1:10) {
        if (!is.null(p)) {
            matplot(f[, i, ], lty = 1, col = adjustcolor(rainbow(10)[i], 
                adj), type = "l", ylim = ran, main = models[i])
            axis(3, at = seq_along(cl), labels = npar(cl, p[1], 
                models[i]))
        }
        else {
            matplot(f[, i, ], lty = 1, col = adjustcolor(rainbow(10)[i], 
                adj), main = models[i], type = "l", ylim = ran)
        }
    }
    for (i in 1:10) {
        boxplot(t(f[, i, ]), lty = 1, col = adjustcolor(rainbow(10)[i], 
            0.4), main = models[i], type = "l", ylim = ran)
    }
    par(op$old.par)
  }



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("n2p")
### * n2p

flush(stderr()); flush(stdout())

### Name: nv2p
### Title: Wrapper function for nMm objs
### Aliases: nc2p
### Keywords: misc

### ** Examples

str(MW213)
nc2p(MW213)



cleanEx()
nameEx("nMm2par")
### * nMm2par

flush(stderr()); flush(stdout())

### Name: nMm2par
### Title: Multivariate Normal Mixture Model to parameter for MLE
### Aliases: nMm2par

### ** Examples

A <- MW24
if(FALSE) # currently fails  __FIXME__
nMm2par(A, trafo = "clr1", model = A$model)



cleanEx()
nameEx("norMmix")
### * norMmix

flush(stderr()); flush(stdout())

### Name: norMmix
### Title: Constructor for Multivariate Normal Mixture Objects
### Aliases: norMmix
### Keywords: distributions

### ** Examples

# TODO



cleanEx()
nameEx("norMmixMLE")
### * norMmixMLE

flush(stderr()); flush(stdout())

### Name: norMmixMLE
### Title: Maximum Likelihood Estimation for Multivariate Normal Mixture
###   Models
### Aliases: norMmixMLE ssClaraL

### ** Examples

str(MW214)
set.seed(105)
x <- rnorMmix(1000, MW214)
## Fitting assuming we know the true parametric model
fm1 <- norMmixMLE(x, k = 6, model = "VII")
if(interactive()) ## Fitting "wrong" overparametrized model: typically need more iterations:
  fmW <- norMmixMLE(x, k = 7, model = "VVV", maxit = 200)# default maxit=100 is often too small



cleanEx()
nameEx("par2nMm")
### * par2nMm

flush(stderr()); flush(stdout())

### Name: par2nMm
### Title: Transform Parameter Vector to Multivariate Normal Mixture
### Aliases: par2nMm

### ** Examples

## TODO: Show to get the list, and then how to get a  norMmix() object from the list



cleanEx()
nameEx("parlen")
### * parlen

flush(stderr()); flush(stdout())

### Name: npar
### Title: Number of Free Parameters of Multivariate Normal Mixture Models
### Aliases: npar

### ** Examples

(m <- eval(formals(npar)$model)) # list of 10 models w/ differing Sigma
# A nice table for a given 'p'  and all models, all k in 1:8
sapply(m, npar, k=setNames(,1:8), p = 20)



cleanEx()
nameEx("plot.norMmix")
### * plot.norMmix

flush(stderr()); flush(stdout())

### Name: plot.norMmix
### Title: Plot Method for "norMmix" Objects
### Aliases: plot.norMmix
### Keywords: hplot

### ** Examples

plot(MW212) ## and add a finite sample realization:
points(rnorMmix(n=500, MW212))

## or:
x <- points(rnorMmix(n=500, MW212))
plot(MW212, x)



cleanEx()
nameEx("sllnorMmix")
### * sllnorMmix

flush(stderr()); flush(stdout())

### Name: sllnorMmix
### Title: Simple wrapper for Log-Likelihood Function or Multivariate
###   Normal Mixture
### Aliases: sllnorMmix

### ** Examples

set.seed(2019)
x <- rnorMmix(400, MW27)
sllnorMmix(x, MW27) # -1986.315



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
