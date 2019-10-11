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
