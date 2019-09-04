% Code in ../R/fit.R
\name{fit.norMmix}
\title{fit normal mixture models to data}
\usage{
fit.norMmix <- function(x, k=1:10, models=1:10, trafo=c("clr1","logit"),ll = c("nmm", "mvt"),...) 
}
\arguments{
    \item{x}{data matrix, rows are observations and columns are variables}
    \item{k}{ vector containing a subset of 1:10,  indicating number of clusters to fit}
    \item{models}{ vector as k, indicating models from list: c("EII","VII","EEI","VEI","EVI","VVI","EEE","VEE","EVV","VVV")}
    \item{trafo, ll, \dots}{ arguments controlling the norMmixMLE function, \dots are passed to norMmixMLE}
}
\description{
    \code{fit.norMmix} returns list of fitted norMmix objects
}
\value{
    \code{fit.norMmix} returns list of 5:
    \item{nMm}{list containing all fitted models}
    \item{models}{vector containing fitted model names}
    \item{n}{number of observations of x}
    \item{p}{number of variables of x}
}
\details{Essentially wrapper function for norMmixMLE, applying norMmixMLE to all k and models specified}
\examples{}
