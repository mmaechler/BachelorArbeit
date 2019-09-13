% Code in ../R/llnorMmix.R
\name{llnorMmix}
\alias{llnorMmix}
\title{Log-likelihood of parameter vector given data}
\usage{
llnorMmix(par., tx, k, trafo = c("clr1", "logit"),
          model = c("EII", "VII", "EEI", "VEI", "EVI",
                    "VVI", "EEE", "VEE", "EVV", "VVV"))
}
\arguments{
    \item{par.}{parameter vector}
    \item{tx}{Transposed sample matrix. Rows are variables and columns are observations}
    \item{k}{number of clusters}
    \item{trafo}{either centered log ratio or logit}
    \item{model}{assumed distribution model of normal mixture}
}
\description{
    \code{llnorMmix} returns scalar log-likelihood
}
\details{
    Calculates log-likelihood of a dataset, tx, given a normal mixture model as 
    specified by a parameter vector.
    A parameter vector can be obtained by applying \code{\link{nMm2par}} to a 
    \code{\link{norMmix}} object.
}
\examples{
}

