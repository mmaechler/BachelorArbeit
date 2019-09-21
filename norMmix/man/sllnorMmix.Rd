% Code in ../R/llnorMmix.R
\name{sllnorMmix}
\alias{sllnorMmix}
\title{Simple wrapper for log-likelihood function}
\usage{
llnorMmix(x, obj)
}
\arguments{
    \item{x}{Data vector; rows observations, columns variables}
    \item{obj}{Object of S3 class \code{\link{norMmix}}}
}
\description{
    \code{llnorMmix} returns scalar. Log-likelihood of data given a normal mixture
}
\details{
    Calculates log-likelihood of a dataset, x, given a normal mixture model.
    Intended as a simplified wrapper for \code{\link{llnorMmix}}.
    Removes functionality in favor of ease of use.
}
\examples{
set.seed(2019); x <- rnorMmix(400, MW27)
sllnorMmix(x, MW27)
# [1] -1986.315
}

