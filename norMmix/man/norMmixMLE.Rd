% Code in ../R/norMmixMLE.R
\name{norMmixMLE}
\title{Maximum likelihood Estimation for normal mixture models}
\alias{norMmixMLE}
\alias{ssClaraL}
\usage{
norMmixMLE(x, k, trafo = c("clr1", "logit"),
          model = c("EII","VII","EEI","VEI","EVI",
                    "VVI","EEE","VEE","EVV","VVV"),
          ini = c("clara", "mclVVV"),
          method = "BFGS", maxit = 100, trace = 2, reltol = sqrt(.Machine$double.eps),
          samples = 128, sampsize = ssClaraL, traceClara = 0,
	  \dots)

ssClaraL(n, k, p)
}
\arguments{
  \item{x}{numeric [n x p] matrix}
  \item{k}{positive number of clusters}
  \item{trafo}{a \code{\link{character}} string, transformation to be used.}
  \item{model}{a \code{\link{character}} string, the model to be assumed.}
  \item{ini}{a \code{\link{character}} string specifying the initialization step.}

  \item{method, maxit, reltol, trace, \dots}{arguments for turning the
    optimizer \code{\link{optim}(*, method=method, control=list(...))}.}

  \item{samples, sampsize, traceClara}{if \code{ini = "clara"}, arguments for
    \code{\link[cluster]{clara}()} (package \CRANpkg{cluster}).  Note that
    \code{clara}'s help page emphasizes that larger and more samples should
    be used typically.  Here, \code{sampsize} may be a number \emph{or} as
    by default a \code{\link{function}(n, k, p)} determining the size of
    the subsamples as a function of the problem dimensionalities.}

  \item{n,p}{matrix dimensions \code{\link{nrow}(x)} and \code{\link{ncol}(x)}.}
}
\description{
  \code{norMmixMLE} returns fitted nMm obj
}
\value{
  \code{norMmixMLE} returns fitted nMm obj
}
\details{
  Uses clara() and one M-step from EM-algorithm to initialize parameters
  after that uses general optimizer optim() to calculate ML.
}
\examples{
}
