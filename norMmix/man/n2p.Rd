% documentation in R/param.R
\name{nv2p}
\title{Wrapper function for nMm objs}
%\alias{n2p}
\alias{nc2p}
\description{
  \code{nc2p} returns same as \code{\link{nMm2par}}, using \code{obj$model} as default.
}
\usage{
nc2p(obj)
}
\arguments{
  \item{obj}{an \code{"\link{norMmix}"} object.}
}
\value{
  the same as \code{\link{nMm2par}(..)}, see there.
}
\examples{
str(MW213)
nc2p(MW213)
}
\keyword{misc}
