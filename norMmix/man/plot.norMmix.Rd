\name{plot.norMmix}
\alias{plot.norMmix}
\title{Plot Method for "norMmix" Objects}
\description{
  This is the S3 method for plotting \code{"\link{norMmix}"} objects.
}
\usage{
\S3method{plot}{norMmix}(x, ...)
}
\arguments{
  \item{x}{an \R object inheriting from \code{"\link{norMmix}"}.}
  \item{\dots}{potentical further arguments to the underlying plotting functions.}
}
\value{
  \code{plot.norMmix} returns invisibly coordinates of bounding ellipses of distribution.
}
\examples{
plot(MW212) ; points(rnorMmix(n=500, MW212))
}
\keyword{hplot}
