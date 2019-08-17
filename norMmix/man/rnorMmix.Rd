\name{rnorMmix}
\alias{rnorMmix}
\title{Random Sample from Multivariate Normal Mixture Distribution}
\description{
  Draw \code{n} (p-dimensional) observations randomly from the multivariate normal
  mixture distribution specified by \code{obj}.
}
\usage{
rnorMmix(n, obj, index = FALSE, permute = TRUE)
}
\arguments{
  \item{n}{sample size, non-negative.}
  \item{obj}{a \code{"\link{norMmix}"} object}
  \item{index}{store the clustering information as first column}
  \item{permute}{logical indicating if the observations should be randomly
    permuted after creation \dQuote{cluster by cluster}.}
}
%% \details{
%% }
\value{
  n p-dimensional observations, as numeric \eqn{n \times p}{n * p} matrix.
}
%% \references{
%% }
\author{Nicolas Trutmann}
\seealso{
  \code{\link{rmultinom}}
}
\examples{
}
