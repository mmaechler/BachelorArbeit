\name{npar.fittednorMmix}
\alias{npar.fittednorMmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Extract degrees of freedom from \code{fittednorMmix}
}
\description{
    Extracts a matrix of equal shape as \code{object}'s BIC values.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
npar.fittednorMmix(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
      object of class \code{fittednorMmix}
%%     ~~Describe \code{obj} here~~
    }
  \item{...}{
    }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
    Nicolas Trutmann
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (obj) 
{
    stopifnot(inherits(obj, "fittednorMmix"))
    k <- obj$k
    p <- obj$p
    models <- obj$models
    val <- matrix(0, length(k), length(models))
    rownames(val) <- k
    colnames(val) <- models
    for (i in seq_along(k)) {
        for (j in seq_along(models)) {
            val[i, j] <- npar(k[i], p, models[j])
        }
    }
    val
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
