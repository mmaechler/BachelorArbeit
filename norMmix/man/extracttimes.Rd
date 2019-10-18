\name{extracttimes}
\alias{extracttimes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Extract system time from \code{fittednorMmix} object
}
\description{
    extracts array of \code{\link{system.time}} values from a 
    \code{fittednorMmix} object.
}
\usage{
extracttimes(object, ...)
}
\arguments{
  \item{object}{
      a \code{fittednorMmix} object.
}
  \item{\dots}{
      currently unused %% TODO: need this??
}
}
\details{
    a
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    Array of size components x models x 5(return values of system.time)
    with dimnames: \code{k}, \code{models} and \code{proc_time}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
    Nicolas Trutmann
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{system.time}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, ...) 
{
    stopifnot(inherits(object, "fittednorMmix"))
    ti <- unlist(object$nMmtime)
    na <- names(ti)[1:5]
    co <- object$k
    mo <- object$models
    ti <- c(matrix(ti, ncol = 5, byrow = TRUE))
    r <- array(ti, lengths(list(co, mo, na)))
    dimnames(r) <- list(k = co, models = mo, proc_time = na)
    class(r) <- "fittednorMmix_time"
    r
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
