\name{compplot}
\alias{compplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    composition plot
%%  ~~function to do ... ~~
}
\description{
    Creates a \code{\link{massplot}} like plot. Takes two \code{\link{massbic}}
    arrays and overlays their plots.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compplot(f, g, main = "unnamed")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
      Array as from \code{\link{massbic}} or \code{\link{massbicm}}
}
  \item{g}{
      Array as from \code{\link{massbic}} or \code{\link{massbicm}}
}
  \item{main}{
      Character string to be used as title of the plot.
}
}
\details{
    The intended use for this function is to run \code{\link{massbic}} and 
    \code{\link{massbicm}} with the same arguments, \code{string} and \code{DIR}
    and feeding the results into \code{compplot}
}
\value{
    No return value. Side effect is the generation of a plot.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
    Nicolas Trutmann
}
\note{
    While possible to be used as standalone, the function \code{\link{epfl}} 
    allows for more complex generation of PDFs and is the intended use for 
    simulations.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{massplot}},
    \code{\link{epfl}},
    \code{\link{massbic}},
    \code{\link{massbicm}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (f, g, main = "unnamed") 
{
    ylim <- extendrange(c(f, g))
    adj <- 0.4
    op <- sfsmisc::mult.fig(mfrow = c(2, 5), main = main, mar = 0.1 + 
        c(2, 4, 4, 1))
    models <- dimnames(f)$models
    for (i in 1:10) {
        matplot(f[, i, ], lty = 1, col = adjustcolor(rainbow(20)[i], 
            adj), main = models[i], type = "l", ylim = ylim)
        matplot(g[, i, ], lty = 1, col = adjustcolor(rainbow(20)[i + 
            10], adj), main = models[i], type = "l", ylim = ylim, 
            add = TRUE)
    }
    par(op$old.par)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
