\name{compplot}
\alias{compplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    composition plot
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compplot(f, g, main = "unnamed")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
%%     ~~Describe \code{f} here~~
}
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
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
