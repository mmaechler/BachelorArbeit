\name{massbic}
\alias{massbic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    extract BIC from .rds files
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
massbic(string, DIR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{string}{
%%     ~~Describe \code{string} here~~
}
  \item{DIR}{
%%     ~~Describe \code{DIR} here~~
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
function (string, DIR) 
{
    nm1 <- readRDS(file = file.path(DIR, string[1]))
    cl <- nm1$fit$k
    mo <- nm1$fit$models
    val <- array(0, lengths(list(cl, mo, string)))
    dims <- vector(mode = "integer", length = length(string))
    for (i in 1:length(string)) {
        nm <- readRDS(file = file.path(DIR, string[i]))
        val[, , i] <- BIC(nm$fit)[[1]]
        dims[i] <- nm$fit$p
    }
    dimnames(val) <- list(clusters = cl, models = mo, simulation = string)
    attr(val, "dims") <- dims
    val
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
