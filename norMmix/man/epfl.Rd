\name{epfl}
\alias{epfl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    evaluate and plot from file list
}
\description{
    From a list of character vectors, will apply \code{\link{massbic}}, 
    \code{\link{massbicm}} followed by \code{\link{massplot}} and 
    \code{\link{compplot}}

}
\usage{
epfl(files, savdir, subt = 11)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{ 
      Expected to be of type \code{\link{list}}, containing character vectors. 
      These are assumed to be RDS filenames produced from \code{\link{saveRDS}}.
}
  \item{savdir}{
      String specifying directory in which \code{\link{files}} are located.
}
  \item{subt}{
      Number of characters to be subtracted from string in \code{\link{files}}.
      The default is intended for format ending in "seed=%0.2d.rds".
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
Nicolas Trutmann
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
function (files, savdir, subt = 11) 
{
    stopifnot(is.list(files), dir.exists(savdir))
    setwd(savdir)
    for (fi in files) {
        if (length(fi) == 0) {
            next
        }
        main <- substring(fi[1], 1, nchar(fi[1]) - subt)
        f <- massbic(fi, savdir)
        g <- massbicm(fi, savdir)
        pdf(file = paste0(main, ".pdf"))
        massplot(f, main = main)
        dev.off()
        pdf(file = paste0(main, "_mcl.pdf"))
        massplot(g, main = paste0(main, "_mcl"))
        dev.off()
        pdf(file = paste0(main, "_comp.pdf"))
        compplot(f, g, main = paste0(main, "_comp"))
        dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
