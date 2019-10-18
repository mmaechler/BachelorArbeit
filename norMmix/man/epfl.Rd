\name{epfl}
\alias{epfl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    evaluate and plot from file list
}
\description{
    From a list of character vectors, will apply \code{\link{massbic}}, 
    \code{\link{massbicm}} followed by \code{\link{massplot}} and 
    \code{\link{compplot}}.
    Saves result of \code{\link{massbicm}}

}
\usage{
epfl(files, savdir, subt = 11)
}
\arguments{
  \item{files}{ 
      Expected to be of type \code{\link{list}}, containing character vectors. 
      These are assumed to be RDS filenames produced from \code{\link{saveRDS}}.
      They are assumed to work with \code{\link{massbic}}, meaning they contain
      a list with named element \code{fit}, which is the return value of 
      \code{\link{fitnMm}}. Furthermore, all are assumed to have the same 
      dimensions.
}
  \item{savdir}{
      String specifying directory in which \code{\link{files}} are located.
}
  \item{subt}{
      Number of characters to be subtracted from string in \code{files}.
      The default is intended for format ending in "seed=%0.2d.rds".
}
}
\details{
    Suppose you have a directory \code{dir} containing RDS files.
    To create a list of sorted filenames, use for example the code provided 
    here: 

    #files <- list.files(dir, pattern=".rds")
    #    p1 <- c("pattern1a", "pattern1b") ## search patterns
    #    p2 <- c("pattern2a", "pattern2b")

    #    filelist <- list()
    #    for (i in seq_along(p1)) {
    #        for (j in seq_along(p2)) {
    #            # for lack of AND matching, OR match everything else and invert
    #            r <- grep(paste(p1[-i], p2[-j], sep="|", collapse="|"), 
    #                      files, value=TRUE, invert=TRUE)
    #            filelist[[paste0(p1[i], p2[j])]] <- r
    #        }
    #    }


    This will create a list of character vectors, matching filenames by p1 and 
    p2.
    % TODO: improve clarity here.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    No return value. Intended side effects:
    Produces \code{3*length(files)} PDFs named as follows:
    From each entry in \code{files} takes first string, subtracts \code{subt}
    characters and appends \code{".pdf", "_mcl.pdf"} and \code{"_comp.pdf"}
}
\references{
}
\author{
    Nicolas Trutmann
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{massplot}},
    \code{\link{compplot}},
    \code{\link{fitnMm}}
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
