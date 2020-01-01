\name{print.ls_eTerm}
\alias{print.ls_eTerm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.ls_eTerm(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
%% ~~Describe \code{x} here~~
}
\item{\dots}{
%% ~~Describe \code{\dots} here~~
}
}
\details{
%% ~~ If necessary, more details than the description above ~~
}
\value{
%% ~Describe the value returned
%% If it is a LIST, use
%% \item{comp1 }{Description of 'comp1'}
%% \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%% ~~who you are~~
}
\note{
%% ~~further notes~~
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
function (x, ...)
{
cat(sprintf("An object of S3 class '\%s', with \%d components:",
class(x), length(names(x))), "\n", sep = "")
cat(sprintf(" $df: a data frame of \%d rows X \%d columns",
dim(x$df)[1], dim(x$df)[2]), "\n", sep = "")
cat(sprintf(" $mat: a data matrix of \%d rows X \%d columns",
dim(x$mat)[1], dim(x$mat)[2]), "\n", sep = "")
cat(sprintf(" $gp: a ggplot object"), "\n", sep = "")
cat("\n--------------------------------------------------\n")
cat("$df:\n")
print(x$df[1:min(2, nrow(x$df)), 1:13], row.names = FALSE)
cat("......\n")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
