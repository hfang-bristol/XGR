\name{ls_eTerm}
\alias{ls_eTerm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ls_eTerm(df, mat, gp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{df}{
%% ~~Describe \code{df} here~~
}
\item{mat}{
%% ~~Describe \code{mat} here~~
}
\item{gp}{
%% ~~Describe \code{gp} here~~
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
function (df, mat, gp)
{
if (class(df) != "data.frame" | class(mat) != "matrix" |
all(class(gp) \%in\% c("ggplot", "gg"))) {
stop("The S3 class 'ls_eTerm' object failed to pass integrity
checks!\n")
}
value <- list(df = df, mat = mat, gp = gp)
class(value) <- "ls_eTerm"
return(value)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
