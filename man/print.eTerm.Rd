\name{print.eTerm}
\alias{print.eTerm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.eTerm(x, ...)
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
cat(sprintf("An object of S3 class '\%s', with \%d components
including:",
class(x), length(names(x))), "\n", sep = "")
cat(sprintf(" $term_info: a data frame of \%d rows X \%d columns",
dim(x$term_info)[1], dim(x$term_info)[2]), "\n", sep = "")
cat(sprintf(" $data: a vector (\%d in total)", length(x$data)),
"\n", sep = "")
cat(sprintf(" $background: a vector (\%d in total)",
length(x$background)),
"\n", sep = "")
cat(sprintf(" $adjp: a vector (\%d in total)", length(x$adjp)),
"\n", sep = "")
cat(sprintf(" $cross: a matrix of \%d X \%d", dim(x$cross)[1],
dim(x$cross)[2]), "\n", sep = "")
cat(sprintf(" $g: an 'igraph' object"), "\n", sep = "")
cat(sprintf(" $g$ontology: '\%s'", x$g$ontology), "\n", sep = "")
cat("\n--------------------------------------------------\n")
cat("xEnrichViewer(eTerm):\n")
print(xEnrichViewer(x), row.names = TRUE)
cat("......\n")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
