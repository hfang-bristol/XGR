\name{eTerm}
\alias{eTerm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eTerm(term_info, annotation, g, data, background, overlap, fc, zscore,
pvalue, adjp, cross)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{term_info}{
%% ~~Describe \code{term_info} here~~
}
\item{annotation}{
%% ~~Describe \code{annotation} here~~
}
\item{g}{
%% ~~Describe \code{g} here~~
}
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{background}{
%% ~~Describe \code{background} here~~
}
\item{overlap}{
%% ~~Describe \code{overlap} here~~
}
\item{fc}{
%% ~~Describe \code{fc} here~~
}
\item{zscore}{
%% ~~Describe \code{zscore} here~~
}
\item{pvalue}{
%% ~~Describe \code{pvalue} here~~
}
\item{adjp}{
%% ~~Describe \code{adjp} here~~
}
\item{cross}{
%% ~~Describe \code{cross} here~~
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
function (term_info, annotation, g, data, background, overlap,
fc, zscore, pvalue, adjp, cross)
{
if (class(term_info) != "data.frame" | class(g) != "igraph") {
stop("The S3 class 'eTerm' object failed to pass integrity checks!\n")
}
value <- list(term_info = term_info, annotation = annotation,
g = g, data = data, background = background, overlap = overlap,
fc = fc, zscore = zscore, pvalue = pvalue, adjp = adjp,
cross = cross)
class(value) <- "eTerm"
return(value)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
