\name{xReport}
\alias{xReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xReport(obj, rmd = NULL, output_format = NULL, output_file = NULL,
output_dir = NULL, quiet = T, verbose = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{obj}{
%% ~~Describe \code{obj} here~~
}
\item{rmd}{
%% ~~Describe \code{rmd} here~~
}
\item{output_format}{
%% ~~Describe \code{output_format} here~~
}
\item{output_file}{
%% ~~Describe \code{output_file} here~~
}
\item{output_dir}{
%% ~~Describe \code{output_dir} here~~
}
\item{quiet}{
%% ~~Describe \code{quiet} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
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
function (obj, rmd = NULL, output_format = NULL, output_file = NULL,
output_dir = NULL, quiet = T, verbose = T, ...)
{
res <- NULL
if (is.null(output_dir)) {
output_dir <- getwd()
}
if (class(obj) == "eTerm") {
if (is.null(rmd)) {
rmd <- system.file("DynamicReport", "eTerm.Rmd",
package = "XGR")
}
if (file.exists(rmd)) {
rmarkdown::render(rmd, output_format = output_format,
output_file = output_file, output_dir = output_dir,
quiet = quiet, ...)
if (is.null(output_file)) {
output_file <- gsub("\\.rmd$", "\\.html", basename(rmd),
ignore.case = T)
}
res <- sprintf("Congratulations! A html file '\%s' created in the
directory '\%s'",
output_file, output_dir)
if (verbose) {
message(res, appendLF = TRUE)
}
}
}
invisible(res)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
