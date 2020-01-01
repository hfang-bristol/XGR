\name{xObjSize}
\alias{xObjSize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xObjSize(obj, type = c("auto", "variable", "file"), units = "auto",
verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{obj}{
%% ~~Describe \code{obj} here~~
}
\item{type}{
%% ~~Describe \code{type} here~~
}
\item{units}{
%% ~~Describe \code{units} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
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
function (obj, type = c("auto", "variable", "file"), units = "auto",
verbose = TRUE)
{
type <- match.arg(type)
res <- NULL
if (type == "auto") {
if (base::file.exists(obj)) {
type <- "file"
}
else {
type <- "variable"
}
}
if (type == "variable") {
if (class(suppressWarnings(try(x <- base::get(obj), T))) \%in\%
"try-error") {
warnings(sprintf("The R variable '\%s' NOT found!",
obj))
return(NULL)
}
else {
res <- base::format(utils::object.size(x), units)
if (verbose) {
message(sprintf("The R variable '\%s' in size: \%s",
obj, res), appendLF = T)
}
}
}
if (type == "file") {
if (base::file.exists(obj)) {
if (verbose) {
message(sprintf("The file '\%s' in size: \%d",
obj, file.info(obj)$size), appendLF = T)
}
}
else {
warnings(sprintf("The file '\%s' NOT found!", obj))
return(NULL)
}
}
invisible(res)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
