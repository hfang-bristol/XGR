\name{xCheckParallel}
\alias{xCheckParallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xCheckParallel(multicores = NULL, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{multicores}{
%% ~~Describe \code{multicores} here~~
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
function (multicores = NULL, verbose = T)
{
flag_parallel <- F
pkgs <- c("doParallel", "foreach")
if (all(pkgs \%in\% rownames(utils::installed.packages()))) {
tmp <- sapply(pkgs, function(pkg) {
requireNamespace(pkg, quietly = T)
})
if (all(tmp)) {
doParallel::registerDoParallel()
cores <- foreach::getDoParWorkers()
if (is.null(multicores)) {
multicores <- max(1, ceiling(cores))
}
else if (is.na(multicores)) {
multicores <- max(1, ceiling(cores))
}
else if (multicores < 1 | multicores > 2 * cores) {
multicores <- max(1, ceiling(cores))
}
else {
multicores <- as.integer(multicores)
}
doParallel::registerDoParallel(cores = multicores)
if (verbose) {
message(sprintf("\tdo parallel computation using \%d cores ...",
multicores, as.character(Sys.time())), appendLF = T)
}
flag_parallel <- T
}
}
return(flag_parallel)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
