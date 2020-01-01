\name{xAuxRdWrap}
\alias{xAuxRdWrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xAuxRdWrap(path = "./XGR/man", remove.dontrun = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{path}{
%% ~~Describe \code{path} here~~
}
\item{remove.dontrun}{
%% ~~Describe \code{remove.dontrun} here~~
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
function (path = "./XGR/man", remove.dontrun = FALSE)
{
for (nm in list.files(path, pattern = "\\.Rd$")) {
input_file <- file.path(path, nm)
message(sprintf(input_file), appendLF = T)
x <- readLines(input_file)
y <- vector()
flag_examples <- F
flag_examples_dontrun <- F
num_examples_brace <- 0
for (i in 1:length(x)) {
if (length(grep("^#", x[i])) == 0 & length(grep("^\\\\title",
x[i])) == 0) {
y <- c(y, strwrap(x[i]))
}
else {
y <- c(y, x[i])
}
if (remove.dontrun == T) {
if (x[i] == "\\examples{") {
flag_examples <- T
}
if (flag_examples == T) {
if (x[i] == "\\dontrun{") {
y <- y[-1 * length(y)]
flag_examples_dontrun <- T
}
if (flag_examples_dontrun == T) {
if (x[i] == "}") {
num_examples_brace <- num_examples_brace +
1
}
if (num_examples_brace == 1) {
y <- y[-1 * length(y)]
flag_examples <- F
}
}
}
}
}
fileConn <- file(input_file)
writeLines(y, fileConn)
close(fileConn)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
