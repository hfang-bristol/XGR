\name{xSM2DF}
\alias{xSM2DF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSM2DF(data, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
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
function (data, verbose = TRUE)
{
if (class(data) == "dgCMatrix" | is.data.frame(data)) {
}
if (is.null(rownames(data))) {
names_row <- 1:nrow(data)
}
else {
names_row <- rownames(data)
}
if (is.null(colnames(data))) {
names_col <- 1:ncol(data)
}
else {
names_col <- colnames(data)
}
if (is.data.frame(data) | class(data) == "matrix" | class(data) ==
"dgeMatrix") {
data <- as.matrix(data)
xy <- which(data != 0, arr.ind = TRUE)
ind <- which(data != 0, arr.ind = FALSE)
if (nrow(xy) > 0) {
res_df <- data.frame(rownames = names_row[xy[, 1]],
colnames = names_col[xy[, 2]], values = data[ind],
stringsAsFactors = FALSE)
res_df <- res_df[order(res_df$rownames, res_df$values,
decreasing = FALSE), ]
if (verbose) {
message(sprintf("A matrix of \%d X \%d has been converted into a data
frame of \%d X 3.",
dim(data)[1], dim(data)[2], nrow(res_df)),
appendLF = T)
}
}
else {
res_df <- NULL
}
}
else if (class(data) == "dgCMatrix" | class(data) == "dsCMatrix") {
ijx <- summary(data)
if (nrow(ijx) > 0) {
res_df <- data.frame(rownames = names_row[ijx[, 1]],
colnames = names_col[ijx[, 2]], values = ijx[,
3], stringsAsFactors = FALSE)
res_df <- res_df[order(res_df$rownames, res_df$values,
decreasing = FALSE), ]
if (verbose) {
message(sprintf("A matrix of \%d X \%d has been converted into a data
frame of \%d X 3.",
dim(data)[1], dim(data)[2], nrow(res_df)),
appendLF = T)
}
}
else {
res_df <- NULL
}
}
invisible(res_df)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
