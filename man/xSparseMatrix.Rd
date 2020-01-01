\name{xSparseMatrix}
\alias{xSparseMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSparseMatrix(input.file, rows = NULL, columns = NULL, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{input.file}{
%% ~~Describe \code{input.file} here~~
}
\item{rows}{
%% ~~Describe \code{rows} here~~
}
\item{columns}{
%% ~~Describe \code{columns} here~~
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
function (input.file, rows = NULL, columns = NULL, verbose = T)
{
if (is.matrix(input.file) | is.data.frame(input.file)) {
if (verbose) {
message(sprintf("Load the input data matrix of \%d X \%d ...",
dim(input.file)[1], dim(input.file)[2]), appendLF = T)
}
if (ncol(input.file) == 2) {
input.file <- cbind(input.file, rep(1, nrow(input.file)))
}
if (is.data.frame(input.file)) {
input <- cbind(as.character(input.file[, 1]), as.character(input.file[,
2]), as.character(input.file[, 3]))
}
else {
input <- input.file
}
}
else if (is.character(input.file) & sum(input.file != "")) {
if (verbose) {
message(sprintf("Read the input file '\%s' ...", input.file),
appendLF = T)
}
input <- utils::read.delim(input.file, header = T, sep = "\t",
colClasses = "character")
if (ncol(input) == 2) {
input <- cbind(input, rep(1, nrow(input)))
}
}
else {
return(NULL)
}
x <- input
if (!is.null(x)) {
if (is.null(rows)) {
x_row <- sort(unique(x[, 1]))
}
else {
x_row <- rows
}
if (is.null(columns)) {
x_col <- sort(unique(x[, 2]))
}
else {
x_col <- columns
}
ind_row <- match(x[, 1], x_row)
ind_col <- match(x[, 2], x_col)
x.sparse <- Matrix::sparseMatrix(i = ind_row, j = ind_col,
x = as.numeric(x[, 3]), dims = c(length(x_row), length(x_col)))
rownames(x.sparse) <- x_row
colnames(x.sparse) <- x_col
if (verbose) {
message(sprintf("There are \%d entries, converted into a sparse matrix
of \%d X \%d.",
nrow(x), dim(x.sparse)[1], dim(x.sparse)[2]),
appendLF = T)
}
}
else {
return(NULL)
}
invisible(x.sparse)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
