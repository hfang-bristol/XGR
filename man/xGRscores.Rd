\name{xGRscores}
\alias{xGRscores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGRscores(data, significance.threshold = 0.05, score.cap = 10, verbose
= T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{significance.threshold}{
%% ~~Describe \code{significance.threshold} here~~
}
\item{score.cap}{
%% ~~Describe \code{score.cap} here~~
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
function (data, significance.threshold = 0.05, score.cap = 10,
verbose = T)
{
if (is.null(data)) {
stop("The input data must be not NULL.\n")
}
else {
if (class(data) == "DataFrame") {
data <- S4Vectors::as.matrix(data)
}
if (is.vector(data)) {
if (length(data) > 1) {
if (is.null(names(data))) {
stop("The input data must have names with attached dbSNP ID.\n")
}
}
else {
data <- utils::read.delim(file = data, header = F,
row.names = NULL, stringsAsFactors = F)
}
}
if (is.vector(data)) {
pval <- data
}
else if (is.matrix(data) | is.data.frame(data)) {
data <- as.matrix(data)
data_list <- split(x = data[, 2], f = as.character(data[,
1]))
res_list <- lapply(data_list, function(x) {
x <- as.numeric(x)
x <- x[!is.na(x)]
if (length(x) > 0) {
min(x)
}
else {
NULL
}
})
pval <- unlist(res_list)
}
tmp <- min(pval[pval != 0])
pval[pval < tmp] <- tmp
}
if (verbose) {
now <- Sys.time()
message(sprintf("A total of \%d GR are input", length(pval)),
appendLF = T)
}
if (is.null(significance.threshold)) {
scores <- log10((1 - pval)/pval)
}
else {
scores <- log10((1 - pval)/pval) - log10((1 -
significance.threshold)/significance.threshold)
}
tmp_max <- max(scores[!is.infinite(scores)])
tmp_min <- min(scores[!is.infinite(scores)])
scores[scores > tmp_max] <- tmp_max
scores[scores < tmp_min] <- tmp_min
ind_remained <- which(scores > 0)
seeds.snps <- scores[ind_remained]
pval <- pval[ind_remained]
if (verbose) {
now <- Sys.time()
message(sprintf("A total of \%d GR are scored positively",
sum(seeds.snps > 0)), appendLF = T)
}
df_GR <- data.frame(GR = names(pval), Score = seeds.snps,
Pval = pval, row.names = NULL, stringsAsFactors = F)
if (!is.null(score.cap)) {
score.cap <- as.numeric(score.cap)
if (score.cap <= max(df_GR$Score)) {
df_GR$Score[df_GR$Score >= score.cap] <- score.cap
if (verbose) {
now <- Sys.time()
message(sprintf("GR score capped to the maximum score \%d.",
score.cap), appendLF = T)
}
}
}
df_GR <- df_GR[order(df_GR$Score, -df_GR$Pval, df_GR$GR,
decreasing = TRUE), ]
invisible(df_GR)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
