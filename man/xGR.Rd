\name{xGR}
\alias{xGR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGR(data, format = c("chr:start-end", "data.frame", "bed", "GRanges"),
build.conversion = c(NA, "hg38.to.hg19", "hg18.to.hg19"), verbose = T,
RData.location = "http://galahad.well.ox.ac.uk/bigdata")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{format}{
%% ~~Describe \code{format} here~~
}
\item{build.conversion}{
%% ~~Describe \code{build.conversion} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
}
\item{RData.location}{
%% ~~Describe \code{RData.location} here~~
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
function (data, format = c("chr:start-end", "data.frame", "bed",
"GRanges"), build.conversion = c(NA, "hg38.to.hg19", "hg18.to.hg19"),
verbose = T, RData.location = "http://galahad.well.ox.ac.uk/bigdata")
{
format <- match.arg(format)
build.conversion <- match.arg(build.conversion)
if (is.matrix(data) | is.data.frame(data) | class(data) ==
"GRanges") {
data <- data
}
else if (!is.null(data) & any(!is.na(data))) {
if (length(data) == 1) {
if (file.exists(data)) {
data <- utils::read.delim(file = data, header = F,
row.names = NULL, stringsAsFactors = F)
data <- unique(data[, 1])
}
else {
data <- data
}
}
else {
data <- data
}
}
else {
stop("The file 'data' must be provided!\n")
}
if (format == "data.frame") {
if (ncol(data) >= 3) {
data <- data
}
else if (ncol(data) == 2) {
data <- cbind(data, data[, 2])
}
else {
stop("Your input 'data.file' is not as expected!\n")
}
ind <- suppressWarnings(which(!is.na(as.numeric(data[,
2])) & !is.na(as.numeric(data[, 3]))))
dGR <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(data[ind,
1]), ranges = IRanges::IRanges(start = as.numeric(data[ind,
2]), end = as.numeric(data[ind, 3])), strand = S4Vectors::Rle(rep("*",
length(ind))))
names(dGR) <- paste(data[, 1], ":", data[, 2], "-", data[,
3], sep = "")
}
else if (format == "chr:start-end") {
data <- unique(data)
input <- do.call(rbind, strsplit(data, ":|-"))
if (ncol(input) >= 3) {
data <- matrix(input[, 1:3], nrow = nrow(input))
}
else if (ncol(input) == 2) {
data <- matrix(input[, c(1, 2, 2)], nrow = nrow(input))
}
else {
stop("Your input 'data' does not meet the format 'chr:start-end'!\n")
}
ind <- suppressWarnings(which(!is.na(as.numeric(data[,
2])) & !is.na(as.numeric(data[, 3]))))
dGR <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(data[ind,
1]), ranges = IRanges::IRanges(start = as.numeric(data[ind,
2]), end = as.numeric(data[ind, 3])), strand = S4Vectors::Rle(rep("*",
length(ind))))
names(dGR) <- paste(data[, 1], ":", data[, 2], "-", data[,
3], sep = "")
}
else if (format == "bed") {
ind <- suppressWarnings(which(!is.na(as.numeric(data[,
2])) & !is.na(as.numeric(data[, 3]))))
dGR <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(data[ind,
1]), ranges = IRanges::IRanges(start = as.numeric(data[ind,
2]) + 1, end = as.numeric(data[ind, 3])), strand =
S4Vectors::Rle(rep("*",
length(ind))))
names(dGR) <- paste(data[, 1], ":", data[, 2], "-", data[,
3], sep = "")
}
else if (format == "GRanges") {
dGR <- data
if (is.null(names(dGR))) {
df <- as.data.frame(dGR, row.names = NULL)
names(dGR) <- paste(df$seqnames, ":", df$start, "-",
df$end, sep = "")
}
}
if (!is.na(build.conversion)) {
if (verbose) {
message(sprintf("\tdata genomic regions: lifted over via genome build
conversion `\%s`",
build.conversion), appendLF = T)
}
dGR <- xLiftOver(data.file = dGR, format.file = "GRanges",
build.conversion = build.conversion, merged = F,
verbose = verbose, RData.location = RData.location)
}
invisible(dGR)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
