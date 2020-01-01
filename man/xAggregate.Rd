\name{xAggregate}
\alias{xAggregate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xAggregate(data, bin = F, nbin = 10, scale.log = T, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{bin}{
%% ~~Describe \code{bin} here~~
}
\item{nbin}{
%% ~~Describe \code{nbin} here~~
}
\item{scale.log}{
%% ~~Describe \code{scale.log} here~~
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
function (data, bin = F, nbin = 10, scale.log = T, verbose = T)
{
startT <- Sys.time()
if (verbose) {
message(paste(c("Start at ", as.character(startT)), collapse = ""),
appendLF = TRUE)
message("", appendLF = TRUE)
}
if (is.matrix(data)) {
data <- as.data.frame(data)
}
data_input <- data
if (bin) {
if (verbose) {
message(sprintf("The input data of \%d X \%d is first binned into \%d
ranks per column (\%s) ...",
nrow(data), ncol(data), nbin, as.character(Sys.time())),
appendLF = T)
}
for (j in 1:ncol(data)) {
x <- data[, j]
y <- x[!is.na(x)]
cut_index <- as.numeric(cut(y, breaks = min(y) +
(max(y) - min(y)) * seq(0, 1, length.out = nbin +
1)))
cut_index[is.na(cut_index)] <- 1
x[!is.na(x)] <- cut_index
data[, j] <- x
}
}
else {
if (verbose) {
message(sprintf("The input data of \%d X \%d is directly used for
aggregation (\%s) ...",
nrow(data), ncol(data), as.character(Sys.time())),
appendLF = T)
}
}
num_nonna <- apply(data, 1, function(x) sum(!is.na(x)))
sum_nonna <- apply(data, 1, function(x) sum(x, na.rm = T))
tmp <- max(sum_nonna[!is.infinite(sum_nonna)])
sum_nonna[is.infinite(sum_nonna)] <- tmp
if (scale.log) {
if (verbose) {
message(sprintf("The per-row sum is log-scaled"),
appendLF = T)
}
sum_nonna[sum_nonna == 0] <- NA
sum_nonna <- log(sum_nonna)
}
else {
if (verbose) {
message(sprintf("The per-row sum is WITHOUT log-scaled"),
appendLF = T)
}
}
scale_sum <- (sum_nonna - min(sum_nonna, na.rm = T))/(max(sum_nonna,
na.rm = T) - min(sum_nonna, na.rm = T)) * 0.9999999
data_input$Aggregate <- num_nonna + scale_sum
endT <- Sys.time()
if (verbose) {
message(paste(c("\nFinish at ", as.character(endT)),
collapse = ""), appendLF = TRUE)
}
runTime <- as.numeric(difftime(strptime(endT, "\%Y-\%m-\%d
\%H:\%M:\%S"),
strptime(startT, "\%Y-\%m-\%d \%H:\%M:\%S"), units = "secs"))
message(paste(c("Runtime in total is: ", runTime, " secs\n"),
collapse = ""), appendLF = TRUE)
invisible(data_input)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
