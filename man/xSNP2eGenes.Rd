\name{xSNP2eGenes}
\alias{xSNP2eGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSNP2eGenes(data, include.eQTL = NA, eQTL.customised = NULL,
cdf.function = c("empirical", "exponential"), plot = FALSE, verbose =
TRUE, RData.location = "http://galahad.well.ox.ac.uk/bigdata", guid =
NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{include.eQTL}{
%% ~~Describe \code{include.eQTL} here~~
}
\item{eQTL.customised}{
%% ~~Describe \code{eQTL.customised} here~~
}
\item{cdf.function}{
%% ~~Describe \code{cdf.function} here~~
}
\item{plot}{
%% ~~Describe \code{plot} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
}
\item{RData.location}{
%% ~~Describe \code{RData.location} here~~
}
\item{guid}{
%% ~~Describe \code{guid} here~~
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
function (data, include.eQTL = NA, eQTL.customised = NULL, cdf.function
= c("empirical",
"exponential"), plot = FALSE, verbose = TRUE, RData.location =
"http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
cdf.function <- match.arg(cdf.function)
data <- gsub("_", ":", data, perl = TRUE)
data <- gsub("imm:", "chr", data, perl = TRUE)
data <- unique(data)
if (verbose) {
now <- Sys.time()
message(sprintf("A total of \%d SNPs are input", length(data)),
appendLF = TRUE)
}
df_SGS <- xDefineEQTL(data = NULL, include.eQTL = include.eQTL,
eQTL.customised = eQTL.customised, verbose = verbose,
RData.location = RData.location, guid = guid)
if (!is.null(df_SGS)) {
uid <- paste(df_SGS[, 1], df_SGS[, 2], sep = "_")
df <- cbind(uid, df_SGS)
res_list <- split(x = df$Sig, f = df$uid)
vec <- unlist(lapply(res_list, min))
raw_score <- -1 * log10(vec)
if (cdf.function == "exponential") {
lambda <- MASS::fitdistr(raw_score, "exponential")$estimate
ind <- match(df_SGS[, 1], data)
df <- data.frame(df_SGS[!is.na(ind), ])
wE <- stats::pexp(-log10(df$Sig), rate = lambda)
if (nrow(df) == 0) {
df_eGenes <- NULL
}
else {
df_eGenes <- data.frame(Gene = df$Gene, SNP = df$SNP,
Sig = df$Sig, Weight = wE, row.names = NULL,
stringsAsFactors = FALSE)
}
if (plot) {
hist(raw_score, breaks = 1000, freq = FALSE,
col = "grey", xlab = "-log10(p-values)", main = "")
curve(stats::dexp(x = raw_score, rate = lambda),
0:max(raw_score), col = 2, add = TRUE)
}
if (verbose) {
now <- Sys.time()
message(sprintf("eQTL weights are CDF of exponential empirical
distributions (parameter lambda=\%f)",
lambda), appendLF = TRUE)
}
}
else if (cdf.function == "empirical") {
my.CDF <- stats::ecdf(raw_score)
ind <- match(df_SGS[, 1], data)
df <- data.frame(df_SGS[!is.na(ind), ])
wE <- my.CDF(-log10(df$Sig))
if (nrow(df) == 0) {
df_eGenes <- NULL
}
else {
df_eGenes <- data.frame(Gene = df$Gene, SNP = df$SNP,
Sig = df$Sig, Weight = wE, row.names = NULL,
stringsAsFactors = FALSE)
df_eGenes <- df_eGenes[order(df_eGenes$Gene,
df_eGenes$Sig, df_eGenes$SNP, decreasing = FALSE),
]
}
if (plot) {
plot(my.CDF, xlab = "-log10(p-values)", ylab = "Empirical CDF (eQTL
weights)",
main = "")
}
if (verbose) {
now <- Sys.time()
message(sprintf("eQTL weights are CDF of empirical distributions"),
appendLF = TRUE)
}
}
if (verbose) {
now <- Sys.time()
message(sprintf("\%d eGenes are defined involving \%d eQTL",
length(unique(df_eGenes$Gene)), length(unique(df_eGenes$SNP))),
appendLF = TRUE)
}
}
else {
df_eGenes <- NULL
if (verbose) {
now <- Sys.time()
message(sprintf("No eQTL genes are defined"), appendLF = TRUE)
}
}
if (!is.null(df_eGenes)) {
ind <- xSymbol2GeneID(df_eGenes$Gene, details = FALSE,
verbose = verbose, RData.location = RData.location,
guid = guid)
df_eGenes <- df_eGenes[!is.na(ind), ]
if (nrow(df_eGenes) == 0) {
df_eGenes <- NULL
}
}
invisible(df_eGenes)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
