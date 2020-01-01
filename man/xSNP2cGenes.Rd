\name{xSNP2cGenes}
\alias{xSNP2cGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSNP2cGenes(data, entity = c("SNP", "chr:start-end", "data.frame",
"bed", "GRanges"), include.HiC = NA, GR.SNP = c("dbSNP_GWAS",
"dbSNP_Common"), cdf.function = c("empirical", "exponential"), plot =
FALSE, verbose = TRUE, RData.location =
"http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{entity}{
%% ~~Describe \code{entity} here~~
}
\item{include.HiC}{
%% ~~Describe \code{include.HiC} here~~
}
\item{GR.SNP}{
%% ~~Describe \code{GR.SNP} here~~
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
function (data, entity = c("SNP", "chr:start-end", "data.frame",
"bed", "GRanges"), include.HiC = NA, GR.SNP = c("dbSNP_GWAS",
"dbSNP_Common"), cdf.function = c("empirical", "exponential"),
plot = FALSE, verbose = TRUE, RData.location =
"http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
entity <- match.arg(entity)
cdf.function <- match.arg(cdf.function)
df_FTS <- xDefineHIC(data = NULL, include.HiC = include.HiC,
verbose = verbose, RData.location = RData.location, guid = guid)
if (is.null(data)) {
return(NULL)
}
PCHiC <- xDefineHIC(data = data, entity = entity, include.HiC =
include.HiC,
GR.SNP = GR.SNP, verbose = verbose, RData.location = RData.location,
guid = guid)
df_data <- PCHiC$df
if (!is.null(df_FTS)) {
uid <- paste(df_FTS[, 1], df_FTS[, 2], sep = "_")
df <- cbind(uid, df_FTS)
res_list <- split(x = df$score, f = df$uid)
raw_score <- unlist(lapply(res_list, max))
uid_data <- paste(df_data[, 1], df_data[, 2], sep = "_")
if (cdf.function == "exponential") {
lambda <- MASS::fitdistr(raw_score, "exponential")$estimate
wE <- stats::pexp(df_data$score, rate = lambda)
if (nrow(df_data) == 0) {
df_cGenes <- NULL
}
else {
Gene <- sapply(1:nrow(df_data), function(i) {
if (df_data$SNP_end[i] == "bait/from") {
df_data$to_genes[i]
}
else {
df_data$from_genes[i]
}
})
df_cGenes_ori <- data.frame(Gene = Gene, SNP = df_data$SNP,
Sig = df_data$score, Weight = wE, row.names = NULL,
stringsAsFactors = FALSE)
ls_tmp <- strsplit(df_cGenes_ori$Gene, ";")
res_ls <- lapply(1:length(ls_tmp), function(i) {
x <- ls_tmp[[i]]
x <- x[x != "."]
data.frame(Gene = x, df_cGenes_ori[rep(i, length(x)),
-1], row.names = NULL, stringsAsFactors = FALSE)
})
df_cGenes <- do.call(rbind, res_ls)
}
if (plot) {
hist(raw_score, breaks = 1000, freq = FALSE,
col = "grey", xlab = "Score", main = "")
curve(stats::dexp(x = raw_score, rate = lambda),
5:max(raw_score), col = 2, add = TRUE)
}
if (verbose) {
now <- Sys.time()
message(sprintf("HiC weights are CDF of exponential empirical
distributions (parameter lambda=\%f)",
lambda), appendLF = TRUE)
}
}
else if (cdf.function == "empirical") {
my.CDF <- stats::ecdf(raw_score)
wE <- my.CDF(df_data$score)
if (nrow(df) == 0) {
df_cGenes <- NULL
}
else {
Gene <- sapply(1:nrow(df_data), function(i) {
if (df_data$SNP_end[i] == "bait/from") {
df_data$to_genes[i]
}
else {
df_data$from_genes[i]
}
})
df_cGenes_ori <- data.frame(Gene = Gene, SNP = df_data$SNP,
Sig = df_data$score, Weight = wE, row.names = NULL,
stringsAsFactors = FALSE)
ls_tmp <- strsplit(df_cGenes_ori$Gene, ";")
res_ls <- lapply(1:length(ls_tmp), function(i) {
x <- ls_tmp[[i]]
x <- x[x != "."]
data.frame(Gene = x, df_cGenes_ori[rep(i, length(x)),
-1], row.names = NULL, stringsAsFactors = FALSE)
})
df_cGenes <- do.call(rbind, res_ls)
df_cGenes <- df_cGenes[order(df_cGenes$Gene,
-df_cGenes$Sig, df_cGenes$SNP, decreasing = FALSE),
]
}
if (plot) {
plot(my.CDF, xlab = "Score", ylab = "Empirical CDF (HiC weights)",
main = "", xlim = c(0, 50))
}
if (verbose) {
now <- Sys.time()
message(sprintf("HiC weights are CDF of empirical distributions"),
appendLF = TRUE)
}
}
if (verbose) {
now <- Sys.time()
message(sprintf("\%d cGenes are defined involving \%d SNP",
length(unique(df_cGenes$Gene)), length(unique(df_cGenes$SNP))),
appendLF = TRUE)
}
}
else {
df_cGenes <- NULL
if (verbose) {
now <- Sys.time()
message(sprintf("No HiC genes are defined"), appendLF = TRUE)
}
}
if (!is.null(df_cGenes)) {
ind <- xSymbol2GeneID(df_cGenes$Gene, details = FALSE,
verbose = verbose, RData.location = RData.location,
guid = guid)
df_cGenes <- df_cGenes[!is.na(ind), ]
if (nrow(df_cGenes) == 0) {
df_cGenes <- NULL
}
else {
if (1) {
ind <- which(!grepl("^HIST", df_cGenes$Gene))
df_cGenes <- df_cGenes[ind, ]
if (nrow(df_cGenes) == 0) {
df_cGenes <- NULL
}
}
}
}
invisible(df_cGenes)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
