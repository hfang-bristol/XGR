\name{xGR2xGeneScores}
\alias{xGR2xGeneScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGR2xGeneScores(data, significance.threshold = NULL, score.cap = NULL,
build.conversion = c(NA, "hg38.to.hg19", "hg18.to.hg19"), crosslink =
c("genehancer", "PCHiC_combined", "GTEx_V6p_combined", "nearby"),
crosslink.customised = NULL, cdf.function = c("original", "empirical"),
scoring.scheme = c("max", "sum", "sequential"), nearby.distance.max =
50000, nearby.decay.kernel = c("rapid", "slow", "linear", "constant"),
nearby.decay.exponent = 2, verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
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
\item{build.conversion}{
%% ~~Describe \code{build.conversion} here~~
}
\item{crosslink}{
%% ~~Describe \code{crosslink} here~~
}
\item{crosslink.customised}{
%% ~~Describe \code{crosslink.customised} here~~
}
\item{cdf.function}{
%% ~~Describe \code{cdf.function} here~~
}
\item{scoring.scheme}{
%% ~~Describe \code{scoring.scheme} here~~
}
\item{nearby.distance.max}{
%% ~~Describe \code{nearby.distance.max} here~~
}
\item{nearby.decay.kernel}{
%% ~~Describe \code{nearby.decay.kernel} here~~
}
\item{nearby.decay.exponent}{
%% ~~Describe \code{nearby.decay.exponent} here~~
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
function (data, significance.threshold = NULL, score.cap = NULL,
build.conversion = c(NA, "hg38.to.hg19", "hg18.to.hg19"),
crosslink = c("genehancer", "PCHiC_combined", "GTEx_V6p_combined",
"nearby"), crosslink.customised = NULL, cdf.function = c("original",
"empirical"), scoring.scheme = c("max", "sum", "sequential"),
nearby.distance.max = 50000, nearby.decay.kernel = c("rapid",
"slow", "linear", "constant"), nearby.decay.exponent = 2,
verbose = T, RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
build.conversion <- match.arg(build.conversion)
cdf.function <- match.arg(cdf.function)
scoring.scheme <- match.arg(scoring.scheme)
nearby.decay.kernel <- match.arg(nearby.decay.kernel)
if (verbose) {
now <- Sys.time()
message(sprintf("\n#######################################################",
appendLF = T))
message(sprintf("'xGRscores' is being called to score GR (\%s):",
as.character(now)), appendLF = T)
message(sprintf("#######################################################",
appendLF = T))
}
df_GR <- xGRscores(data = data, significance.threshold =
significance.threshold,
score.cap = score.cap, verbose = verbose)
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("'xGRscores' has been finished (\%s)!",
as.character(now)), appendLF = T)
message(sprintf("#######################################################\n",
appendLF = T))
}
if (verbose) {
now <- Sys.time()
message(sprintf("\n#######################################################",
appendLF = T))
message(sprintf("'xGR2xGenes' is being called to define crosslinked
genes (\%s):",
as.character(now)), appendLF = T)
message(sprintf("#######################################################",
appendLF = T))
}
df_xGenes <- xGR2xGenes(data = df_GR$GR, format = "chr:start-end",
build.conversion = build.conversion, crosslink = crosslink,
crosslink.customised = crosslink.customised, cdf.function =
cdf.function,
scoring = F, nearby.distance.max = nearby.distance.max,
nearby.decay.kernel = nearby.decay.kernel, nearby.decay.exponent =
nearby.decay.exponent,
verbose = verbose, silent = !verbose, RData.location = RData.location,
guid = guid)
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("'xGR2xGenes' has been finished (\%s)!",
as.character(now)), appendLF = T)
message(sprintf("#######################################################\n",
appendLF = T))
}
if (1) {
ind <- match(df_xGenes$GR, df_GR$GR)
score <- df_xGenes$Score * df_GR$Score[ind]
df_link <- data.frame(GR = df_xGenes$GR, Gene = df_xGenes$Gene,
Score = score, Score_GR = df_GR$Score[ind], Score_link =
df_xGenes$Score,
stringsAsFactors = FALSE)
df_link <- df_link[order(-df_link$Score, decreasing = FALSE),
]
ind <- xGRsort(df_link$GR)
df_link <- df_link[ind, ]
Gene2GR <- data.frame(Gene = df_xGenes$Gene, GR = df_xGenes$GR,
Score = score, stringsAsFactors = FALSE)
Gene2GR <- Gene2GR[order(Gene2GR$Gene, -Gene2GR$Score,
decreasing = FALSE), ]
Gene <- Score <- NULL
if (scoring.scheme == "max") {
summaryFun <- max
}
else if (scoring.scheme == "sum") {
summaryFun <- sum
}
else if (scoring.scheme == "sequential") {
summaryFun <- function(x) {
base::sum(x/base::rank(-x, ties.method = "min"))
}
}
tmp <- as.data.frame(df_xGenes \%>\% dplyr::group_by(Gene) \%>\%
dplyr::summarise(Score = summaryFun(Score)))
seeds.genes <- tmp$Score
names(seeds.genes) <- tmp$Gene
}
rescaleFun <- function(x) {
0.100001 + 0.9 * 0.99999888888 * (x - min(x))/(max(x) -
min(x))
}
x <- rescaleFun(seeds.genes)
pval <- 10^(-x * 10)
df_Gene <- data.frame(Gene = names(seeds.genes), Score = seeds.genes,
Pval = pval, row.names = NULL, stringsAsFactors = F)
df_Gene <- df_Gene[order(df_Gene$Score, decreasing = TRUE),
]
if (verbose) {
now <- Sys.time()
message(sprintf("In summary, \%d Genes are defined as seeds from a list
of \%d (out of \%d) genomic regions and scored using '\%s' scoring
scheme",
length(seeds.genes), length(unique(df_xGenes$GR)),
nrow(df_GR), scoring.scheme, as.character(now)),
appendLF = T)
}
df_GR <- df_GR[order(df_GR$Score, df_GR$GR, decreasing = TRUE),
]
mSeed <- list(GR = df_GR, Gene = df_Gene, Link = df_link)
class(mSeed) <- "mSeed"
invisible(mSeed)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
