\name{xGR2GeneScores}
\alias{xGR2GeneScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGR2GeneScores(data, significance.threshold = 5e-05, score.cap = 10,
build.conversion = c(NA, "hg38.to.hg19", "hg18.to.hg19"), distance.max
= 50000, decay.kernel = c("slow", "linear", "rapid", "constant"),
decay.exponent = 2, GR.Gene = c("UCSC_knownGene",
"UCSC_knownCanonical"), scoring.scheme = c("max", "sum", "sequential"),
verbose = T, RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
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
\item{distance.max}{
%% ~~Describe \code{distance.max} here~~
}
\item{decay.kernel}{
%% ~~Describe \code{decay.kernel} here~~
}
\item{decay.exponent}{
%% ~~Describe \code{decay.exponent} here~~
}
\item{GR.Gene}{
%% ~~Describe \code{GR.Gene} here~~
}
\item{scoring.scheme}{
%% ~~Describe \code{scoring.scheme} here~~
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
function (data, significance.threshold = 5e-05, score.cap = 10,
build.conversion = c(NA, "hg38.to.hg19", "hg18.to.hg19"),
distance.max = 50000, decay.kernel = c("slow", "linear",
"rapid", "constant"), decay.exponent = 2, GR.Gene = c("UCSC_knownGene",
"UCSC_knownCanonical"), scoring.scheme = c("max", "sum",
"sequential"), verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
build.conversion <- match.arg(build.conversion)
decay.kernel <- match.arg(decay.kernel)
scoring.scheme <- match.arg(scoring.scheme)
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
message(sprintf("'xGR2nGenes' is being called to define nearby genes
(\%s):",
as.character(now)), appendLF = T)
message(sprintf("#######################################################",
appendLF = T))
}
df_nGenes <- xGR2nGenes(data = df_GR$GR, format = "chr:start-end",
build.conversion = build.conversion, distance.max = distance.max,
decay.kernel = decay.kernel, decay.exponent = decay.exponent,
GR.Gene = GR.Gene, scoring = F, verbose = verbose, RData.location =
RData.location,
guid = guid)
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("'xGR2nGenes' has been finished (\%s)!",
as.character(now)), appendLF = T)
message(sprintf("#######################################################\n",
appendLF = T))
}
if (1) {
ind <- match(df_nGenes$GR, df_GR$GR)
score <- df_nGenes$Weight * df_GR$Score[ind]
Gene2GR <- data.frame(Gene = df_nGenes$Gene, GR = df_nGenes$GR,
Score = score, stringsAsFactors = FALSE)
Gene2GR <- Gene2GR[order(Gene2GR$Gene, -Gene2GR$Score,
decreasing = FALSE), ]
ls_gene <- split(x = Gene2GR$Score, f = Gene2GR$Gene)
if (scoring.scheme == "max") {
seeds.genes <- sapply(ls_gene, max)
}
else if (scoring.scheme == "sum") {
seeds.genes <- sapply(ls_gene, sum)
}
else if (scoring.scheme == "sequential") {
seeds.genes <- sapply(ls_gene, function(x) {
base::sum(x/base::rank(-x, ties.method = "min"))
})
}
}
else {
allGenes <- sort(df_nGenes$Gene)
allGR <- sort(df_GR$GR)
G2S_n <- xSparseMatrix(df_nGenes[, c("Gene", "GR", "Weight")],
rows = allGenes, columns = allGR, verbose = verbose)
G2S <- G2S_n
ind <- match(colnames(G2S), df_GR$GR)
df_GR <- df_GR[ind, ]
GR_score <- df_GR$Score
names(GR_score) <- colnames(G2S)
mat_GR_score <- matrix(rep(GR_score, each = nrow(G2S)),
nrow = nrow(G2S))
G2S_score <- G2S * mat_GR_score
if (scoring.scheme == "max") {
seeds.genes <- apply(G2S_score, 1, function(x) {
base::max(x)
})
}
else if (scoring.scheme == "sum") {
seeds.genes <- apply(G2S_score, 1, function(x) {
base::sum(x)
})
}
else if (scoring.scheme == "sequential") {
seeds.genes <- apply(G2S_score, 1, function(x) {
base::sum(base::sort(x, decreasing = T)/(1:length(x)))
})
}
}
x <- seeds.genes
pval <- 10^(-1 * x)
df_Gene <- data.frame(Gene = names(seeds.genes), Score = seeds.genes,
Pval = pval, row.names = NULL, stringsAsFactors = F)
df_Gene <- df_Gene[order(df_Gene$Score, decreasing = TRUE),
]
if (verbose) {
now <- Sys.time()
message(sprintf("In summary, \%d Genes are defined as seeds and scored
using '\%s' scoring scheme (\%s)",
length(seeds.genes), scoring.scheme, as.character(now)),
appendLF = T)
}
df_GR <- df_GR[order(df_GR$Score, df_GR$GR, decreasing = TRUE),
]
mSeed <- list(GR = df_GR, Gene = df_Gene, Call = match.call())
class(mSeed) <- "mSeed"
invisible(mSeed)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
