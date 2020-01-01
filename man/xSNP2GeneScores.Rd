\name{xSNP2GeneScores}
\alias{xSNP2GeneScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSNP2GeneScores(data, include.LD = NA, LD.customised = NULL, LD.r2 =
0.8, significance.threshold = 5e-05, score.cap = 10, distance.max =
50000, decay.kernel = c("slow", "linear", "rapid", "constant"),
decay.exponent = 2, GR.SNP = c("dbSNP_GWAS", "dbSNP_Common",
"dbSNP_Single"), GR.Gene = c("UCSC_knownGene", "UCSC_knownCanonical"),
include.TAD = c("none", "GM12878", "IMR90", "MSC", "TRO", "H1", "MES",
"NPC"), scoring.scheme = c("max", "sum", "sequential"), verbose = T,
RData.location = "http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{include.LD}{
%% ~~Describe \code{include.LD} here~~
}
\item{LD.customised}{
%% ~~Describe \code{LD.customised} here~~
}
\item{LD.r2}{
%% ~~Describe \code{LD.r2} here~~
}
\item{significance.threshold}{
%% ~~Describe \code{significance.threshold} here~~
}
\item{score.cap}{
%% ~~Describe \code{score.cap} here~~
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
\item{GR.SNP}{
%% ~~Describe \code{GR.SNP} here~~
}
\item{GR.Gene}{
%% ~~Describe \code{GR.Gene} here~~
}
\item{include.TAD}{
%% ~~Describe \code{include.TAD} here~~
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
function (data, include.LD = NA, LD.customised = NULL, LD.r2 = 0.8,
significance.threshold = 5e-05, score.cap = 10, distance.max = 50000,
decay.kernel = c("slow", "linear", "rapid", "constant"),
decay.exponent = 2, GR.SNP = c("dbSNP_GWAS", "dbSNP_Common",
"dbSNP_Single"), GR.Gene = c("UCSC_knownGene", "UCSC_knownCanonical"),
include.TAD = c("none", "GM12878", "IMR90", "MSC", "TRO",
"H1", "MES", "NPC"), scoring.scheme = c("max", "sum",
"sequential"), verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
decay.kernel <- match.arg(decay.kernel)
include.TAD <- match.arg(include.TAD)
scoring.scheme <- match.arg(scoring.scheme)
if (verbose) {
now <- Sys.time()
message(sprintf("\n#######################################################",
appendLF = T))
message(sprintf("'xSNPscores' is being called to score SNPs (\%s):",
as.character(now)), appendLF = T)
message(sprintf("#######################################################",
appendLF = T))
}
df_SNP <- xSNPscores(data = data, include.LD = include.LD,
LD.customised = LD.customised, LD.r2 = LD.r2, significance.threshold =
significance.threshold,
score.cap = score.cap, verbose = verbose, RData.location =
RData.location,
guid = guid)
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("'xSNPscores' has been finished (\%s)!",
as.character(now)), appendLF = T)
message(sprintf("#######################################################\n",
appendLF = T))
}
if (verbose) {
now <- Sys.time()
message(sprintf("\n#######################################################",
appendLF = T))
message(sprintf("'xSNP2nGenes' is being called to define nearby genes
(\%s):",
as.character(now)), appendLF = T)
message(sprintf("#######################################################",
appendLF = T))
}
df_nGenes <- xSNP2nGenes(data = df_SNP$SNP, distance.max =
distance.max,
decay.kernel = decay.kernel, decay.exponent = decay.exponent,
GR.SNP = GR.SNP, GR.Gene = GR.Gene, include.TAD = include.TAD,
verbose = verbose, RData.location = RData.location, guid = guid)
if (include.TAD != "none") {
TAD <- NULL
df_nGenes <- base::subset(df_nGenes, TAD != "Excluded")
}
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("'xSNP2nGenes' has been finished (\%s)!",
as.character(now)), appendLF = T)
message(sprintf("#######################################################\n",
appendLF = T))
}
allGenes <- sort(df_nGenes$Gene)
allSNPs <- sort(df_SNP$SNP)
G2S_n <- xSparseMatrix(df_nGenes[, c("Gene", "SNP", "Weight")],
rows = allGenes, columns = allSNPs, verbose = verbose)
G2S <- G2S_n
ind <- match(colnames(G2S), df_SNP$SNP)
df_SNP <- df_SNP[ind, ]
SNP_score <- df_SNP$Score
names(SNP_score) <- colnames(G2S)
mat_SNP_score <- matrix(rep(SNP_score, each = nrow(G2S)),
nrow = nrow(G2S))
G2S_score <- G2S * mat_SNP_score
Gene2SNP <- xSM2DF(data = G2S_score, verbose = FALSE)
colnames(Gene2SNP) <- c("Gene", "SNP", "Score")
ls_gene <- split(x = Gene2SNP$Score, f = Gene2SNP$Gene)
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
x <- seeds.genes
pval <- 10^(-1 * x)
df_Gene <- data.frame(Gene = names(seeds.genes), Score = seeds.genes,
Pval = pval, row.names = NULL, stringsAsFactors = F)
df_Gene <- df_Gene[order(df_Gene$Score, decreasing = TRUE),
]
if (verbose) {
now <- Sys.time()
message(sprintf("In summary, \%d Genes are defined as seeds and scored
using '\%s' scoring scheme",
length(seeds.genes), scoring.scheme, as.character(now)),
appendLF = T)
}
df_SNP <- df_SNP[order(df_SNP$Flag, df_SNP$Score, -df_SNP$Pval,
df_SNP$SNP, decreasing = TRUE), ]
mSeed <- list(SNP = df_SNP, Gene = df_Gene, Call = match.call())
class(mSeed) <- "mSeed"
invisible(mSeed)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
