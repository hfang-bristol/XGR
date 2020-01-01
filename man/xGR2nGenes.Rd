\name{xGR2nGenes}
\alias{xGR2nGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGR2nGenes(data, format = c("chr:start-end", "data.frame", "bed",
"GRanges"), build.conversion = c(NA, "hg38.to.hg19", "hg18.to.hg19"),
distance.max = 50000, decay.kernel = c("rapid", "slow", "linear",
"constant"), decay.exponent = 2, GR.Gene = c("UCSC_knownGene",
"UCSC_knownCanonical"), scoring = F, scoring.scheme = c("max", "sum",
"sequential"), scoring.rescale = F, verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
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
\item{scoring}{
%% ~~Describe \code{scoring} here~~
}
\item{scoring.scheme}{
%% ~~Describe \code{scoring.scheme} here~~
}
\item{scoring.rescale}{
%% ~~Describe \code{scoring.rescale} here~~
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
function (data, format = c("chr:start-end", "data.frame", "bed",
"GRanges"), build.conversion = c(NA, "hg38.to.hg19", "hg18.to.hg19"),
distance.max = 50000, decay.kernel = c("rapid", "slow", "linear",
"constant"), decay.exponent = 2, GR.Gene = c("UCSC_knownGene",
"UCSC_knownCanonical"), scoring = F, scoring.scheme = c("max",
"sum", "sequential"), scoring.rescale = F, verbose = T,
RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
format <- match.arg(format)
build.conversion <- match.arg(build.conversion)
decay.kernel <- match.arg(decay.kernel)
scoring.scheme <- match.arg(scoring.scheme)
dGR <- xGR(data = data, format = format, build.conversion =
build.conversion,
verbose = verbose, RData.location = RData.location, guid = guid)
if (verbose) {
now <- Sys.time()
message(sprintf("Load positional information for Genes (\%s) ...",
as.character(now)), appendLF = T)
}
if (class(GR.Gene) == "GRanges") {
gr_Gene <- xGR(GR.Gene, format = "GRanges", build.conversion =
build.conversion,
verbose = verbose, RData.location = RData.location,
guid = guid)
}
else {
gr_Gene <- xRDataLoader(RData.customised = GR.Gene[1],
verbose = verbose, RData.location = RData.location,
guid = guid)
if (is.null(gr_Gene)) {
GR.Gene <- "UCSC_knownGene"
if (verbose) {
message(sprintf("Instead, \%s will be used", GR.Gene),
appendLF = T)
}
gr_Gene <- xRDataLoader(RData.customised = GR.Gene,
verbose = verbose, RData.location = RData.location,
guid = guid)
}
}
if (verbose) {
now <- Sys.time()
message(sprintf("Define nearby genes (\%s) ...", as.character(now)),
appendLF = T)
}
maxgap <- distance.max - 1
minoverlap <- 0L
subject <- gr_Gene
query <- dGR
q2r <-
as.matrix(as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query
= query,
subject = subject, maxgap = maxgap, minoverlap = minoverlap,
type = "any", select = "all", ignore.strand = T))))
if (length(q2r) > 0) {
if (verbose) {
now <- Sys.time()
message(sprintf("Calculate distance (\%s) ...", as.character(now)),
appendLF = T)
}
if (1) {
x <- subject[q2r[, 2], ]
y <- query[q2r[, 1], ]
dists <- GenomicRanges::distance(x, y, select = "all",
ignore.strand = T)
df_nGenes <- data.frame(Gene = names(x), GR = names(y),
Dist = dists, stringsAsFactors = F)
}
else {
list_gene <- split(x = q2r[, 1], f = q2r[, 2])
ind_gene <- as.numeric(names(list_gene))
res_list <- lapply(1:length(ind_gene), function(i) {
x <- subject[ind_gene[i], ]
y <- query[list_gene[[i]], ]
dists <- GenomicRanges::distance(x, y, select = "all",
ignore.strand = T)
res <- data.frame(Gene = rep(names(x), length(dists)),
GR = names(y), Dist = dists, stringsAsFactors = F)
})
df_nGenes <- do.call(rbind, res_list)
}
if (distance.max == 0) {
x <- df_nGenes$Dist
}
else {
x <- df_nGenes$Dist/distance.max
}
if (decay.kernel == "slow") {
y <- 1 - (x)^decay.exponent
}
else if (decay.kernel == "rapid") {
y <- (1 - x)^decay.exponent
}
else if (decay.kernel == "linear") {
y <- 1 - x
}
else {
y <- 1
}
df_nGenes$Weight <- y
if (verbose) {
now <- Sys.time()
message(sprintf("\%d Genes are defined as nearby genes within \%d(bp)
genomic distance window using '\%s' decay kernel (\%s)",
length(unique(df_nGenes$Gene)), distance.max,
decay.kernel, as.character(now)), appendLF = T)
}
df_nGenes <- df_nGenes[order(df_nGenes$Gene, df_nGenes$Dist,
decreasing = FALSE), ]
if (scoring) {
if (0) {
G2S_score <- xSparseMatrix(df_nGenes[, c("Gene",
"GR", "Weight")], verbose = verbose)
}
ls_gene <- split(x = df_nGenes$Weight, f = df_nGenes$Gene)
if (scoring.scheme == "max") {
seeds.genes <- sapply(ls_gene, max)
if (0) {
ind <- !duplicated(df_nGenes[, 1])
seeds.genes <- df_nGenes[ind, 4]
names(seeds.genes) <- df_nGenes[ind, 1]
system.time({
seeds.genes <- apply(G2S_score, 1, function(x) {
base::max(x)
})
})
}
}
else if (scoring.scheme == "sum") {
seeds.genes <- sapply(ls_gene, sum)
if (0) {
seeds.genes <- base::rowSums(as.matrix(G2S_score))
seeds.genes <- apply(G2S_score, 1, function(x) {
base::sum(x)
})
}
}
else if (scoring.scheme == "sequential") {
seeds.genes <- sapply(ls_gene, function(x) {
base::sum(x/base::rank(-x, ties.method = "min"))
})
if (0) {
seeds.genes <- apply(G2S_score, 1, function(x) {
base::sum(base::sort(x, decreasing = T)/(1:length(x)))
})
}
}
if (verbose) {
now <- Sys.time()
message(sprintf("In summary, \%d Genes are defined as seeds and scored
using '\%s' scoring scheme (\%s)",
length(seeds.genes), scoring.scheme, as.character(now)),
appendLF = T)
}
if (scoring.rescale) {
if (verbose) {
now <- Sys.time()
message(sprintf("Also rescale score into the [0,1] range (\%s)",
as.character(now)), appendLF = T)
}
seeds.genes <- (seeds.genes - min(seeds.genes))/(max(seeds.genes) -
min(seeds.genes))
}
df_Gene <- data.frame(Gene = names(seeds.genes),
Score = seeds.genes, stringsAsFactors = F)
rownames(df_Gene) <- NULL
df_Gene <- df_Gene[order(df_Gene$Score, decreasing = TRUE),
]
invisible(df_Gene)
}
else {
invisible(df_nGenes)
}
}
else {
df_nGenes <- NULL
if (verbose) {
now <- Sys.time()
message(sprintf("No nearby genes are defined"), appendLF = T)
}
invisible(df_nGenes)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
