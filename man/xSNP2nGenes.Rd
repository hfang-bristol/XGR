\name{xSNP2nGenes}
\alias{xSNP2nGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSNP2nGenes(data, distance.max = 2e+05, decay.kernel = c("rapid",
"slow", "linear", "constant"), decay.exponent = 2, GR.SNP =
c("dbSNP_GWAS", "dbSNP_Common", "dbSNP_Single"), GR.Gene =
c("UCSC_knownGene", "UCSC_knownCanonical"), include.TAD = c("none",
"GM12878", "IMR90", "MSC", "TRO", "H1", "MES", "NPC"), verbose = T,
RData.location = "http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
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
function (data, distance.max = 2e+05, decay.kernel = c("rapid",
"slow", "linear", "constant"), decay.exponent = 2, GR.SNP =
c("dbSNP_GWAS",
"dbSNP_Common", "dbSNP_Single"), GR.Gene = c("UCSC_knownGene",
"UCSC_knownCanonical"), include.TAD = c("none", "GM12878",
"IMR90", "MSC", "TRO", "H1", "MES", "NPC"), verbose = T,
RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
decay.kernel <- match.arg(decay.kernel)
include.TAD <- match.arg(include.TAD)
data <- gsub("_", ":", data, perl = T)
data <- gsub("imm:", "chr", data, perl = T)
data <- unique(data)
gr_SNP <- xSNPlocations(data = data, GR.SNP = GR.SNP, verbose =
verbose,
RData.location = RData.location, guid = guid)
if (verbose) {
now <- Sys.time()
message(sprintf("Load positional information for Genes (\%s) ...",
as.character(now)), appendLF = T)
}
if (class(GR.Gene) == "GRanges") {
gr_Gene <- GR.Gene
}
else {
gr_Gene <- xRDataLoader(GR.Gene[1], verbose = verbose,
RData.location = RData.location, guid = guid)
if (is.null(gr_Gene)) {
GR.Gene <- "UCSC_knownGene"
if (verbose) {
message(sprintf("Instead, \%s will be used", GR.Gene),
appendLF = T)
}
gr_Gene <- xRDataLoader(GR.Gene, verbose = verbose,
RData.location = RData.location, guid = guid)
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
query <- gr_SNP
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
df_y <- GenomicRanges::as.data.frame(y, row.names = NULL)
df_x <- GenomicRanges::as.data.frame(x, row.names = NULL)
if (0) {
df_interval <- data.frame(seqnames = df_y$seqnames,
start = df_y$start, end = df_y$end, stringsAsFactors = FALSE)
ind <- df_y$start < df_x$start
df_interval[ind, ] <- data.frame(seqnames = df_y$seqnames[ind],
start = df_y$start[ind], end = df_x$start[ind],
stringsAsFactors = FALSE)
ind <- df_y$start > df_x$end
df_interval[ind, ] <- data.frame(seqnames = df_y$seqnames[ind],
start = df_x$end[ind], end = df_y$start[ind],
stringsAsFactors = FALSE)
}
else {
df_interval <- data.frame(seqnames = df_y$seqnames,
start = df_y$start, end = df_y$end, stringsAsFactors = FALSE)
ind <- df_x$strand == "+"
df_interval[ind, ] <- data.frame(seqnames = df_y$seqnames[ind],
start = df_y$start[ind], end = df_x$start[ind],
stringsAsFactors = FALSE)
ind <- df_x$strand == "-"
df_interval[ind, ] <- data.frame(seqnames = df_y$seqnames[ind],
start = df_y$start[ind], end = df_x$end[ind],
stringsAsFactors = FALSE)
ind <- df_x$strand == "*"
df_interval[ind, ] <- data.frame(seqnames = df_y$seqnames[ind],
start = df_y$start[ind], end = (df_x$start[ind] +
df_x$end[ind])/2, stringsAsFactors = FALSE)
ind <- df_interval$start > df_interval$end
df_interval[ind, ] <- data.frame(seqnames = df_interval$seqnames[ind],
start = df_interval$end[ind], end = df_interval$start[ind],
stringsAsFactors = FALSE)
}
df_interval[dists == 0, "start"] <- df_y$start[dists ==
0]
df_interval[dists == 0, "end"] <- df_y$end[dists ==
0]
vec_interval <- paste0(df_interval$seqnames, ":",
as.character(df_interval$start), "-", as.character(df_interval$end))
df_nGenes <- data.frame(Gene = names(x), SNP = names(y),
Dist = dists, Gap = vec_interval, stringsAsFactors = F)
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
df_nGenes <- df_nGenes[, c("Gene", "SNP", "Dist", "Weight",
"Gap")]
df_nGenes <- df_nGenes[order(df_nGenes$Gene, df_nGenes$Dist,
decreasing = FALSE), ]
}
else {
df_nGenes <- NULL
if (verbose) {
now <- Sys.time()
message(sprintf("No nearby genes are defined"), appendLF = T)
}
}
default.include.TAD <- c("GM12878", "IMR90", "MSC", "TRO",
"H1", "MES", "MES")
ind <- match(default.include.TAD, include.TAD)
include.TAD <- default.include.TAD[!is.na(ind)]
if (length(include.TAD) > 0) {
if (verbose) {
now <- Sys.time()
message(sprintf("Inclusion of TAD boundary regions is based on '\%s'",
include.TAD), appendLF = T)
}
df_nGenes$TAD <- rep("Excluded", nrow(df_nGenes))
TAD <- xRDataLoader(paste0("TAD.", include.TAD), RData.location =
RData.location,
guid = guid, verbose = verbose)
iGR <- xGR(data = df_nGenes$Gap, format = "chr:start-end",
RData.location = RData.location, guid = guid)
q2r <-
as.matrix(as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query
= iGR,
subject = TAD, type = "within", select = "all", ignore.strand = T))))
q2r <- q2r[!duplicated(q2r[, 1]), ]
df_nGenes$TAD[q2r[, 1]] <- GenomicRanges::mcols(TAD)[q2r[,
2], ]
df_nGenes$TAD[df_nGenes$Dist == 0] <- df_nGenes$Gap[df_nGenes$Dist ==
0]
if (verbose) {
now <- Sys.time()
message(sprintf("\t\%d out of \%d SNP-nGene pairs are within the same
TAD boundary regions",
sum(df_nGenes$TAD != "Excluded"), length(iGR)),
appendLF = T)
message(sprintf("\t\%d out of \%d genes are defined as nearby genes
after considering TAD boundary regions",
length(unique(df_nGenes[df_nGenes$TAD != "Excluded",
"Gene"])), length(unique(df_nGenes$Gene))),
appendLF = T)
}
}
if (!is.null(df_nGenes)) {
ind <- xSymbol2GeneID(df_nGenes$Gene, details = FALSE,
verbose = verbose, RData.location = RData.location,
guid = guid)
df_nGenes <- df_nGenes[!is.na(ind), ]
if (nrow(df_nGenes) == 0) {
df_nGenes <- NULL
}
}
invisible(df_nGenes)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
