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
c("dbSNP_GWAS", "dbSNP_Common"), GR.Gene = c("UCSC_knownGene",
"UCSC_knownCanonical"), verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata")
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
function (data, distance.max = 2e+05, decay.kernel = c("rapid",
"slow", "linear", "constant"), decay.exponent = 2, GR.SNP =
c("dbSNP_GWAS",
"dbSNP_Common"), GR.Gene = c("UCSC_knownGene", "UCSC_knownCanonical"),
verbose = T, RData.location = "http://galahad.well.ox.ac.uk/bigdata")
{
decay.kernel <- match.arg(decay.kernel)
data <- gsub("_", ":", data, perl = T)
data <- gsub("imm:", "chr", data, perl = T)
data <- unique(data)
gr_SNP <- xSNPlocations(data = data, GR.SNP = GR.SNP, verbose =
verbose,
RData.location = RData.location)
if (verbose) {
now <- Sys.time()
message(sprintf("Load positional information for Genes (\%s) ...",
as.character(now)), appendLF = T)
}
if (class(GR.Gene) == "GRanges") {
gr_Gene <- GR.Gene
}
else {
gr_Gene <- xRDataLoader(RData.customised = GR.Gene[1],
verbose = verbose, RData.location = RData.location)
if (is.null(gr_Gene)) {
GR.Gene <- "UCSC_knownGene"
if (verbose) {
message(sprintf("Instead, \%s will be used", GR.Gene),
appendLF = T)
}
gr_Gene <- xRDataLoader(RData.customised = GR.Gene,
verbose = verbose, RData.location = RData.location)
}
}
if (verbose) {
now <- Sys.time()
message(sprintf("Define nearby genes (\%s) ...", as.character(now)),
appendLF = T)
}
maxgap <- distance.max
minoverlap <- 1L
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
df_nGenes <- data.frame(Gene = names(x), SNP = names(y),
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
SNP = names(y), Dist = dists, stringsAsFactors = F)
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
}
else {
df_nGenes <- NULL
if (verbose) {
now <- Sys.time()
message(sprintf("No nearby genes are defined"), appendLF = T)
}
}
invisible(df_nGenes)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
