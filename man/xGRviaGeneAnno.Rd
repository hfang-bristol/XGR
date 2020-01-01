\name{xGRviaGeneAnno}
\alias{xGRviaGeneAnno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGRviaGeneAnno(data.file, background.file = NULL, format.file =
c("data.frame", "bed", "chr:start-end", "GRanges"), build.conversion =
c(NA, "hg38.to.hg19", "hg18.to.hg19"), gap.max = 0, GR.Gene =
c("UCSC_knownGene", "UCSC_knownCanonical"), ontology = NA, size.range =
c(10, 2000), min.overlap = 5, which.distance = NULL, test = c("fisher",
"hypergeo", "binomial"), background.annotatable.only = NULL, p.tail =
c("one-tail", "two-tails"), p.adjust.method = c("BH", "BY",
"bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm =
c("none", "pc", "elim", "lea"), elim.pvalue = 0.01, lea.depth = 2,
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
true.path.rule = F, verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data.file}{
%% ~~Describe \code{data.file} here~~
}
\item{background.file}{
%% ~~Describe \code{background.file} here~~
}
\item{format.file}{
%% ~~Describe \code{format.file} here~~
}
\item{build.conversion}{
%% ~~Describe \code{build.conversion} here~~
}
\item{gap.max}{
%% ~~Describe \code{gap.max} here~~
}
\item{GR.Gene}{
%% ~~Describe \code{GR.Gene} here~~
}
\item{ontology}{
%% ~~Describe \code{ontology} here~~
}
\item{size.range}{
%% ~~Describe \code{size.range} here~~
}
\item{min.overlap}{
%% ~~Describe \code{min.overlap} here~~
}
\item{which.distance}{
%% ~~Describe \code{which.distance} here~~
}
\item{test}{
%% ~~Describe \code{test} here~~
}
\item{background.annotatable.only}{
%% ~~Describe \code{background.annotatable.only} here~~
}
\item{p.tail}{
%% ~~Describe \code{p.tail} here~~
}
\item{p.adjust.method}{
%% ~~Describe \code{p.adjust.method} here~~
}
\item{ontology.algorithm}{
%% ~~Describe \code{ontology.algorithm} here~~
}
\item{elim.pvalue}{
%% ~~Describe \code{elim.pvalue} here~~
}
\item{lea.depth}{
%% ~~Describe \code{lea.depth} here~~
}
\item{path.mode}{
%% ~~Describe \code{path.mode} here~~
}
\item{true.path.rule}{
%% ~~Describe \code{true.path.rule} here~~
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
function (data.file, background.file = NULL, format.file =
c("data.frame",
"bed", "chr:start-end", "GRanges"), build.conversion = c(NA,
"hg38.to.hg19", "hg18.to.hg19"), gap.max = 0, GR.Gene =
c("UCSC_knownGene",
"UCSC_knownCanonical"), ontology = NA, size.range = c(10,
2000), min.overlap = 5, which.distance = NULL, test = c("fisher",
"hypergeo", "binomial"), background.annotatable.only = NULL,
p.tail = c("one-tail", "two-tails"), p.adjust.method = c("BH",
"BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm =
c("none",
"pc", "elim", "lea"), elim.pvalue = 0.01, lea.depth = 2,
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
true.path.rule = F, verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
startT <- Sys.time()
message(paste(c("Start at ", as.character(startT)), collapse = ""),
appendLF = T)
message("", appendLF = T)
format.file <- match.arg(format.file)
build.conversion <- match.arg(build.conversion)
test <- match.arg(test)
p.tail <- match.arg(p.tail)
p.adjust.method <- match.arg(p.adjust.method)
ontology.algorithm <- match.arg(ontology.algorithm)
path.mode <- match.arg(path.mode)
if (verbose) {
now <- Sys.time()
message(sprintf("First, import the files formatted as '\%s' (\%s) ...",
format.file, as.character(now)), appendLF = T)
}
if (verbose) {
now <- Sys.time()
message(sprintf("\timport the data file (\%s) ...", as.character(now)),
appendLF = T)
}
if (is.matrix(data.file) | is.data.frame(data.file) | class(data.file)
==
"GRanges") {
data <- data.file
}
else if (!is.null(data.file) & any(!is.na(data.file))) {
if (length(data.file) == 1) {
data <- utils::read.delim(file = data.file, header = F,
row.names = NULL, stringsAsFactors = F)
}
else {
data <- data.file
}
}
else {
stop("The file 'data.file' must be provided!\n")
}
if (verbose) {
now <- Sys.time()
message(sprintf("\timport the background file (\%s) ...",
as.character(now)), appendLF = T)
}
if (is.matrix(background.file) | is.data.frame(background.file) |
class(background.file) == "GRanges") {
background <- background.file
}
else if (!is.null(background.file)) {
if (length(background.file) == 1) {
background <- utils::read.delim(file = background.file,
header = F, row.names = NULL, stringsAsFactors = F)
background <- unique(background[, 1])
}
else {
background <- background.file
}
}
else {
background <- NULL
}
if (verbose) {
now <- Sys.time()
message(sprintf("Second, construct GenomicRanges object (\%s) ...",
as.character(now)), appendLF = T)
}
if (format.file == "data.frame") {
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
data <- data[ind, ]
dGR <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(data[,
1]), ranges = IRanges::IRanges(start = as.numeric(data[,
2]), end = as.numeric(data[, 3])), strand = S4Vectors::Rle(rep("*",
nrow(data))))
if (!is.null(background)) {
if (ncol(background) >= 3) {
background <- background
}
else if (ncol(background) == 2) {
background <- cbind(background, background[,
2])
}
else {
stop("Your input 'background.file' is not as expected!\n")
}
ind <- suppressWarnings(which(!is.na(as.numeric(background[,
2])) & !is.na(as.numeric(background[, 3]))))
background <- background[ind, ]
bGR <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(background[,
1]), ranges = IRanges::IRanges(start = as.numeric(background[,
2]), end = as.numeric(background[, 3])), strand =
S4Vectors::Rle(rep("*",
nrow(background))))
}
else {
bGR <- NULL
}
}
else if (format.file == "chr:start-end") {
if (is.null(dim(data))) {
input <- do.call(rbind, strsplit(data, ":|-"))
}
else {
input <- do.call(rbind, strsplit(data[, 1], ":|-"))
}
if (ncol(input) >= 3) {
data <- input
}
else if (ncol(input) == 2) {
data <- cbind(input, input[, 2])
}
else {
stop("Your input 'data.file' does not meet the format
'chr:start-end'!\n")
}
ind <- suppressWarnings(which(!is.na(as.numeric(data[,
2])) & !is.na(as.numeric(data[, 3]))))
data <- data[ind, ]
dGR <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(data[,
1]), ranges = IRanges::IRanges(start = as.numeric(data[,
2]), end = as.numeric(data[, 3])), strand = S4Vectors::Rle(rep("*",
nrow(data))))
if (!is.null(background)) {
input <- do.call(rbind, strsplit(background[, 1],
":|-"))
if (ncol(input) >= 3) {
background <- input
}
else if (ncol(input) == 2) {
background <- cbind(input, input[, 2])
}
else {
stop("Your input 'background.file' does not meet the format
'chr:start-end'!\n")
}
ind <- suppressWarnings(which(!is.na(as.numeric(background[,
2])) & !is.na(as.numeric(background[, 3]))))
background <- background[ind, ]
bGR <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(background[,
1]), ranges = IRanges::IRanges(start = as.numeric(background[,
2]), end = as.numeric(background[, 3])), strand =
S4Vectors::Rle(rep("*",
nrow(data))))
}
else {
bGR <- NULL
}
}
else if (format.file == "bed") {
ind <- suppressWarnings(which(!is.na(as.numeric(data[,
2])) & !is.na(as.numeric(data[, 3]))))
data <- data[ind, ]
dGR <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(data[,
1]), ranges = IRanges::IRanges(start = as.numeric(data[,
2]) + 1, end = as.numeric(data[, 3])), strand = S4Vectors::Rle(rep("*",
nrow(data))))
if (!is.null(background)) {
ind <- suppressWarnings(which(!is.na(as.numeric(background[,
2])) & !is.na(as.numeric(background[, 3]))))
background <- background[ind, ]
bGR <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(background[,
1]), ranges = IRanges::IRanges(start = as.numeric(background[,
2]) + 1, end = as.numeric(background[, 3])),
strand = S4Vectors::Rle(rep("*", nrow(data))))
}
else {
bGR <- NULL
}
}
else if (format.file == "GRanges") {
dGR <- data
if (!is.null(background)) {
bGR <- background
}
else {
bGR <- NULL
}
}
if (verbose) {
now <- Sys.time()
message(sprintf("Third, define nearby genes of interest and genes as
the background (\%s) ...",
as.character(now)), appendLF = T)
}
if (verbose) {
now <- Sys.time()
message(sprintf("\tload positional information for Genes (\%s) ...",
as.character(now)), appendLF = T)
}
gr_Gene <- xRDataLoader(GR.Gene[1], verbose = verbose, RData.location =
RData.location,
guid = guid)
if (is.null(gr_Gene)) {
GR.Gene <- "UCSC_knownGene"
if (verbose) {
message(sprintf("\tinstead, \%s will be used", GR.Gene),
appendLF = T)
}
gr_Gene <- xRDataLoader(GR.Gene, verbose = verbose, RData.location =
RData.location,
guid = guid)
}
if (!is.na(build.conversion)) {
if (verbose) {
message(sprintf("\tdata genomic regions: lifted over via genome build
conversion `\%s`",
build.conversion), appendLF = T)
}
dGR <- xLiftOver(data.file = dGR, format.file = "GRanges",
build.conversion = build.conversion, merged = F,
verbose = verbose, RData.location = RData.location,
guid = guid)
}
maxgap <- gap.max - 1
minoverlap <- 0L
subject <- gr_Gene
query <- dGR
hits <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query =
query,
subject = subject, maxgap = maxgap, minoverlap = minoverlap,
type = "any", select = "all", ignore.strand = T)))
dGR_genes <- unique(names(gr_Gene[hits[, 2]]))
if (verbose) {
now <- Sys.time()
message(sprintf("\t\%d nearby genes within \%d distance are defined
(\%s) ...",
length(dGR_genes), gap.max, as.character(now)), appendLF = T)
}
if (!is.null(bGR)) {
if (!is.na(build.conversion)) {
if (verbose) {
message(sprintf("\tbackground genomic regions: lifted over via genome
build conversion `\%s`",
build.conversion), appendLF = T)
}
bGR <- xLiftOver(data.file = bGR, format.file = "GRanges",
build.conversion = build.conversion, merged = F,
verbose = verbose, RData.location = RData.location,
guid = guid)
}
maxgap <- gap.max - 1
minoverlap <- 0L
subject <- gr_Gene
query <- bGR
hits <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query =
query,
subject = subject, maxgap = maxgap, minoverlap = minoverlap,
type = "any", select = "all", ignore.strand = T)))
bGR_genes <- unique(names(gr_Gene[hits[, 2]]))
if (verbose) {
now <- Sys.time()
message(sprintf("\t\%d nearby genes with \%d distance are defined as
the background (\%s) ...",
length(bGR_genes), gap.max, as.character(now)),
appendLF = T)
}
}
else {
bGR_genes <- NULL
}
if (verbose) {
now <- Sys.time()
message(sprintf("\n#######################################################",
appendLF = T))
message(sprintf("'xEnricherGenes' is being called (\%s):",
as.character(now)), appendLF = T)
message(sprintf("#######################################################",
appendLF = T))
}
eTerm <- xEnricherGenes(data = dGR_genes, background = bGR_genes,
ontology = ontology, size.range = size.range, min.overlap =
min.overlap,
which.distance = which.distance, test = test,
background.annotatable.only = background.annotatable.only,
p.tail = p.tail, p.adjust.method = p.adjust.method, ontology.algorithm
= ontology.algorithm,
elim.pvalue = elim.pvalue, lea.depth = lea.depth, path.mode =
path.mode,
true.path.rule = true.path.rule, verbose = verbose, RData.location =
RData.location,
guid = guid)
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("'xEnricherGenes' has been finished (\%s)!",
as.character(now)), appendLF = T)
message(sprintf("#######################################################\n",
appendLF = T))
}
endT <- Sys.time()
message(paste(c("\nEnd at ", as.character(endT)), collapse = ""),
appendLF = T)
runTime <- as.numeric(difftime(strptime(endT, "\%Y-\%m-\%d
\%H:\%M:\%S"),
strptime(startT, "\%Y-\%m-\%d \%H:\%M:\%S"), units = "secs"))
message(paste(c("Runtime in total is: ", runTime, " secs\n"),
collapse = ""), appendLF = T)
invisible(eTerm)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
