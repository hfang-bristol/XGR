\name{xGRviaGenomicAnno}
\alias{xGRviaGenomicAnno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGRviaGenomicAnno(data.file, annotation.file = NULL, background.file =
NULL, format.file = c("data.frame", "bed", "chr:start-end", "GRanges"),
build.conversion = c(NA, "hg38.to.hg19", "hg18.to.hg19"), resolution =
c("bases", "regions", "hybrid"), background.annotatable.only = T,
p.tail = c("one-tail", "two-tails"), p.adjust.method = c("BH", "BY",
"bonferroni", "holm", "hochberg", "hommel"), GR.annotation = NA,
verbose = T, RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data.file}{
%% ~~Describe \code{data.file} here~~
}
\item{annotation.file}{
%% ~~Describe \code{annotation.file} here~~
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
\item{resolution}{
%% ~~Describe \code{resolution} here~~
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
\item{GR.annotation}{
%% ~~Describe \code{GR.annotation} here~~
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
function (data.file, annotation.file = NULL, background.file = NULL,
format.file = c("data.frame", "bed", "chr:start-end", "GRanges"),
build.conversion = c(NA, "hg38.to.hg19", "hg18.to.hg19"),
resolution = c("bases", "regions", "hybrid"),
background.annotatable.only = T,
p.tail = c("one-tail", "two-tails"), p.adjust.method = c("BH",
"BY", "bonferroni", "holm", "hochberg", "hommel"), GR.annotation = NA,
verbose = T, RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
startT <- Sys.time()
message(paste(c("Start at ", as.character(startT)), collapse = ""),
appendLF = T)
message("", appendLF = T)
format.file <- match.arg(format.file)
build.conversion <- match.arg(build.conversion)
resolution <- match.arg(resolution)
p.adjust.method <- match.arg(p.adjust.method)
p.tail <- match.arg(p.tail)
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
message(sprintf("\timport the annotation file (\%s) ...",
as.character(now)), appendLF = T)
}
if (is.matrix(annotation.file) | is.data.frame(annotation.file) |
class(annotation.file) == "list") {
annotation <- annotation.file
}
else if (!is.null(annotation.file)) {
annotation <- utils::read.delim(file = annotation.file,
header = F, row.names = NULL, stringsAsFactors = F)
}
else {
if (verbose) {
message("\t\tThe file 'annotation.file' is not provided, so built-in
RData will be used instead!")
}
annotation <- NULL
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
if (!is.null(annotation)) {
if (ncol(annotation) >= 4) {
annotation <- annotation
}
else {
stop("Your input 'annotation.file' is not as expected!\n")
}
anno_ls <- split(x = annotation[, -4], f = annotation[,
4])
aGR <- base::lapply(anno_ls, function(x) {
ind <- suppressWarnings(which(!is.na(as.numeric(x[,
2])) & !is.na(as.numeric(x[, 3]))))
x <- x[ind, ]
gr <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(x[,
1]), ranges = IRanges::IRanges(start = as.numeric(x[,
2]), end = as.numeric(x[, 3])), strand = S4Vectors::Rle(rep("*",
nrow(x))))
})
}
else {
aGRL <- xDefineGenomicAnno(GR.annotation, verbose = verbose,
RData.location = RData.location, guid = guid)
aGR <- lapply(aGRL, function(x) x)
}
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
input <- do.call(rbind, strsplit(data[, 1], ":|-"))
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
if (!is.null(annotation)) {
input <- do.call(rbind, strsplit(annotation[, 1],
":|-"))
if (ncol(input) >= 3) {
annotation <- cbind(input[, 1:3], annotation[,
2])
}
else if (ncol(input) == 2) {
annotation <- cbind(input[, c(1, 2, 2)], annotation[,
2])
}
else {
stop("Your input 'annotation.file' does not meet the format
'chr:start-end'!\n")
}
anno_ls <- split(x = annotation[, -4], f = annotation[,
4])
aGR <- base::lapply(anno_ls, function(x) {
ind <- suppressWarnings(which(!is.na(as.numeric(x[,
2])) & !is.na(as.numeric(x[, 3]))))
x <- x[ind, ]
gr <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(x[,
1]), ranges = IRanges::IRanges(start = as.numeric(x[,
2]), end = as.numeric(x[, 3])), strand = S4Vectors::Rle(rep("*",
nrow(x))))
})
}
else {
aGRL <- xDefineGenomicAnno(GR.annotation, verbose = verbose,
RData.location = RData.location, guid = guid)
aGR <- lapply(aGRL, function(x) x)
}
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
if (!is.null(annotation)) {
anno_ls <- split(x = annotation[, -4], f = annotation[,
4])
aGR <- base::lapply(anno_ls, function(x) {
ind <- suppressWarnings(which(!is.na(as.numeric(x[,
2])) & !is.na(as.numeric(x[, 3]))))
x <- x[ind, ]
gr <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(x[,
1]), ranges = IRanges::IRanges(start = as.numeric(x[,
2]) + 1, end = as.numeric(x[, 3])), strand = S4Vectors::Rle(rep("*",
nrow(x))))
})
}
else {
aGRL <- xDefineGenomicAnno(GR.annotation, verbose = verbose,
RData.location = RData.location, guid = guid)
aGR <- lapply(aGRL, function(x) x)
}
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
if (!is.null(annotation)) {
aGR <- annotation
}
else {
aGRL <- xDefineGenomicAnno(GR.annotation, verbose = verbose,
RData.location = RData.location, guid = guid)
aGR <- lapply(aGRL, function(x) x)
}
if (!is.null(background)) {
bGR <- background
}
else {
bGR <- NULL
}
}
mergeOverlaps <- function(qGR, sGR, maxgap = -1L, minoverlap = 0L) {
hits <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query =
qGR,
subject = sGR, maxgap = maxgap, minoverlap = minoverlap,
type = "any", select = "all", ignore.strand = T)))
qhits <- qGR[hits[, 1]]
shits <- sGR[hits[, 2]]
oGR <- IRanges::pintersect(qhits, shits, ignore.strand = T)
IRanges::reduce(oGR)
}
doBinomialTest <- function(X, K, M, N, p.tail) {
N <- max(N, M)
if (K == 0 || M == 0 || N == 0) {
p.value <- 1
}
else {
if (p.tail == "one-tail") {
p.value <- stats::pbinom(X, K, M/N, lower.tail = F,
log.p = F)
}
else {
if (X >= K * M/N) {
p.value <- stats::pbinom(X, K, M/N, lower.tail = F,
log.p = F)
}
else {
p.value <- stats::pbinom(X, K, M/N, lower.tail = T,
log.p = F)
}
}
}
return(p.value)
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
if (!is.null(annotation.file)) {
if (verbose) {
message(sprintf("\tannotation genomic regions: lifted over via genome
build conversion `\%s`",
build.conversion), appendLF = T)
}
aGR <- lapply(aGR, function(gr) {
xLiftOver(data.file = gr, format.file = "GRanges",
build.conversion = build.conversion, merged = F,
verbose = verbose, RData.location = RData.location,
guid = guid)
})
}
if (!is.null(bGR)) {
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
}
if (verbose) {
now <- Sys.time()
message(sprintf("Third, define the background (\%s) ...",
as.character(now)), appendLF = T)
}
dGR_reduced <- IRanges::reduce(dGR)
aGR_reduced <- base::lapply(aGR, function(x) {
IRanges::reduce(x)
})
if (is.null(bGR)) {
if (verbose) {
now <- Sys.time()
message(sprintf("\tall annotatable regions (by default) are used as the
background (\%s) ...",
as.character(now)), appendLF = T)
}
aGRL <- GenomicRanges::GRangesList(aGR_reduced)
bGR_reduced <- IRanges::reduce(BiocGenerics::unlist(aGRL))
}
else {
bGR_reduced <- IRanges::reduce(bGR)
aGR_reduced <- base::lapply(aGR_reduced, function(gr) {
mergeOverlaps(qGR = gr, sGR = bGR_reduced, maxgap = -1L,
minoverlap = 0L)
})
if (background.annotatable.only) {
if (verbose) {
now <- Sys.time()
message(sprintf("\tthe given background regions but restricted to the
annotatable are used as the background (\%s) ...",
as.character(now)), appendLF = T)
}
aGRL <- GenomicRanges::GRangesList(aGR_reduced)
bGR_reduced <- IRanges::reduce(BiocGenerics::unlist(aGRL))
}
else {
if (verbose) {
now <- Sys.time()
message(sprintf("\tthe given background regions are used as the
background (\%s) ...",
as.character(now)), appendLF = T)
}
}
}
dGR_reduced <- mergeOverlaps(qGR = dGR_reduced, sGR = bGR_reduced,
maxgap = -1L, minoverlap = 0L)
oGR_reduced <- base::lapply(aGR_reduced, function(gr) {
mergeOverlaps(qGR = gr, sGR = dGR_reduced, maxgap = -1L,
minoverlap = 0L)
})
if (verbose) {
now <- Sys.time()
message(sprintf("Forth, perform enrichment analysis at '\%s' resolution
with '\%s' p-values (\%s) ...",
resolution, p.tail, as.character(now)), appendLF = T)
}
if (resolution == "bases") {
data_nBases <- sum(as.numeric(IRanges::width(dGR_reduced)))
overlap_nBases <- base::sapply(oGR_reduced, function(gr) {
sum(as.numeric(IRanges::width(gr)))
})
annotation_nBases <- base::sapply(aGR_reduced, function(gr) {
sum(as.numeric(IRanges::width(gr)))
})
background_nBases <- sum(as.numeric(IRanges::width(bGR_reduced)))
if (verbose) {
now <- Sys.time()
message(sprintf("\tthe number of bases: data (\%d)",
data_nBases), appendLF = T)
message(sprintf("\tthe number of annotations: \%d",
length(annotation_nBases)), appendLF = T)
}
}
else if (resolution == "regions") {
data_nBases <- length(dGR_reduced)
overlap_nBases <- base::sapply(oGR_reduced, length)
annotation_nBases <- base::sapply(aGR_reduced, length)
background_nBases <- max(length(bGR_reduced), max(annotation_nBases))
if (verbose) {
now <- Sys.time()
message(sprintf("\tthe number of regions: data (\%d)",
data_nBases), appendLF = T)
message(sprintf("\tthe number of annotations: \%d",
length(annotation_nBases)), appendLF = T)
}
}
else if (resolution == "hybrid") {
data_nBases <- length(dGR_reduced)
overlap_nBases <- base::sapply(oGR_reduced, length)
annotation_nBases <- base::sapply(aGR_reduced, function(gr) {
sum(as.numeric(IRanges::width(gr)))
})
background_nBases <- sum(as.numeric(IRanges::width(bGR_reduced)))
if (verbose) {
now <- Sys.time()
message(sprintf("\tthe number of regions: data (\%d)",
data_nBases), appendLF = T)
message(sprintf("\tthe number of annotations: \%d",
length(annotation_nBases)), appendLF = T)
}
}
res_ls <- base::lapply(1:length(overlap_nBases), function(i) {
X <- as.numeric(overlap_nBases[i])
K <- data_nBases
M <- as.numeric(annotation_nBases[i])
N <- background_nBases
x.exp <- K * M/N
var.exp <- K * M/N * (N - M)/N * (N - K)/(N - 1)
if (is.na(var.exp)) {
z.score <- 0
}
else {
if (var.exp != 0) {
suppressWarnings(z.score <- (X - x.exp)/sqrt(var.exp))
}
else {
z.score <- 0
}
}
if (is.na(z.score)) {
z.score <- 0
}
p.value <- doBinomialTest(X, K, M, N, p.tail)
if (all(resolution \%in\% c("bases", "regions"))) {
cTab <- matrix(c(X, K - X, M - X, N - M - K + X),
nrow = 2, dimnames = list(c("anno", "notAnno"),
c("group", "notGroup")))
if (class(suppressWarnings(try(res <- stats::fisher.test(cTab),
T))) == "try-error") {
or <- CIl <- CIu <- NA
}
else {
res <- stats::fisher.test(cTab)
or <- as.vector(res$estimate)
CIl <- as.vector(res$conf.int)[1]
CIu <- as.vector(res$conf.int)[2]
}
}
else {
or <- CIl <- CIu <- NA
}
c(X, K, M, N, X/K, M/N, (X/K)/(M/N), z.score, p.value,
or, CIl, CIu)
})
res_df <- do.call(rbind, res_ls)
enrichment_df <- data.frame(names(overlap_nBases), res_df,
stringsAsFactors = F)
colnames(enrichment_df) <- c("name", "nOverlap", "nData",
"nAnno", "nBG", "obsProb", "expProb", "fc", "zscore",
"pvalue", "or", "CIl", "CIu")
p.adjust.method = c("BH", "BY", "bonferroni", "holm", "hochberg",
"hommel")[1]
pvals <- enrichment_df$pvalue
adjpvals <- stats::p.adjust(pvals, method = p.adjust.method)
enrichment_df$adjp <- adjpvals
enrichment_df$zscore <- signif(enrichment_df$zscore, digits = 3)
pvals <- enrichment_df$pvalue
adjpvals <- enrichment_df$adjp
pvals <- signif(pvals, digits = 2)
adjpvals <- signif(adjpvals, digits = 2)
pvals <- base::sapply(pvals, function(x) {
if (x < 0.1 & x != 0) {
as.numeric(format(x, scientific = T))
}
else {
x
}
})
adjpvals <- base::sapply(adjpvals, function(x) {
if (x < 0.1 & x != 0) {
as.numeric(format(x, scientific = T))
}
else {
x
}
})
enrichment_df$pvalue <- pvals
enrichment_df$adjp <- adjpvals
endT <- Sys.time()
message(paste(c("\nEnd at ", as.character(endT)), collapse = ""),
appendLF = T)
runTime <- as.numeric(difftime(strptime(endT, "\%Y-\%m-\%d
\%H:\%M:\%S"),
strptime(startT, "\%Y-\%m-\%d \%H:\%M:\%S"), units = "secs"))
message(paste(c("Runtime in total is: ", runTime, " secs\n"),
collapse = ""), appendLF = T)
res_df <- enrichment_df[, c("name", "nAnno", "nOverlap",
"fc", "zscore", "pvalue", "adjp", "or", "CIl", "CIu",
"expProb", "obsProb")]
eTerm <- xEnrichViewer(res_df, top_num = nrow(res_df), sortBy =
"zscore")
invisible(eTerm)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
