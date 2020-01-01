\name{xGRviaGenomicAnnoAdv}
\alias{xGRviaGenomicAnnoAdv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGRviaGenomicAnnoAdv(data.file, annotation.file = NULL, background.file
= NULL, format.file = c("data.frame", "bed", "chr:start-end",
"GRanges"), build.conversion = c(NA, "hg38.to.hg19", "hg18.to.hg19"),
background.annotatable.only = F, num.samples = 1000, gap.max = 50000,
max.distance = NULL, p.adjust.method = c("BH", "BY", "bonferroni",
"holm", "hochberg", "hommel"), GR.annotation = NA, parallel = TRUE,
multicores = NULL, verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
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
\item{background.annotatable.only}{
%% ~~Describe \code{background.annotatable.only} here~~
}
\item{num.samples}{
%% ~~Describe \code{num.samples} here~~
}
\item{gap.max}{
%% ~~Describe \code{gap.max} here~~
}
\item{max.distance}{
%% ~~Describe \code{max.distance} here~~
}
\item{p.adjust.method}{
%% ~~Describe \code{p.adjust.method} here~~
}
\item{GR.annotation}{
%% ~~Describe \code{GR.annotation} here~~
}
\item{parallel}{
%% ~~Describe \code{parallel} here~~
}
\item{multicores}{
%% ~~Describe \code{multicores} here~~
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
background.annotatable.only = F, num.samples = 1000, gap.max = 50000,
max.distance = NULL, p.adjust.method = c("BH", "BY", "bonferroni",
"holm", "hochberg", "hommel"), GR.annotation = NA, parallel = TRUE,
multicores = NULL, verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
startT <- Sys.time()
message(paste(c("Start at ", as.character(startT)), collapse = ""),
appendLF = T)
message("", appendLF = T)
format.file <- match.arg(format.file)
build.conversion <- match.arg(build.conversion)
p.adjust.method <- match.arg(p.adjust.method)
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
message("\t\tThe file 'annotation.file' is not provided, so built-in
RData will be used instead!")
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
2] + 1), end = as.numeric(background[, 3])),
strand = S4Vectors::Rle(rep("*", nrow(data))))
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
mergeOverlaps <- function(qGR, sGR, out.format = c("GR",
"counts"), maxgap = -1L, minoverlap = 0L) {
out.format <- match.arg(out.format)
hits <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query =
qGR,
subject = sGR, maxgap = maxgap, minoverlap = minoverlap,
type = "any", select = "all", ignore.strand = T)))
qhits <- qGR[hits[, 1]]
shits <- sGR[hits[, 2]]
gr <- IRanges::pintersect(qhits, shits, ignore.strand = T)
if (out.format == "GR") {
gr
}
else {
sum(as.numeric(IRanges::width(gr)))
}
}
mergeOverlaps_GRL <- function(qGR, sGRL, maxgap = -1L, minoverlap = 0L,
out.format = c("GR", "counts")) {
out.format <- match.arg(out.format)
hits_GRL <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query =
qGR,
subject = sGRL, maxgap = maxgap, minoverlap = minoverlap,
type = "any", select = "all", ignore.strand = T)))
qhits_GRL <- hits_GRL[, 1]
shits_GRL <- hits_GRL[, 2]
res_ls <- split(x = qhits_GRL, f = shits_GRL)
names_ls <- as.numeric(names(res_ls))
if (1) {
res <- lapply(1:length(res_ls), function(i) {
s_ind <- names_ls[i]
sGR <- sGRL[[s_ind]]
mergeOverlaps(qGR = qGR, sGR = sGR, out.format = out.format,
maxgap = maxgap, minoverlap = minoverlap)
})
}
else {
sGRL_o <- sGRL[names_ls]
res <- lapply(sGRL_o, function(x) {
mergeOverlaps(qGR = qGR, sGR = x, out.format = out.format,
maxgap = maxgap, minoverlap = minoverlap)
})
}
if (out.format == "GR") {
out <- vector("list", length(sGRL))
names(out) <- names(sGRL)
out[names_ls] <- res
out <- lapply(out, function(x) {
if (is.null(x)) {
x <- GenomicRanges::GRanges()
}
else {
x
}
})
}
else {
out <- as.list(rep(0, length(sGRL)))
names(out) <- names(sGRL)
out[names_ls] <- res
}
out
}
progress_indicate <- function(i, B, step, flag = F) {
if (i\%\%ceiling(B/step) == 0 | i == B | i == 1) {
if (flag & verbose) {
message(sprintf("\t\%d out of \%d (\%s)", i, B,
as.character(Sys.time())), appendLF = T)
}
}
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
aGRL <- GenomicRanges::GRangesList(aGR_reduced)
oGR_reduced <- mergeOverlaps_GRL(qGR = dGR_reduced, sGRL = aGRL,
maxgap = -1L, minoverlap = 0L)
if (verbose) {
now <- Sys.time()
message(sprintf("Forth, generate null distribution via doing \%d
sampling (\%s) ...",
num.samples, as.character(now)), appendLF = T)
}
if (verbose) {
now <- Sys.time()
message(sprintf("\n#######################################################",
appendLF = T))
message(sprintf("'xGRsampling' is being called (\%s):",
as.character(now)), appendLF = T)
message(sprintf("#######################################################",
appendLF = T))
}
sGR_list <- xGRsampling(GR.data = dGR_reduced, GR.background =
bGR_reduced,
num.samples = num.samples, gap.max = gap.max, max.distance =
max.distance,
verbose = verbose)
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("'xGRsampling' has been finished (\%s)!",
as.character(now)), appendLF = T)
message(sprintf("#######################################################\n",
appendLF = T))
}
if (verbose) {
now <- Sys.time()
message(sprintf("Fifth, perform enrichment analysis (\%s) ...",
as.character(now)), appendLF = T)
}
data_nBases <- sum(as.numeric(IRanges::width(dGR_reduced)))
annotation_nBases <- base::sapply(aGR_reduced, function(gr) {
sum(as.numeric(IRanges::width(gr)))
})
background_nBases <- sum(as.numeric(IRanges::width(bGR_reduced)))
overlap_nBases <- base::sapply(oGR_reduced, function(gr) {
sum(as.numeric(IRanges::width(gr)))
})
if (verbose) {
now <- Sys.time()
message(sprintf("\tthe number of nucleotides/bases: data (\%d),
background (\%d)",
data_nBases, background_nBases), appendLF = T)
message(sprintf("\tthe number of annotations: \%d",
length(annotation_nBases)),
appendLF = T)
}
flag_parallel <- F
if (parallel == TRUE) {
flag_parallel <- xCheckParallel(multicores = multicores,
verbose = verbose)
if (flag_parallel) {
i <- 1
b2f <- foreach::`\%dopar\%`(foreach::foreach(i = 1:length(sGR_list),
.inorder = T, .combine = rbind), {
progress_indicate(i, length(sGR_list), 10, flag = T)
res <- mergeOverlaps_GRL(qGR = sGR_list[[i]],
sGRL = aGRL, out.format = "counts")
unlist(res)
})
}
}
if (flag_parallel == F) {
b_overlap_nBases <- base::lapply(1:length(sGR_list),
function(i) {
progress_indicate(i, length(sGR_list), 10, flag = T)
res <- mergeOverlaps_GRL(qGR = sGR_list[[i]],
sGRL = aGRL, out.format = "counts")
unlist(res)
})
b2f <- do.call(rbind, b_overlap_nBases)
}
obs <- overlap_nBases
exp_mean <- apply(b2f, 2, mean)
exp_std <- apply(b2f, 2, stats::sd)
ratio <- obs/exp_mean
zscore <- (obs - exp_mean)/exp_std
obs_matrix <- matrix(rep(obs, each = num.samples), nrow = num.samples)
pvalue <- apply((obs_matrix - b2f) <= 0, 2, sum)/num.samples
zscore[is.na(zscore)] <- 0
zscore[is.infinite(zscore)] <- max(zscore[!is.infinite(zscore)])
pvalue[is.na(ratio)] <- 1
ratio[is.na(ratio)] <- 1
enrichment_df <- data.frame(names(overlap_nBases), annotation_nBases,
data_nBases, background_nBases, obs, exp_mean, ratio,
zscore, pvalue, row.names = NULL, stringsAsFactors = F)
colnames(enrichment_df) <- c("name", "nAnno", "nData", "nBG",
"nOverlap", "nExpect", "fc", "zscore", "pvalue")
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
"fc", "zscore", "pvalue", "adjp", "nData", "nBG")]
invisible(res_df)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
