\name{xCircos}
\alias{xCircos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xCircos(g, entity = c("SNP", "Gene", "Both"), top_num = 50, colormap =
c("yr", "bwr", "jet", "gbr", "wyr", "br", "rainbow", "wb",
"lightyellow-orange"), rescale = T, nodes.query = NULL, ideogram = T,
chr.exclude = "auto", entity.label.cex = 0.7, entity.label.side =
c("in", "out"), entity.label.track = 1, entity.label.query = NULL,
GR.SNP = c("dbSNP_GWAS", "dbSNP_Common", "dbSNP_Single"), GR.Gene =
c("UCSC_knownGene", "UCSC_knownCanonical"), verbose = T, RData.location
= "http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{g}{
%% ~~Describe \code{g} here~~
}
\item{entity}{
%% ~~Describe \code{entity} here~~
}
\item{top_num}{
%% ~~Describe \code{top_num} here~~
}
\item{colormap}{
%% ~~Describe \code{colormap} here~~
}
\item{rescale}{
%% ~~Describe \code{rescale} here~~
}
\item{nodes.query}{
%% ~~Describe \code{nodes.query} here~~
}
\item{ideogram}{
%% ~~Describe \code{ideogram} here~~
}
\item{chr.exclude}{
%% ~~Describe \code{chr.exclude} here~~
}
\item{entity.label.cex}{
%% ~~Describe \code{entity.label.cex} here~~
}
\item{entity.label.side}{
%% ~~Describe \code{entity.label.side} here~~
}
\item{entity.label.track}{
%% ~~Describe \code{entity.label.track} here~~
}
\item{entity.label.query}{
%% ~~Describe \code{entity.label.query} here~~
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
function (g, entity = c("SNP", "Gene", "Both"), top_num = 50,
colormap = c("yr", "bwr", "jet", "gbr", "wyr", "br", "rainbow",
"wb", "lightyellow-orange"), rescale = T, nodes.query = NULL,
ideogram = T, chr.exclude = "auto", entity.label.cex = 0.7,
entity.label.side = c("in", "out"), entity.label.track = 1,
entity.label.query = NULL, GR.SNP = c("dbSNP_GWAS", "dbSNP_Common",
"dbSNP_Single"), GR.Gene = c("UCSC_knownGene", "UCSC_knownCanonical"),
verbose = T, RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
entity <- match.arg(entity)
entity.label.side <- match.arg(entity.label.side)
flag_package <- F
pkgs <- c("RCircos")
if (all(pkgs \%in\% rownames(utils::installed.packages()))) {
tmp <- sapply(pkgs, function(pkg) {
requireNamespace(pkg, quietly = T)
})
if (all(tmp)) {
flag_package <- T
}
}
if (!flag_package) {
stop("The package 'RCircos' is not available.\n")
}
if (class(g) != "igraph") {
stop("The function must apply to a 'igraph' object.\n")
}
df <- igraph::get.data.frame(g, what = "edges")
if (is.null(df$weight)) {
df$weight <- rep(1, nrow(df))
rescale <- F
}
else {
df$weight <- as.numeric(df$weight)
}
df <- df[with(df, order(-weight)), ]
if (!is.null(nodes.query)) {
ind_from <- which(!is.na(match(df[, 1], nodes.query)))
ind_to <- which(!is.na(match(df[, 2], nodes.query)))
ind <- union(ind_from, ind_to)
if (length(ind) > 0) {
df <- df[ind, ]
if (verbose) {
ind <- match(nodes.query, union(df[, 1], df[,
2]))
nodes.query <- nodes.query[!is.na(ind)]
now <- Sys.time()
message(sprintf("Circos plot restricted to nodes '\%s' (\%s) ...",
paste(nodes.query, collapse = ","), as.character(now)),
appendLF = T)
}
}
}
if (is.null(top_num)) {
top_num <- nrow(df)
}
if (top_num > nrow(df)) {
top_num <- nrow(df)
}
top_num <- as.integer(top_num)
df <- df[1:top_num, ]
if (verbose) {
now <- Sys.time()
message(sprintf("Loading positional information for \%s (\%s) ...",
entity, as.character(now)), appendLF = T)
}
if (entity == "SNP" | entity == "Both") {
if (class(GR.SNP) == "GRanges") {
pos_snp <- GR.SNP
}
else {
pos_snp <- xRDataLoader(RData.customised = GR.SNP[1],
verbose = verbose, RData.location = RData.location,
guid = guid)
if (is.null(pos_snp)) {
GR.SNP <- "dbSNP_GWAS"
if (verbose) {
message(sprintf("Instead, \%s will be used",
GR.SNP), appendLF = T)
}
pos_snp <- xRDataLoader(RData.customised = GR.SNP,
verbose = verbose, RData.location = RData.location,
guid = guid)
}
}
}
if (entity == "Gene" | entity == "Both") {
if (class(GR.Gene) == "GRanges") {
pos_gene <- GR.Gene
}
else {
pos_gene <- xRDataLoader(RData.customised = GR.Gene[1],
verbose = verbose, RData.location = RData.location,
guid = guid)
if (is.null(pos_gene)) {
GR.Gene <- "UCSC_knownGene"
if (verbose) {
message(sprintf("Instead, \%s will be used",
GR.Gene), appendLF = T)
}
pos_gene <- xRDataLoader(RData.customised = GR.Gene,
verbose = verbose, RData.location = RData.location,
guid = guid)
}
}
}
if (entity == "SNP") {
pos <- pos_snp
}
else if (entity == "Gene") {
pos <- pos_gene
}
else if (entity == "Both") {
GenomicRanges::mcols(pos_gene) <- NULL
GenomicRanges::mcols(pos_snp) <- NULL
pos <- c(pos_gene, pos_snp)
}
allnames <- names(pos)
A <- match(df$from, allnames)
B <- match(df$to, allnames)
flag <- !is.na(A) & !is.na(B)
AA <- A[flag]
BB <- B[flag]
input.data.A <- GenomicRanges::as.data.frame(pos[AA], row.names = NULL)
input.data.B <- GenomicRanges::as.data.frame(pos[BB], row.names = NULL)
input.data <- cbind.data.frame(input.data.A[, 1:3], input.data.B[,
1:3])
if (is.null(df$weight)) {
input.data$similarity <- rep(1, sum(flag))
}
else {
input.data$similarity <- as.numeric(as.character(df$weight[flag]))
}
label.data <- rbind(input.data.A[, 1:3], input.data.B[, 1:3])
label.data$Name <- c(df$from[flag], df$to[flag])
if (!is.null(chr.exclude)) {
chr.exclude <- chr.exclude[!is.na(chr.exclude)]
if (length(chr.exclude) == 0) {
chr.exclude <- NULL
}
else if (sum(chr.exclude == "auto") > 0) {
flag <- levels(label.data$seqnames) \%in\%
as.character(unique(label.data$seqnames))
chr.exclude <- levels(label.data$seqnames)[!flag]
}
}
if (verbose) {
now <- Sys.time()
message(sprintf("Loading human chromosome banding information (hg19)
(\%s) ...",
as.character(now)), appendLF = T)
}
eval(parse(text = "data(UCSC.HG19.Human.CytoBandIdeogram)"))
cyto.info <- ""
eval(parse(text = paste("cyto.info <-
UCSC.HG19.Human.CytoBandIdeogram",
sep = "")))
if (ideogram == F) {
cyto.info$Stain <- rep("gpos100", nrow(cyto.info))
}
if (verbose) {
now <- Sys.time()
message(sprintf("Initialising RCircos Core Components (\%s) ...",
as.character(now)), appendLF = T)
}
num.inside <- 1
num.outside <- 1
suppressMessages(RCircos.Set.Core.Components(cyto.info, chr.exclude,
num.inside, num.outside))
params <- RCircos.Get.Plot.Parameters()
if (0) {
params$track.padding <- 0
params$track.height <- 0.05
params$chr.ideog.pos <- 1
params$highlight.pos <- 1.1
params$chr.name.pos <- 1.1
params$track.out.start <- 1.2
params$highlight.width <- 0
}
params$text.size <- entity.label.cex
RCircos.Reset.Plot.Parameters(params)
if (verbose) {
now <- Sys.time()
message(sprintf("Plotting chromosome ideogram (\%s) ...",
as.character(now)), appendLF = T)
}
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
if (verbose) {
now <- Sys.time()
message(sprintf("Plotting link data (\%s) ...", as.character(now)),
appendLF = T)
}
if (rescale) {
sim <- input.data$similarity
if (verbose) {
now <- Sys.time()
message(sprintf("Also rescale similarity into the [0,1] range (\%s)",
as.character(now)), appendLF = T)
}
input.data$similarity <- (sim - min(sim))/(max(sim) -
min(sim))
}
palette.name <- xColormap(colormap = colormap)
cut_index <- as.numeric(cut(input.data$similarity, breaks =
min(input.data$similarity) +
(max(input.data$similarity) - min(input.data$similarity)) *
seq(0, 1, 0.05)))
cut_index[is.na(cut_index)] <- 1
input.data$PlotColor <- palette.name(20)[cut_index]
input.data <- input.data[order(input.data$similarity, decreasing = F),
]
RCircos.Link.Plot(input.data, track.num = 1, FALSE)
if (verbose) {
now <- Sys.time()
message(sprintf("Adding SNP and/or gene names (\%s) ...",
as.character(now)), appendLF = T)
}
name.col <- "Name"
side <- entity.label.side
track.num <- entity.label.track
label.data <- label.data[!duplicated(label.data$Name), ]
if (!is.null(entity.label.query)) {
ind <- match(label.data$Name, entity.label.query)
if (sum(!is.na(ind)) >= 1) {
label.data <- label.data[!is.na(ind), ]
}
}
if (verbose) {
RCircos.Gene.Name.Plot(label.data, name.col, track.num,
side)
}
else {
suppressMessages(RCircos.Gene.Name.Plot(label.data, name.col,
track.num, side))
}
invisible()
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
