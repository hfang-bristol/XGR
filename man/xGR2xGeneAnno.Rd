\name{xGR2xGeneAnno}
\alias{xGR2xGeneAnno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGR2xGeneAnno(data, background = NULL, format = c("chr:start-end",
"data.frame", "bed", "GRanges"), build.conversion = c(NA,
"hg38.to.hg19", "hg18.to.hg19"), crosslink = c("genehancer",
"PCHiC_combined", "GTEx_V6p_combined", "nearby"), crosslink.customised
= NULL, crosslink.top = NULL, nearby.distance.max = 50000,
nearby.decay.kernel = c("rapid", "slow", "linear", "constant"),
nearby.decay.exponent = 2, ontology = NA, size.range = c(10, 2000),
min.overlap = 5, which.distance = NULL, test = c("hypergeo", "fisher",
"binomial"), background.annotatable.only = NULL, p.tail = c("one-tail",
"two-tails"), p.adjust.method = c("BH", "BY", "bonferroni", "holm",
"hochberg", "hommel"), ontology.algorithm = c("none", "pc", "elim",
"lea"), elim.pvalue = 0.01, lea.depth = 2, path.mode = c("all_paths",
"shortest_paths", "all_shortest_paths"), true.path.rule = F,
out.evidence = T, out.evidence.plot = F, verbose = T, silent = F,
RData.location = "http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{background}{
%% ~~Describe \code{background} here~~
}
\item{format}{
%% ~~Describe \code{format} here~~
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
\item{crosslink.top}{
%% ~~Describe \code{crosslink.top} here~~
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
\item{out.evidence}{
%% ~~Describe \code{out.evidence} here~~
}
\item{out.evidence.plot}{
%% ~~Describe \code{out.evidence.plot} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
}
\item{silent}{
%% ~~Describe \code{silent} here~~
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
function (data, background = NULL, format = c("chr:start-end",
"data.frame", "bed", "GRanges"), build.conversion = c(NA,
"hg38.to.hg19", "hg18.to.hg19"), crosslink = c("genehancer",
"PCHiC_combined", "GTEx_V6p_combined", "nearby"), crosslink.customised
= NULL,
crosslink.top = NULL, nearby.distance.max = 50000, nearby.decay.kernel
= c("rapid",
"slow", "linear", "constant"), nearby.decay.exponent = 2,
ontology = NA, size.range = c(10, 2000), min.overlap = 5,
which.distance = NULL, test = c("hypergeo", "fisher", "binomial"),
background.annotatable.only = NULL, p.tail = c("one-tail",
"two-tails"), p.adjust.method = c("BH", "BY", "bonferroni",
"holm", "hochberg", "hommel"), ontology.algorithm = c("none",
"pc", "elim", "lea"), elim.pvalue = 0.01, lea.depth = 2,
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
true.path.rule = F, out.evidence = T, out.evidence.plot = F,
verbose = T, silent = F, RData.location =
"http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
startT <- Sys.time()
if (!silent) {
message(paste(c("Start at ", as.character(startT)), collapse = ""),
appendLF = TRUE)
message("", appendLF = TRUE)
}
else {
verbose <- FALSE
}
format <- match.arg(format)
build.conversion <- match.arg(build.conversion)
nearby.decay.kernel <- match.arg(nearby.decay.kernel)
test <- match.arg(test)
p.tail <- match.arg(p.tail)
p.adjust.method <- match.arg(p.adjust.method)
ontology.algorithm <- match.arg(ontology.algorithm)
path.mode <- match.arg(path.mode)
if (verbose) {
now <- Sys.time()
message(sprintf("First, import the data/background formatted as '\%s'
(\%s) ...",
format, as.character(now)), appendLF = T)
}
dGR <- xGR(data = data, format = format, build.conversion =
build.conversion,
verbose = verbose, RData.location = RData.location, guid = guid)
bGR <- xGR(data = background, format = format, build.conversion =
build.conversion,
verbose = verbose, RData.location = RData.location, guid = guid)
if (verbose) {
now <- Sys.time()
message(sprintf("Second, define crosslinked genes based on '\%s' (\%s)
...",
crosslink, as.character(now)), appendLF = T)
}
df_xGenes_data <- xGR2xGenes(data = dGR, format = "GRanges",
crosslink = crosslink, crosslink.customised = crosslink.customised,
cdf.function = "original", scoring = TRUE, scoring.scheme = "max",
scoring.rescale = F, nearby.distance.max = nearby.distance.max,
nearby.decay.kernel = nearby.decay.kernel, nearby.decay.exponent =
nearby.decay.exponent,
verbose = verbose, silent = !verbose, RData.location = RData.location,
guid = guid)
df_xGenes_background <- xGR2xGenes(data = bGR, format = "GRanges",
crosslink = crosslink, crosslink.customised = crosslink.customised,
cdf.function = "original", scoring = TRUE, scoring.scheme = "max",
scoring.rescale = F, nearby.distance.max = nearby.distance.max,
nearby.decay.kernel = nearby.decay.kernel, nearby.decay.exponent =
nearby.decay.exponent,
verbose = verbose, silent = !verbose, RData.location = RData.location,
guid = guid)
Score <- Gene <- NULL
df_xGenes_data <- df_xGenes_data \%>\% dplyr::arrange(-Score)
if (is.null(crosslink.top)) {
crosslink.top <- nrow(df_xGenes_data)
}
if (crosslink.top > nrow(df_xGenes_data)) {
crosslink.top <- nrow(df_xGenes_data)
}
crosslink.top <- as.integer(crosslink.top)
crosslink.cutoff <- df_xGenes_data[crosslink.top, "Score"]
dGR_genes <- df_xGenes_data$Gene[df_xGenes_data$Score >=
crosslink.cutoff]
if (!is.null(df_xGenes_background)) {
bGR_genes <- (df_xGenes_background \%>\% dplyr::arrange(-Score))$Gene
}
else {
bGR_genes <- NULL
}
if (verbose) {
if (is.null(bGR_genes)) {
message(sprintf("\t\%d (out of \%d crosslinked genes) are used.",
length(dGR_genes), nrow(df_xGenes_data), as.character(Sys.time())),
appendLF = T)
}
else {
message(sprintf("\t\%d (out of \%d crosslinked genes) and \%d
background genes are used.",
length(dGR_genes), nrow(df_xGenes_data), length(bGR_genes),
as.character(Sys.time())), appendLF = T)
}
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
true.path.rule = true.path.rule, verbose = verbose, silent = !verbose,
RData.location = RData.location, guid = guid)
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("'xEnricherGenes' has been finished (\%s)!",
as.character(now)), appendLF = T)
message(sprintf("#######################################################\n",
appendLF = T))
}
if (!is.null(eTerm)) {
ind <- match(df_xGenes_data$Gene, dGR_genes)
eTerm$crosslink <- df_xGenes_data[!is.na(ind), c("Gene",
"Score", "Pval")]
if (out.evidence) {
df_evidence <- xGR2xGenes(data = dGR, format = "GRanges",
crosslink = crosslink, crosslink.customised = crosslink.customised,
cdf.function = "original", scoring = FALSE, scoring.scheme = "max",
scoring.rescale = F, nearby.distance.max = nearby.distance.max,
nearby.decay.kernel = nearby.decay.kernel, nearby.decay.exponent =
nearby.decay.exponent,
verbose = verbose, silent = !verbose, RData.location = RData.location,
guid = guid)
ind <- match(df_evidence$Gene, dGR_genes)
evidence <- df_evidence[!is.na(ind), c("GR", "Gene",
"Score")]
eTerm$evidence <- evidence
if (out.evidence.plot) {
Gene <- Score <- NULL
mat_evidence <- tidyr::spread(evidence, key = Gene,
value = Score)
mat <- mat_evidence[, -1]
rownames(mat) <- mat_evidence[, 1]
ind <- xGRsort(rownames(mat))
mat <- mat[ind, ]
rowsep <- xGRsep(rownames(mat))
rowsep <- nrow(mat) - rowsep
if (ncol(mat) >= 0) {
reorder <- "none"
}
else {
reorder <- "col"
}
gp_evidence <- xHeatmap(mat, reorder = reorder,
colormap = "spectral", ncolors = 64, barwidth = 0.4,
x.rotate = 90, shape = 19, size = 2, x.text.size = 6,
y.text.size = 6, na.color = "transparent")
gp_evidence <- gp_evidence + theme(legend.title = element_text(size =
8),
legend.position = "left") + scale_y_discrete(position = "right")
gp_evidence <- gp_evidence + geom_hline(yintercept = rowsep +
0.5, color = "grey90", size = 0.5)
eTerm$gp_evidence <- gp_evidence
}
}
}
endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, "\%Y-\%m-\%d
\%H:\%M:\%S"),
strptime(startT, "\%Y-\%m-\%d \%H:\%M:\%S"), units = "secs"))
if (!silent) {
message(paste(c("\nEnd at ", as.character(endT)), collapse = ""),
appendLF = TRUE)
message(paste(c("Runtime in total (xGR2xGeneAnno): ",
runTime, " secs\n"), collapse = ""), appendLF = TRUE)
}
invisible(eTerm)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
