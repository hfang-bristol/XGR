\name{xEnricherSNPs}
\alias{xEnricherSNPs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnricherSNPs(data, background = NULL, ontology = c("EF", "EF_disease",
"EF_phenotype", "EF_bp"), include.LD = NA, LD.r2 = 0.8, size.range =
c(10, 2000), min.overlap = 5, which.distance = NULL, test = c("fisher",
"hypergeo", "binomial"), background.annotatable.only = NULL, p.tail =
c("one-tail", "two-tails"), p.adjust.method = c("BH", "BY",
"bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm =
c("none", "pc", "elim", "lea"), elim.pvalue = 0.01, lea.depth = 2,
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
true.path.rule = T, verbose = T, silent = FALSE, RData.location =
"http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{background}{
%% ~~Describe \code{background} here~~
}
\item{ontology}{
%% ~~Describe \code{ontology} here~~
}
\item{include.LD}{
%% ~~Describe \code{include.LD} here~~
}
\item{LD.r2}{
%% ~~Describe \code{LD.r2} here~~
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
function (data, background = NULL, ontology = c("EF", "EF_disease",
"EF_phenotype", "EF_bp"), include.LD = NA, LD.r2 = 0.8, size.range =
c(10,
2000), min.overlap = 5, which.distance = NULL, test = c("fisher",
"hypergeo", "binomial"), background.annotatable.only = NULL,
p.tail = c("one-tail", "two-tails"), p.adjust.method = c("BH",
"BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm =
c("none",
"pc", "elim", "lea"), elim.pvalue = 0.01, lea.depth = 2,
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
true.path.rule = T, verbose = T, silent = FALSE, RData.location =
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
ontology <- match.arg(ontology)
test <- match.arg(test)
p.tail <- match.arg(p.tail)
p.adjust.method <- match.arg(p.adjust.method)
ontology.algorithm <- match.arg(ontology.algorithm)
path.mode <- match.arg(path.mode)
p.tail <- match.arg(p.tail)
if (length(data) == 0) {
return(NULL)
}
if (is.vector(data)) {
data <- unique(data)
}
else {
stop("The input data must be a vector.\n")
}
data <- as.character(data)
if (!is.na(ontology)) {
if (verbose) {
now <- Sys.time()
message(sprintf("Load the ontology \%s and its SNP annotations (\%s)
...",
ontology, as.character(now)), appendLF = T)
}
ig <- xRDataLoader(RData = paste("ig.EF", sep = ""),
RData.location = RData.location, guid = guid, verbose = verbose)
V(ig)$namespace <- ontology
if (ontology != "EF") {
if (ontology == "EF_disease") {
node <- "EFO:0000408"
}
else if (ontology == "EF_phenotype") {
node <- "EFO:0000651"
}
else if (ontology == "EF_bp") {
node <- "GO:0008150"
}
neighs.out <- igraph::neighborhood(ig, order = vcount(ig),
nodes = node, mode = "out")
nodeInduced <- V(ig)[unique(unlist(neighs.out))]$name
g <- igraph::induced.subgraph(ig, vids = nodeInduced)
V(g)$namespace <- ontology
}
else {
g <- ig
}
if (is.null(V(g)$term_namespace)) {
V(g)$term_namespace <- ontology
}
anno <- xRDataLoader(RData = paste("GWAS2EF", sep = ""),
RData.location = RData.location, guid = guid, verbose = verbose)
if (LD.r2 >= 0.8 & LD.r2 <= 1) {
default.include.LD <- c("ACB", "AFR", "AMR", "ASW",
"BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "EAS",
"ESN", "EUR", "FIN", "GBR", "GIH", "GWD", "IBS",
"ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL",
"PJL", "PUR", "SAS", "STU", "TSI", "YRI")
ind <- match(default.include.LD, include.LD)
include.LD <- default.include.LD[!is.na(ind)]
if (length(include.LD) > 0) {
GWAS_LD <- xRDataLoader(RData = "GWAS_LD", RData.location =
RData.location,
guid = guid, verbose = verbose)
ld.list <- lapply(include.LD, function(x) {
data_ld <- ""
eval(parse(text = paste("data_ld <- GWAS_LD$",
x, sep = "")))
ind <- match(colnames(data_ld), data)
ld <- rownames(data_ld)[which(Matrix::rowSums(data_ld[,
!is.na(ind)] >= LD.r2) > 0)]
ld
})
data <- base::union(unlist(ld.list), data)
}
}
}
else {
stop("There is no input for the ontology.\n")
}
if (verbose) {
now <- Sys.time()
message(sprintf("\n#######################################################",
appendLF = T))
message(sprintf("'xEnricher' is being called (\%s):",
as.character(now)), appendLF = T)
message(sprintf("#######################################################",
appendLF = T))
}
eTerm <- xEnricher(data = data, annotation = anno, g = g,
background = background, size.range = size.range, min.overlap =
min.overlap,
which.distance = which.distance, test = test,
background.annotatable.only = background.annotatable.only,
p.tail = p.tail, p.adjust.method = p.adjust.method, ontology.algorithm
= ontology.algorithm,
elim.pvalue = elim.pvalue, lea.depth = lea.depth, path.mode =
path.mode,
true.path.rule = true.path.rule, verbose = verbose)
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("'xEnricher' has been finished (\%s)!",
as.character(now)), appendLF = T)
message(sprintf("#######################################################\n",
appendLF = T))
}
endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, "\%Y-\%m-\%d
\%H:\%M:\%S"),
strptime(startT, "\%Y-\%m-\%d \%H:\%M:\%S"), units = "secs"))
if (!silent) {
message(paste(c("\nEnd at ", as.character(endT)), collapse = ""),
appendLF = TRUE)
message(paste(c("Runtime in total (xEnricherSNPs): ",
runTime, " secs\n"), collapse = ""), appendLF = TRUE)
}
invisible(eTerm)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
