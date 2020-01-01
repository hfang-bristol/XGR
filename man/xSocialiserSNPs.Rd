\name{xSocialiserSNPs}
\alias{xSocialiserSNPs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSocialiserSNPs(data, ontology = c("EF", "EF_disease", "EF_phenotype",
"EF_bp"), include.LD = NA, LD.r2 = 0.8, measure = c("BM.average",
"BM.max", "BM.complete", "average", "max"), method.term = c("Resnik",
"Lin", "Schlicker", "Jiang", "Pesquita"), rescale = TRUE, force = TRUE,
fast = TRUE, parallel = TRUE, multicores = NULL, path.mode =
c("all_paths", "shortest_paths", "all_shortest_paths"), true.path.rule
= T, verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
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
\item{measure}{
%% ~~Describe \code{measure} here~~
}
\item{method.term}{
%% ~~Describe \code{method.term} here~~
}
\item{rescale}{
%% ~~Describe \code{rescale} here~~
}
\item{force}{
%% ~~Describe \code{force} here~~
}
\item{fast}{
%% ~~Describe \code{fast} here~~
}
\item{parallel}{
%% ~~Describe \code{parallel} here~~
}
\item{multicores}{
%% ~~Describe \code{multicores} here~~
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
function (data, ontology = c("EF", "EF_disease", "EF_phenotype",
"EF_bp"), include.LD = NA, LD.r2 = 0.8, measure = c("BM.average",
"BM.max", "BM.complete", "average", "max"), method.term = c("Resnik",
"Lin", "Schlicker", "Jiang", "Pesquita"), rescale = TRUE,
force = TRUE, fast = TRUE, parallel = TRUE, multicores = NULL,
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
true.path.rule = T, verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
startT <- Sys.time()
message(paste(c("Start at ", as.character(startT)), collapse = ""),
appendLF = T)
message("", appendLF = T)
ontology <- match.arg(ontology)
measure <- match.arg(measure)
method.term <- match.arg(method.term)
path.mode <- match.arg(path.mode)
if (is.vector(data)) {
data <- unique(data)
}
else {
stop("The input data must be a vector.\n")
}
if (!is.na(ontology)) {
if (verbose) {
now <- Sys.time()
message(sprintf("Load the ontology \%s and its SNP annotations (\%s)
...",
ontology, as.character(now)), appendLF = T)
}
ig <- xRDataLoader(RData = paste("ig.EF", sep = ""),
RData.location = RData.location, guid = guid, verbose = verbose)
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
}
else {
g <- ig
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
message(sprintf("'xSocialiser' is being called (\%s):",
as.character(now)), appendLF = T)
message(sprintf("#######################################################",
appendLF = T))
}
res <- xSocialiser(data = data, annotation = anno, g = g,
measure = measure, method.term = method.term, rescale = rescale,
force = force, fast = fast, parallel = parallel, multicores =
multicores,
path.mode = path.mode, true.path.rule = true.path.rule,
verbose = verbose)
if (!is.null(res)) {
dag <- res$dag
sim_ig <- res
tEdges <- igraph::get.data.frame(sim_ig, what = "edges")
sEdges <- tEdges[with(tEdges, order(-weight, from, to)),
]
relations <- apply(sEdges, 1, function(x) {
if (x[1] > x[2]) {
return(cbind(x[2], x[1], x[3]))
}
else {
x
}
})
relations <- t(relations)
colnames(relations) <- c("from", "to", "weight")
res <- igraph::graph.data.frame(d = relations, directed = F)
if (class(res) == "igraph") {
if (!is.null(E(res)$weight)) {
E(res)$weight <- as.numeric(E(res)$weight)
}
}
res$dag <- dag
}
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("'xSocialiser' has been finished (\%s)!",
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
invisible(res)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
