\name{xEnricherGenes}
\alias{xEnricherGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnricherGenes(data, background = NULL, check.symbol.identity = F,
ontology = NA, ontology.customised = NULL, size.range = c(10, 2000),
min.overlap = 5, which.distance = NULL, test = c("fisher", "hypergeo",
"binomial"), background.annotatable.only = NULL, p.tail = c("one-tail",
"two-tails"), p.adjust.method = c("BH", "BY", "bonferroni", "holm",
"hochberg", "hommel"), ontology.algorithm = c("none", "pc", "elim",
"lea"), elim.pvalue = 0.01, lea.depth = 2, path.mode = c("all_paths",
"shortest_paths", "all_shortest_paths"), true.path.rule = F, verbose =
T, silent = F, RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{background}{
%% ~~Describe \code{background} here~~
}
\item{check.symbol.identity}{
%% ~~Describe \code{check.symbol.identity} here~~
}
\item{ontology}{
%% ~~Describe \code{ontology} here~~
}
\item{ontology.customised}{
%% ~~Describe \code{ontology.customised} here~~
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
function (data, background = NULL, check.symbol.identity = F,
ontology = NA, ontology.customised = NULL, size.range = c(10,
2000), min.overlap = 5, which.distance = NULL, test = c("fisher",
"hypergeo", "binomial"), background.annotatable.only = NULL,
p.tail = c("one-tail", "two-tails"), p.adjust.method = c("BH",
"BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm =
c("none",
"pc", "elim", "lea"), elim.pvalue = 0.01, lea.depth = 2,
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
true.path.rule = F, verbose = T, silent = F, RData.location =
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
ontology <- ontology[1]
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
warnings("The input data must be a vector.\n")
return(NULL)
}
data <- as.character(data)
aOnto <- xDefineOntology(ontology, ontology.customised =
ontology.customised,
verbose = verbose, RData.location = RData.location, guid = guid)
g <- aOnto$g
anno <- aOnto$anno
if (is.null(g)) {
warnings("There is no input for the ontology.\n")
return(NULL)
}
if (is.null(ontology.customised)) {
if (verbose) {
now <- Sys.time()
message(sprintf("Do gene mapping from Symbols to EntrezIDs (\%s) ...",
as.character(now)), appendLF = T)
}
data <- xSymbol2GeneID(data, check.symbol.identity =
check.symbol.identity,
verbose = verbose, RData.location = RData.location,
guid = guid)
data <- data[!is.na(data)]
if (length(background) > 0) {
background <- xSymbol2GeneID(background, check.symbol.identity =
check.symbol.identity,
verbose = verbose, RData.location = RData.location,
guid = guid)
background <- background[!is.na(background)]
}
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
if (is.null(ontology.customised) & class(eTerm) == "eTerm") {
EG <- xRDataLoader(RData.customised = paste("org.Hs.eg",
sep = ""), RData.location = RData.location, guid = guid,
verbose = verbose)
allGeneID <- EG$gene_info$GeneID
allSymbol <- as.vector(EG$gene_info$Symbol)
overlap <- eTerm$overlap
overlap_symbols <- lapply(overlap, function(x) {
ind <- match(x, allGeneID)
allSymbol[ind]
})
eTerm$overlap <- overlap_symbols
eTerm$data <- allSymbol[match(eTerm$data, allGeneID)]
eTerm$background <- allSymbol[match(eTerm$background,
allGeneID)]
annotation <- eTerm$annotation
annotation_symbols <- lapply(annotation, function(x) {
ind <- match(x, allGeneID)
allSymbol[ind]
})
eTerm$annotation <- annotation_symbols
}
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
message(paste(c("Runtime in total (xEnricherGenes): ",
runTime, " secs\n"), collapse = ""), appendLF = TRUE)
}
invisible(eTerm)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
