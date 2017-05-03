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
ontology = c("GOBP", "GOMF", "GOCC", "PS", "PS2", "SF", "Pfam", "DO",
"HPPA", "HPMI", "HPCM", "HPMA", "MP", "EF", "MsigdbH", "MsigdbC1",
"MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP", "MsigdbC2KEGG",
"MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR",
"MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC",
"MsigdbC6", "MsigdbC7", "DGIdb", "GTExV4", "GTExV6", "CreedsDisease",
"CreedsDiseaseUP", "CreedsDiseaseDN", "CreedsDrug", "CreedsDrugUP",
"CreedsDrugDN", "CreedsGene", "CreedsGeneUP",
"CreedsGeneDN"), size.range = c(10, 2000), min.overlap = 3,
which.distance = NULL, test = c("hypergeo", "fisher", "binomial"),
p.adjust.method = c("BH", "BY", "bonferroni", "holm", "hochberg",
"hommel"), ontology.algorithm = c("none", "pc", "elim", "lea"),
elim.pvalue = 0.01, lea.depth = 2, path.mode = c("all_paths",
"shortest_paths", "all_shortest_paths"), true.path.rule = F, verbose =
T, RData.location = "http://galahad.well.ox.ac.uk/bigdata")
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
ontology = c("GOBP", "GOMF", "GOCC", "PS", "PS2", "SF", "Pfam",
"DO", "HPPA", "HPMI", "HPCM", "HPMA", "MP", "EF", "MsigdbH",
"MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP",
"MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA",
"MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM",
"MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6",
"MsigdbC7", "DGIdb", "GTExV4", "GTExV6", "CreedsDisease",
"CreedsDiseaseUP", "CreedsDiseaseDN", "CreedsDrug", "CreedsDrugUP",
"CreedsDrugDN", "CreedsGene", "CreedsGeneUP", "CreedsGeneDN"),
size.range = c(10, 2000), min.overlap = 3, which.distance = NULL,
test = c("hypergeo", "fisher", "binomial"), p.adjust.method = c("BH",
"BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm =
c("none",
"pc", "elim", "lea"), elim.pvalue = 0.01, lea.depth = 2,
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
true.path.rule = F, verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata")
{
startT <- Sys.time()
message(paste(c("Start at ", as.character(startT)), collapse = ""),
appendLF = T)
message("", appendLF = T)
ontology <- match.arg(ontology)
test <- match.arg(test)
p.adjust.method <- match.arg(p.adjust.method)
ontology.algorithm <- match.arg(ontology.algorithm)
path.mode <- match.arg(path.mode)
if (length(data) == 0) {
return(FALSE)
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
message(sprintf("Load the ontology \%s and its gene annotations (\%s)
...",
ontology, as.character(now)), appendLF = T)
}
flag_PS2 <- FALSE
if (ontology == "PS2") {
flag_PS2 <- TRUE
ontology <- "PS"
}
GS <- xRDataLoader(RData.customised = paste("org.Hs.eg",
ontology, sep = ""), RData.location = RData.location,
verbose = verbose)
if (flag_PS2) {
tmp <- as.character(unique(GS$set_info$name))
inds <- sapply(tmp, function(x) which(GS$set_info$name ==
x))
set_info <- data.frame()
for (i in 1:length(inds)) {
set_info <- rbind(set_info, as.matrix(GS$set_info[max(inds[[i]]),
]))
}
gs <- list()
for (i in 1:length(inds)) {
gs[[i]] <- unlist(GS$gs[inds[[i]]], use.names = F)
}
names(gs) <- rownames(set_info)
GS$set_info <- set_info
GS$gs <- gs
}
anno <- GS$gs
all.ontologies <- c("GOBP", "GOMF", "GOCC", "DO", "HPPA",
"HPMI", "HPCM", "HPMA", "MP", "EF")
flag_ontology <- ontology \%in\% all.ontologies
if (flag_ontology) {
g <- xRDataLoader(RData.customised = paste("ig.",
ontology, sep = ""), RData.location = RData.location,
verbose = verbose)
}
else {
ontology.algorithm <- "none"
nodes <- data.frame(name = as.character(GS$set_info$setID),
term_id = as.character(GS$set_info$setID), term_name =
as.character(GS$set_info$name),
term_distance = as.character(GS$set_info$distance),
stringsAsFactors = F)
nodes <- rbind(nodes, c("root", "root", "root", "root"))
relations <- data.frame(from = "root", to = nodes$name)
g <- igraph::graph.data.frame(d = relations, directed = T,
vertices = nodes)
}
}
else {
stop("There is no input for the ontology.\n")
}
if (verbose) {
now <- Sys.time()
message(sprintf("Do gene mapping from Symbols to EntrezIDs (\%s) ...",
as.character(now)), appendLF = T)
}
data <- xSymbol2GeneID(data, check.symbol.identity =
check.symbol.identity,
verbose = verbose, RData.location = RData.location)
data <- data[!is.na(data)]
if (length(background) > 0) {
background <- xSymbol2GeneID(background, check.symbol.identity =
check.symbol.identity,
verbose = verbose, RData.location = RData.location)
background <- background[!is.na(background)]
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
which.distance = which.distance, test = test, p.adjust.method =
p.adjust.method,
ontology.algorithm = ontology.algorithm, elim.pvalue = elim.pvalue,
lea.depth = lea.depth, path.mode = path.mode, true.path.rule =
true.path.rule,
verbose = verbose)
if (1 & class(eTerm) == "eTerm") {
EG <- xRDataLoader(RData.customised = paste("org.Hs.eg",
sep = ""), RData.location = RData.location, verbose = verbose)
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
