\name{xEnricher}
\alias{xEnricher}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnricher(data, annotation, g, background = NULL, size.range = c(10,
2000), min.overlap = 3, which.distance = NULL, test = c("hypergeo",
"fisher", "binomial"), p.adjust.method = c("BH", "BY", "bonferroni",
"holm", "hochberg", "hommel"), ontology.algorithm = c("none", "pc",
"elim", "lea"), elim.pvalue = 0.01, lea.depth = 2, path.mode =
c("all_paths", "shortest_paths", "all_shortest_paths"), true.path.rule
= TRUE, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{annotation}{
%% ~~Describe \code{annotation} here~~
}
\item{g}{
%% ~~Describe \code{g} here~~
}
\item{background}{
%% ~~Describe \code{background} here~~
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
function (data, annotation, g, background = NULL, size.range = c(10,
2000), min.overlap = 3, which.distance = NULL, test = c("hypergeo",
"fisher", "binomial"), p.adjust.method = c("BH", "BY", "bonferroni",
"holm", "hochberg", "hommel"), ontology.algorithm = c("none",
"pc", "elim", "lea"), elim.pvalue = 0.01, lea.depth = 2,
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
true.path.rule = TRUE, verbose = T)
{
test <- match.arg(test)
p.adjust.method <- match.arg(p.adjust.method)
ontology.algorithm <- match.arg(ontology.algorithm)
path.mode <- match.arg(path.mode)
if (is.vector(data)) {
data <- unique(data)
}
else {
stop("The input data must be a vector.\n")
}
if (class(annotation) == "GS") {
originAnnos <- annotation$gs
}
else if (class(annotation) == "list") {
originAnnos <- annotation
}
else if (class(annotation) == "dgCMatrix") {
D <- annotation
originAnnos <- sapply(1:ncol(D), function(j) {
names(which(D[, j] != 0))
})
names(originAnnos) <- colnames(annotation)
}
else {
stop("The input annotation must be either 'GS' or 'list' or 'dgCMatrix'
object.\n")
}
annotation <- originAnnos
ig <- g
if (class(ig) != "igraph") {
stop("The function must apply to the 'igraph' object.\n")
}
else {
if (verbose) {
now <- Sys.time()
message(sprintf("First, generate a subgraph induced (via '\%s' mode) by
the annotation data (\%s) ...",
path.mode, as.character(now)), appendLF = T)
}
subg <- xDAGanno(g = ig, annotation = annotation, path.mode =
path.mode,
true.path.rule = true.path.rule, verbose = verbose)
gs <- V(subg)$anno
names(gs) <- V(subg)$name
gs.distance <- V(subg)$term_distance
names(gs.distance) <- V(subg)$name
}
if (1) {
if (is.vector(background)) {
background <- base::unique(background)
background <- background[!is.null(background)]
background <- background[!is.na(background)]
}
if (length(background) > 0) {
if (1) {
background <- base::union(background, data)
}
gs <- lapply(gs, function(x) {
ind <- match(x, background)
x[!is.na(ind)]
})
}
}
if (!is.null(which.distance) & sum(is.na(gs.distance)) ==
0) {
distance_filtered <- lapply(which.distance, function(x) {
names(gs)[(gs.distance == as.integer(x))]
})
distance_filtered <- unlist(distance_filtered)
}
else {
distance_filtered <- names(gs)
}
ind.distance <- match(distance_filtered, names(gs))
gs.length <- sapply(gs, length)
ind.length <- which(gs.length >= size.range[1] & gs.length <=
size.range[2])
ind <- intersect(ind.distance, ind.length)
gs <- gs[ind]
if (length(gs) == 0) {
warnings("There are no terms being used.\n")
return(F)
}
doFisherTest <- function(genes.group, genes.term, genes.universe) {
genes.hit <- intersect(genes.group, genes.term)
X <- length(genes.hit)
K <- length(genes.group)
M <- length(genes.term)
N <- length(genes.universe)
cTab <- matrix(c(X, K - X, M - X, N - M - K + X), nrow = 2,
dimnames = list(c("anno", "notAnno"), c("group",
"notGroup")))
p.value <- ifelse(all(cTab == 0), 1, stats::fisher.test(cTab,
alternative = "greater")$p.value)
return(p.value)
}
doHypergeoTest <- function(genes.group, genes.term, genes.universe) {
genes.hit <- intersect(genes.group, genes.term)
X <- length(genes.hit)
K <- length(genes.group)
M <- length(genes.term)
N <- length(genes.universe)
x <- X
m <- M
n <- N - M
k <- K
p.value <- ifelse(m == 0 || k == 0, 1, stats::phyper(x,
m, n, k, lower.tail = F, log.p = F))
return(p.value)
}
doBinomialTest <- function(genes.group, genes.term, genes.universe) {
genes.hit <- intersect(genes.group, genes.term)
X <- length(genes.hit)
K <- length(genes.group)
M <- length(genes.term)
N <- length(genes.universe)
p.value <- ifelse(K == 0 || M == 0 || N == 0, 1, stats::pbinom(X,
K, M/N, lower.tail = F, log.p = F))
return(p.value)
}
zscoreHyper <- function(genes.group, genes.term, genes.universe) {
genes.hit <- intersect(genes.group, genes.term)
X <- length(genes.hit)
K <- length(genes.group)
M <- length(genes.term)
N <- length(genes.universe)
if (1) {
x.exp <- K * M/N
var.exp <- K * M/N * (N - M)/N * (N - K)/(N - 1)
if (is.na(var.exp)) {
z <- NA
}
else {
if (var.exp != 0) {
suppressWarnings(z <- (X - x.exp)/sqrt(var.exp))
}
else {
z <- NA
}
}
}
else {
x <- X
m <- M
n <- N - M
k <- K
suppressWarnings(d <- stats::dhyper(x, m, n, k, log = TRUE) -
log(2))
suppressWarnings(pupper <- stats::phyper(x, m, n,
k, lower.tail = FALSE, log.p = TRUE))
suppressWarnings(plower <- stats::phyper(x - 1, m,
n, k, lower.tail = TRUE, log.p = TRUE))
d[is.na(d)] <- -Inf
pupper[is.na(pupper)] <- -Inf
plower[is.na(plower)] <- -Inf
a <- pupper
b <- d - pupper
a[b > 0] <- d[b > 0]
b <- -abs(b)
pmidupper <- a + log1p(exp(b))
pmidupper[is.infinite(a)] <- a[is.infinite(a)]
a <- plower
b <- d - plower
a[b > 0] <- d[b > 0]
b <- -abs(b)
pmidlower <- a + log1p(exp(b))
pmidlower[is.infinite(a)] <- a[is.infinite(a)]
up <- pmidupper < pmidlower
if (any(up))
z <- stats::qnorm(pmidupper, lower.tail = FALSE,
log.p = TRUE)
if (any(!up))
z <- stats::qnorm(pmidlower, lower.tail = TRUE,
log.p = TRUE)
}
return(z)
}
fcHyper <- function(genes.group, genes.term, genes.universe) {
genes.hit <- intersect(genes.group, genes.term)
X <- length(genes.hit)
K <- length(genes.group)
M <- length(genes.term)
N <- length(genes.universe)
x.exp <- K * M/N
fc <- X/x.exp
return(fc)
}
if (verbose) {
now <- Sys.time()
message(sprintf("Next, prepare enrichment analysis (\%s) ...",
as.character(now)), appendLF = T)
}
terms <- names(gs)
genes.universe <- unique(unlist(gs[terms]))
genes.group <- intersect(genes.universe, data)
if (length(genes.group) == 0) {
warnings("There is no gene being used.\n")
return(F)
}
else {
if (verbose) {
now <- Sys.time()
message(sprintf("\tThere are \%d genes/SNPs of interest tested against
\%d genes/SNPs as the background (\%s)",
length(genes.group), length(genes.universe),
as.character(now)), appendLF = T)
}
}
subg <- dnet::dDAGinduce(g = subg, nodes_query = terms, path.mode =
path.mode)
set_info <- data.frame(id = V(subg)$term_id, name = V(subg)$term_name,
distance = V(subg)$term_distance, row.names = V(subg)$name)
if (ontology.algorithm == "none") {
if (verbose) {
now <- Sys.time()
message(sprintf("Third, perform enrichment analysis using '\%s' test
(\%s) ...",
test, as.character(now)), appendLF = T)
if (is.null(which.distance)) {
message(sprintf("\tThere are \%d terms being used, each restricted
within [\%s] annotations",
length(terms), paste(size.range, collapse = ",")),
appendLF = T)
}
else {
message(sprintf("\tThere are \%d terms being used, each restricted
within [\%s] annotations and [\%s] distance",
length(terms), paste(size.range, collapse = ","),
paste(which.distance, collapse = ",")), appendLF = T)
}
}
pvals <- sapply(terms, function(term) {
genes.term <- unique(unlist(gs[term]))
p.value <- switch(test, fisher = doFisherTest(genes.group,
genes.term, genes.universe), hypergeo = doHypergeoTest(genes.group,
genes.term, genes.universe), binomial = doBinomialTest(genes.group,
genes.term, genes.universe))
})
zscores <- sapply(terms, function(term) {
genes.term <- unique(unlist(gs[term]))
zscoreHyper(genes.group, genes.term, genes.universe)
})
fcs <- sapply(terms, function(term) {
genes.term <- unique(unlist(gs[term]))
fcHyper(genes.group, genes.term, genes.universe)
})
}
else if (ontology.algorithm == "pc" || ontology.algorithm ==
"elim" || ontology.algorithm == "lea") {
if (verbose) {
now <- Sys.time()
message(sprintf("Third, perform enrichment analysis based on '\%s'
test, and also using '\%s' algorithm to respect ontology structure
(\%s) ...",
test, ontology.algorithm, as.character(now)),
appendLF = T)
}
if (verbose) {
message(sprintf("\tThere are \%d terms being used",
length(V(subg))), appendLF = T)
}
level2node <- dnet::dDAGlevel(subg, level.mode = "longest_path",
return.mode = "level2node")
level2node.Hash <- list2env(level2node)
nLevels <- length(level2node)
node2pval.Hash <- new.env(hash = T, parent = emptyenv())
node2zscore.Hash <- new.env(hash = T, parent = emptyenv())
node2fc.Hash <- new.env(hash = T, parent = emptyenv())
if (ontology.algorithm == "pc") {
for (i in nLevels:2) {
currNodes <- get(as.character(i), envir = level2node.Hash,
mode = "character")
for (currNode in currNodes) {
genes.term <- unique(unlist(gs[currNode]))
pvalue_whole <- switch(test, fisher = doFisherTest(genes.group,
genes.term, genes.universe), hypergeo = doHypergeoTest(genes.group,
genes.term, genes.universe), binomial = doBinomialTest(genes.group,
genes.term, genes.universe))
zscore_whole <- zscoreHyper(genes.group, genes.term,
genes.universe)
fc_whole <- fcHyper(genes.group, genes.term,
genes.universe)
neighs.in <- igraph::neighborhood(subg, order = 1,
nodes = currNode, mode = "in")
adjNodes <- setdiff(V(subg)[unlist(neighs.in)]$name,
currNode)
genes.parent <- unique(unlist(gs[adjNodes]))
genes.group.parent <- intersect(genes.group,
genes.parent)
genes.term.parent <- intersect(genes.term,
genes.parent)
pvalue_relative <- switch(test, fisher =
doFisherTest(genes.group.parent,
genes.term.parent, genes.parent), hypergeo =
doHypergeoTest(genes.group.parent,
genes.term.parent, genes.parent), binomial =
doBinomialTest(genes.group.parent,
genes.term.parent, genes.parent))
zscore_relative <- zscoreHyper(genes.group.parent,
genes.term.parent, genes.parent)
fc_relative <- fcHyper(genes.group.parent,
genes.term.parent, genes.parent)
pvalue <- max(pvalue_whole, pvalue_relative)
assign(currNode, pvalue, envir = node2pval.Hash)
zscore <- ifelse(pvalue_whole > pvalue_relative,
zscore_whole, zscore_relative)
assign(currNode, zscore, envir = node2zscore.Hash)
fc <- ifelse(pvalue_whole > pvalue_relative,
fc_whole, fc_relative)
assign(currNode, fc, envir = node2fc.Hash)
}
if (verbose) {
message(sprintf("\tAt level \%d, there are \%d nodes/terms",
i, length(currNodes), appendLF = T))
}
}
root <- dnet::dDAGroot(subg)
assign(root, 1, envir = node2pval.Hash)
assign(root, 0, envir = node2zscore.Hash)
}
else if (ontology.algorithm == "elim") {
sigNode2pval.Hash <- new.env(hash = T, parent = emptyenv())
ancNode2gene.Hash <- new.env(hash = T, parent = emptyenv())
if (is.null(elim.pvalue) || is.na(elim.pvalue) ||
elim.pvalue > 1 || elim.pvalue < 0) {
elim.pvalue <- 0.01
}
pval.cutoff <- elim.pvalue
for (i in nLevels:1) {
currNodes <- get(as.character(i), envir = level2node.Hash,
mode = "character")
currAnno <- gs[currNodes]
for (currNode in currNodes) {
genes.term <- unique(unlist(gs[currNode]))
if (exists(currNode, envir = ancNode2gene.Hash,
mode = "numeric")) {
genes.elim <- get(currNode, envir = ancNode2gene.Hash,
mode = "numeric")
genes.term <- setdiff(genes.term, genes.elim)
}
pvalue <- switch(test, fisher = doFisherTest(genes.group,
genes.term, genes.universe), hypergeo = doHypergeoTest(genes.group,
genes.term, genes.universe), binomial = doBinomialTest(genes.group,
genes.term, genes.universe))
zscore <- zscoreHyper(genes.group, genes.term,
genes.universe)
fc <- fcHyper(genes.group, genes.term, genes.universe)
assign(currNode, pvalue, envir = node2pval.Hash)
assign(currNode, zscore, envir = node2zscore.Hash)
assign(currNode, fc, envir = node2fc.Hash)
if (pvalue < pval.cutoff) {
assign(currNode, pvalue, envir = sigNode2pval.Hash)
elimGenesID <- currAnno[[currNode]]
dag.ancestors <- dnet::dDAGinduce(subg, currNode,
path.mode = "all_paths")
ancestors <- setdiff(V(dag.ancestors)$name,
currNode)
oldAncestors2GenesID <- sapply(ancestors,
function(ancestor) {
if (exists(ancestor, envir = ancNode2gene.Hash,
mode = "numeric")) {
get(ancestor, envir = ancNode2gene.Hash,
mode = "numeric")
}
})
newAncestors2GenesID <- lapply(oldAncestors2GenesID,
function(oldGenes) {
base::union(oldGenes, elimGenesID)
})
if (length(newAncestors2GenesID) > 0) {
sapply(names(newAncestors2GenesID), function(ancestor) {
assign(ancestor, newAncestors2GenesID[[ancestor]],
envir = ancNode2gene.Hash)
})
}
}
}
if (verbose) {
num.signodes <- length(ls(sigNode2pval.Hash))
num.ancnodes <- length(ls(ancNode2gene.Hash))
num.elimgenes <- length(unique(unlist(as.list(ancNode2gene.Hash))))
message(sprintf("\tAt level \%d, there are \%d nodes/terms: up to \%d
significant nodes, \%d ancestral nodes changed (\%d genes/SNPs
eliminated)",
i, length(currNodes), num.signodes, num.ancnodes,
num.elimgenes), appendLF = T)
}
}
}
else if (ontology.algorithm == "lea") {
node2pvalo.Hash <- new.env(hash = T, parent = emptyenv())
if (is.null(lea.depth) || is.na(lea.depth) || lea.depth <
0) {
lea.depth <- 2
}
depth.cutoff <- as.integer(lea.depth)
for (i in nLevels:1) {
currNodes <- get(as.character(i), envir = level2node.Hash,
mode = "character")
currAnno <- gs[currNodes]
num.recalculate <- 0
for (currNode in currNodes) {
genes.term <- unique(unlist(gs[currNode]))
pvalue.old <- switch(test, fisher = doFisherTest(genes.group,
genes.term, genes.universe), hypergeo = doHypergeoTest(genes.group,
genes.term, genes.universe), binomial = doBinomialTest(genes.group,
genes.term, genes.universe))
zscore.old <- zscoreHyper(genes.group, genes.term,
genes.universe)
fc.old <- fcHyper(genes.group, genes.term,
genes.universe)
assign(currNode, pvalue.old, envir = node2pvalo.Hash)
neighs.out <- igraph::neighborhood(subg, order = depth.cutoff,
nodes = currNode, mode = "out")
adjNodes <- setdiff(V(subg)[unlist(neighs.out)]$name,
currNode)
if (length(adjNodes) != 0) {
if (1) {
pvalue.children <- sapply(adjNodes, function(child) {
if (exists(child, envir = node2pvalo.Hash,
mode = "numeric")) {
get(child, envir = node2pvalo.Hash,
mode = "numeric")
}
})
}
else {
pvalue.children <- sapply(adjNodes, function(child) {
if (exists(child, envir = node2pval.Hash,
mode = "numeric")) {
get(child, envir = node2pval.Hash,
mode = "numeric")
}
})
}
chNodes <- names(pvalue.children[pvalue.children <
pvalue.old])
if (length(chNodes) > 0) {
num.recalculate <- num.recalculate + 1
genes.elim <- unique(unlist(gs[chNodes]))
genes.term.new <- setdiff(genes.term, genes.elim)
pvalue.new <- switch(test, fisher = doFisherTest(genes.group,
genes.term.new, genes.universe), hypergeo = doHypergeoTest(genes.group,
genes.term.new, genes.universe), binomial = doBinomialTest(genes.group,
genes.term.new, genes.universe))
zscore.new <- zscoreHyper(genes.group,
genes.term.new, genes.universe)
fc.new <- fcHyper(genes.group, genes.term.new,
genes.universe)
pvalue <- max(pvalue.new, pvalue.old)
zscore <- ifelse(pvalue.new > pvalue.old,
zscore.new, zscore.old)
fc <- ifelse(pvalue.new > pvalue.old, fc.new,
fc.old)
}
else {
pvalue <- pvalue.old
zscore <- zscore.old
fc <- fc.old
}
}
else {
pvalue <- pvalue.old
zscore <- zscore.old
fc <- fc.old
}
assign(currNode, pvalue, envir = node2pval.Hash)
assign(currNode, zscore, envir = node2zscore.Hash)
assign(currNode, fc, envir = node2fc.Hash)
}
if (verbose) {
message(sprintf("\tAt level \%d, there are \%d nodes/terms and \%d
being recalculated",
i, length(currNodes), num.recalculate), appendLF = T)
}
}
}
pvals <- unlist(as.list(node2pval.Hash))
zscores <- unlist(as.list(node2zscore.Hash))
fcs <- unlist(as.list(node2fc.Hash))
}
overlaps <- lapply(names(gs), function(term) {
genes.term <- unique(unlist(gs[term]))
x <- intersect(genes.group, genes.term)
x
})
names(overlaps) <- names(gs)
flag_filter <- sapply(overlaps, function(x) ifelse(length(x) >=
min.overlap, T, F))
if (sum(flag_filter) == 0) {
warnings("It seems there are no terms meeting the specified
'size.range' and 'min.overlap'.\n")
return(F)
}
gs <- gs[flag_filter]
overlaps <- overlaps[flag_filter]
common <- intersect(names(gs), names(zscores))
ind_gs <- match(common, names(gs))
ind_zscores <- match(common, names(zscores))
gs <- gs[ind_gs[!is.na(ind_gs)]]
overlaps <- overlaps[ind_gs[!is.na(ind_gs)]]
zscores <- zscores[ind_zscores[!is.na(ind_zscores)]]
fcs <- fcs[ind_zscores[!is.na(ind_zscores)]]
pvals <- pvals[ind_zscores[!is.na(ind_zscores)]]
flag <- !is.na(zscores)
gs <- gs[flag]
overlaps <- overlaps[flag]
zscores <- zscores[flag]
fcs <- fcs[flag]
pvals <- pvals[flag]
if (length(pvals) == 0) {
warnings("There are no pvals being calcualted.\n")
return(F)
}
ind <- match(rownames(set_info), names(pvals))
set_info <- set_info[!is.na(ind), ]
zscores <- signif(zscores, digits = 3)
fcs <- signif(fcs, digits = 3)
pvals <- sapply(pvals, function(x) min(x, 1))
if (verbose) {
now <- Sys.time()
message(sprintf("Last, adjust the p-values for \%d terms (with \%d
minimum overlaps) using the \%s method (\%s) ...",
length(pvals), min.overlap, p.adjust.method, as.character(now)),
appendLF = T)
}
adjpvals <- stats::p.adjust(pvals, method = p.adjust.method)
pvals <- signif(pvals, digits = 2)
adjpvals <- sapply(adjpvals, function(x) min(x, 1))
adjpvals <- signif(adjpvals, digits = 2)
if (0) {
tmp <- as.numeric(format(.Machine)["double.xmin"])
tmp <- signif(tmp, digits = 2)
pvals[pvals < tmp] <- tmp
adjpvals[adjpvals < tmp] <- tmp
}
pvals <- sapply(pvals, function(x) {
if (x < 0.1 & x != 0) {
as.numeric(format(x, scientific = T))
}
else {
x
}
})
adjpvals <- sapply(adjpvals, function(x) {
if (x < 0.1 & x != 0) {
as.numeric(format(x, scientific = T))
}
else {
x
}
})
cross <- matrix(0, nrow = length(overlaps), ncol = length(overlaps))
if (length(overlaps) >= 2) {
for (i in seq(1, length(overlaps) - 1)) {
x1 <- overlaps[[i]]
for (j in seq(i + 1, length(overlaps))) {
x2 <- overlaps[[j]]
cross[i, j] <- length(intersect(x1, x2))
cross[j, i] <- length(intersect(x1, x2))
}
}
colnames(cross) <- rownames(cross) <- names(overlaps)
diag(cross) <- sapply(overlaps, length)
}
eTerm <- list(term_info = set_info, annotation = gs, g = subg,
data = genes.group, background = genes.universe, overlap = overlaps,
fc = fcs, zscore = zscores, pvalue = pvals, adjp = adjpvals,
cross = cross, call = match.call())
class(eTerm) <- "eTerm"
invisible(eTerm)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
