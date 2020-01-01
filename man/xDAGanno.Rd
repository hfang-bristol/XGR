\name{xDAGanno}
\alias{xDAGanno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xDAGanno(g, annotation, path.mode = c("all_paths", "shortest_paths",
"all_shortest_paths"), true.path.rule = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{g}{
%% ~~Describe \code{g} here~~
}
\item{annotation}{
%% ~~Describe \code{annotation} here~~
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
function (g, annotation, path.mode = c("all_paths", "shortest_paths",
"all_shortest_paths"), true.path.rule = TRUE, verbose = TRUE)
{
path.mode <- match.arg(path.mode)
ig <- g
if (class(ig) != "igraph") {
stop("The function must apply to the 'igraph' object.\n")
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
if (is.list(originAnnos)) {
originNodes <- names(originAnnos)
ind <- match(originNodes, V(ig)$name)
nodes_mapped <- originNodes[!is.na(ind)]
if (length(nodes_mapped) == 0) {
stop("The input annotation data do not contain terms matched to the
nodes/terms in the input graph.\n")
}
}
dag <- dnet::dDAGinduce(ig, originNodes, path.mode = path.mode)
allNodes <- V(dag)$name
node2domain.HoH <- new.env(hash = T, parent = emptyenv())
lapply(allNodes, function(node) {
e <- new.env(hash = T, parent = emptyenv())
if (node \%in\% originNodes) {
sapply(originAnnos[[node]], function(x) {
assign(as.character(x), "o", envir = e)
})
}
assign(node, e, envir = node2domain.HoH)
})
if (true.path.rule) {
level2node <- dnet::dDAGlevel(dag, level.mode = "longest_path",
return.mode = "level2node")
level2node.Hash <- list2env(level2node)
nLevels <- length(level2node)
for (i in nLevels:1) {
currNodes <- get(as.character(i), envir = level2node.Hash,
mode = "character")
adjNodesList <- lapply(currNodes, function(node) {
neighs.in <- igraph::neighborhood(dag, order = 1,
nodes = node, mode = "in")
setdiff(V(dag)[unlist(neighs.in)]$name, node)
})
names(adjNodesList) <- currNodes
lapply(currNodes, function(node) {
domainsID <- ls(get(node, envir = node2domain.HoH,
mode = "environment"))
lapply(adjNodesList[[node]], function(adjNode) {
adjEnv <- get(adjNode, envir = node2domain.HoH,
mode = "environment")
sapply(domainsID, function(domainID) {
assign(domainID, "i", envir = adjEnv)
})
})
})
if (verbose) {
message(sprintf("\tAt level \%d, there are \%d nodes, and \%d incoming
neighbors.",
i, length(currNodes), length(unique(unlist(adjNodesList)))),
appendLF = T)
}
}
}
node2domains <- as.list(node2domain.HoH)[allNodes]
domain_annotation <- lapply(node2domains, function(node) {
vec <- unlist(as.list(node))
res <- names(vec)
names(res) <- vec
res
})
V(dag)$anno <- domain_annotation
counts <- sapply(V(dag)$anno, length)
IC <- -1 * log10(counts/max(counts))
if (1) {
IC[is.infinite(IC)] <- 0
}
V(dag)$IC <- IC
return(dag)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
