\name{xDAGpropagate}
\alias{xDAGpropagate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xDAGpropagate(g, annotation, path.mode = c("all_paths",
"shortest_paths", "all_shortest_paths"), propagation = c("all", "min",
"max"), verbose = TRUE)
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
\item{propagation}{
%% ~~Describe \code{propagation} here~~
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
"all_shortest_paths"), propagation = c("all", "min", "max"),
verbose = TRUE)
{
path.mode <- match.arg(path.mode)
propagation <- match.arg(propagation)
ig <- g
if (class(ig) != "igraph") {
stop("The function must apply to the 'igraph' object.\n")
}
if (class(annotation) == "list") {
oAnnos <- annotation
}
else {
if (class(annotation) == "data.frame") {
annotation <- xSparseMatrix(annotation, verbose = FALSE)
}
if (class(annotation) == "dgCMatrix") {
D <- annotation
oAnnos <- lapply(1:ncol(D), function(j) {
ind <- which(D[, j] != 0)
x <- D[ind, j]
names(x) <- names(ind)
x
})
names(oAnnos) <- colnames(annotation)
}
else {
stop("The input annotation must be either 'GS' or 'list' or 'dgCMatrix'
object.\n")
}
}
if (is.list(oAnnos)) {
originNodes <- names(oAnnos)
ind <- match(originNodes, V(ig)$name)
nodes_mapped <- originNodes[!is.na(ind)]
if (length(nodes_mapped) == 0) {
stop("The input annotation data do not contain terms matched to the
nodes/terms in the input graph.\n")
}
}
dag <- dnet::dDAGinduce(ig, originNodes, path.mode = path.mode)
allNodes <- V(dag)$name
if (verbose) {
now <- Sys.time()
message(sprintf("Do propagation via '\%s' operation (\%s) ...",
propagation, as.character(now)), appendLF = T)
}
if (1) {
pAnnos <- oAnnos[allNodes]
names(pAnnos) <- allNodes
level2node <- dnet::dDAGlevel(dag, level.mode = "longest_path",
return.mode = "level2node")
nLevels <- length(level2node)
for (i in nLevels:2) {
currNodes <- level2node[[i]]
adjNodesList <- lapply(currNodes, function(node) {
neighs.in <- igraph::neighborhood(dag, order = 1,
nodes = node, mode = "in")
setdiff(V(dag)[unlist(neighs.in)]$name, node)
})
names(adjNodesList) <- currNodes
for (k in 1:length(currNodes)) {
node <- currNodes[k]
nowDomain <- pAnnos[[node]]
nowDomain_mat <- cbind(names(nowDomain), as.numeric(nowDomain))
adjNodes <- adjNodesList[[node]]
res <- lapply(adjNodes, function(adjNode) {
adjDomain <- pAnnos[[adjNode]]
adjDomain_mat <- cbind(names(adjDomain), as.numeric(adjDomain))
all_mat <- rbind(nowDomain_mat, adjDomain_mat)
all_list <- base::split(x = as.numeric(all_mat[,
2]), f = all_mat[, 1])
output_list <- lapply(1:length(all_list), function(i) {
x <- all_list[[i]]
if (propagation == "max") {
y <- max(x)
}
else if (propagation == "min") {
y <- min(x)
}
else {
y <- unique(x)
}
data.frame(name = rep(names(all_list)[i],
length(y)), value = y, stringsAsFactors = FALSE)
})
if (0) {
x_mat <- base::do.call(base::rbind, output_list)
output <- as.vector(x_mat)
names(output) <- rownames(x_mat)
}
else {
x_mat <- base::do.call(base::rbind, output_list)
output <- x_mat$value
names(output) <- x_mat$name
}
return(output)
})
pAnnos[adjNodes] <- res
}
if (verbose) {
message(sprintf("\tAt level \%d, there are \%d nodes, and \%d incoming
neighbors (\%s).",
i, length(currNodes), length(unique(unlist(adjNodesList))),
as.character(Sys.time())), appendLF = T)
}
}
}
ind <- match(V(dag)$name, names(pAnnos))
V(dag)$anno <- pAnnos[ind]
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
