\name{xConverter}
\alias{xConverter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xConverter(obj, from = c("igraph", "dgCMatrix", "dtree", "lol",
"json"), to = c("dgCMatrix", "igraph", "dtree", "lol", "json",
"igraph_tree"), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{obj}{
%% ~~Describe \code{obj} here~~
}
\item{from}{
%% ~~Describe \code{from} here~~
}
\item{to}{
%% ~~Describe \code{to} here~~
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
function (obj, from = c("igraph", "dgCMatrix", "dtree", "lol",
"json"), to = c("dgCMatrix", "igraph", "dtree", "lol", "json",
"igraph_tree"), verbose = TRUE)
{
from <- match.arg(from)
to <- match.arg(to)
if (from != "igraph" & to != "igraph") {
}
if (from == to) {
warnings(sprintf("Since the class '\%s' converted from is the same as
the class '\%s' converted to, it will return exactly what you
input.\n",
from, to))
return(obj)
}
if (from == "igraph") {
if (to == "dgCMatrix") {
if ("weight" \%in\% list.edge.attributes(obj)) {
E(obj)$weight <- as.numeric(E(obj)$weight)
objConverted <- igraph::as_adjacency_matrix(obj,
type = "both", attr = "weight", edges = F,
names = T, sparse = getIgraphOpt("sparsematrices"))
}
else {
objConverted <- igraph::as_adjacency_matrix(obj,
type = "both", attr = NULL, edges = F, names = T,
sparse = getIgraphOpt("sparsematrices"))
}
}
else if (to \%in\% c("dtree", "lol", "json", "igraph_tree")) {
df_edges <- igraph::get.data.frame(obj, what = "edges")
df <- df_edges
root <- setdiff(df[, 1], df[, 2])
if (length(root) != 1) {
df_root <- data.frame("ROOT", unique(root), stringsAsFactors = F)
if (ncol(df_edges) >= 3) {
for (i in 3:ncol(df_edges)) {
df_root[, i] <- NA
}
}
colnames(df_root) <- colnames(df_edges)
df <- rbind(df_root, df)
}
if (to == "igraph_tree") {
objConverted <- igraph::graph_from_data_frame(d = df,
directed = T)
}
else {
dtree <- data.tree::FromDataFrameNetwork(df)
if (to \%in\% c("lol", "json")) {
lol <- data.tree::ToListExplicit(dtree, unname = TRUE)
if (0) {
func.igraph2list <- function(g, thisNode) {
nm <- igraph::vertex_attr(g, "name", thisNode)
childNodes <- V(g)[which(igraph::shortest.paths(g,
thisNode, mode = "out") == 1)]
if (length(childNodes) == 0) {
return(list(name = nm))
}
list(name = nm, children = unname(lapply(childNodes,
func.igraph2list, g = g)))
}
ig <- igraph::graph_from_data_frame(d = df,
directed = T)
root <- setdiff(df[, 1], df[, 2])
lol <- func.igraph2list(ig, V(ig)[root])
}
if (to == "json") {
json <- jsonlite::toJSON(lol)
objConverted <- json
}
else {
objConverted <- lol
}
}
else {
objConverted <- dtree
}
}
}
}
else if (from == "dgCMatrix") {
if (to == "igraph") {
nodes <- data.frame(name = rownames(obj))
nodenames <- rownames(obj)
adjM <- obj
tmp <- which(as.matrix(adjM != 0), arr.ind = T)
if (from == "dgCMatrix") {
ind <- which(tmp[, 1] < tmp[, 2])
ttmp <- matrix(0, nrow = length(ind), ncol = 2)
ttmp[1:length(ind), ] <- tmp[ind, ]
tmp <- ttmp
}
weight_flag <- T
if (all(adjM[tmp] == 1)) {
weight_flag <- F
}
if (weight_flag) {
relations <- data.frame(from = nodenames[tmp[,
1]], to = nodenames[tmp[, 2]], weight = adjM[tmp])
}
else {
relations <- data.frame(from = nodenames[tmp[,
1]], to = nodenames[tmp[, 2]])
}
if (from == "dgCMatrix") {
objConverted <- igraph::graph_from_data_frame(d = relations,
directed = F, vertices = nodes)
}
}
else {
warnings(sprintf("Conversion between '\%s' and '\%s' is not supported;
instead converted to the 'igraph' object first.\n",
from, to))
return(obj)
}
}
else if (from == "dtree") {
if (to \%in\% c("lol", "json")) {
lol <- data.tree::ToListExplicit(obj, unname = TRUE)
if (to == "json") {
json <- jsonlite::toJSON(lol)
objConverted <- json
}
else {
objConverted <- lol
}
}
else if (to == "igraph") {
df_edges <- data.tree::ToDataFrameNetwork(obj)[,
c(1, 2)]
df_edges[, 1] <- gsub(".*/", "", df_edges[, 1])
df_edges[, 2] <- gsub(".*/", "", df_edges[, 2])
ig <- igraph::graph_from_data_frame(d = df_edges,
directed = T)
objConverted <- ig
}
else {
warnings(sprintf("Conversion between '\%s' and '\%s' is not supported;
instead converted to the 'igraph' object first.\n",
from, to))
return(obj)
}
}
else if (from == "lol") {
if (to == "json") {
json <- jsonlite::toJSON(obj)
objConverted <- json
}
else if (to == "dtree") {
dtree <- data.tree::FromListExplicit(obj)
objConverted <- dtree
}
else {
warnings(sprintf("Conversion between '\%s' and '\%s' is not supported;
instead converted to the 'data.tree' object first.\n",
from, to))
return(obj)
}
}
else if (from == "json") {
if (to \%in\% c("lol", "dtree")) {
lol <- jsonlite::fromJSON(obj, simplifyDataFrame = F)
if (to == "dtree") {
dtree <- data.tree::FromListExplicit(lol)
objConverted <- dtree
}
else {
objConverted <- lol
}
}
else {
warnings(sprintf("Conversion between '\%s' and '\%s' is not supported;
instead converted to the 'data.tree' object first.\n",
from, to))
return(obj)
}
}
if (verbose) {
message(sprintf("Your input object of class '\%s' has been converted
into an object of class '\%s'.",
from, to), appendLF = T)
}
return(objConverted)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
