\name{xAddCoords}
\alias{xAddCoords}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xAddCoords(g, node.attr = NULL, glayout = layout_with_kk,
edge.color.alternative = c("grey70", "grey95"), seed = 825, verbose =
TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{g}{
%% ~~Describe \code{g} here~~
}
\item{node.attr}{
%% ~~Describe \code{node.attr} here~~
}
\item{glayout}{
%% ~~Describe \code{glayout} here~~
}
\item{edge.color.alternative}{
%% ~~Describe \code{edge.color.alternative} here~~
}
\item{seed}{
%% ~~Describe \code{seed} here~~
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
function (g, node.attr = NULL, glayout = layout_with_kk,
edge.color.alternative = c("grey70",
"grey95"), seed = 825, verbose = TRUE)
{
if (class(g) == "graphNEL") {
ig <- igraph.from.graphNEL(g)
}
else {
ig <- g
}
if (class(ig) != "igraph") {
stop("The function must apply to either 'igraph' or 'graphNEL'
object.\n")
}
node.attrs <- igraph::vertex_attr_names(ig)
flag <- F
if (is.null(node.attr)) {
flag <- T
}
else if (!(node.attr \%in\% node.attrs)) {
flag <- T
}
if (flag) {
message("The igraph object must have a vertex attribute of
interest.\n")
return(NULL)
}
if (verbose) {
now <- Sys.time()
message(sprintf("The igraph object has \%d nodes and \%d edges (\%s)
...",
vcount(ig), ecount(ig), as.character(now)), appendLF = T)
}
type <- community <- contribution <- name <- NULL
df_nodes <- igraph::get.data.frame(ig, what = "vertices")
df_nodes$node.attr <- df_nodes[, node.attr]
if (all(c("type", node.attr, "contribution") \%in\%
colnames(df_nodes))) {
df_nodes <- df_nodes \%>\% dplyr::arrange(type, node.attr,
-contribution)
}
else {
df_nodes <- df_nodes \%>\% dplyr::arrange(node.attr)
}
if (1) {
ls_tmp <- split(x = df_nodes$name, f = df_nodes$node.attr)
ls_ig <- lapply(ls_tmp, function(x) {
g <- dnet::dNetInduce(ig, nodes_query = x, knn = 0,
largest.comp = TRUE)
})
if (1) {
layouts <- lapply(ls_ig, glayout)
}
else {
layouts <- lapply(ls_ig, function(g) {
igraph::layout_as_bipartite(g, types = as.logical(V(g)$type ==
"xnode"))
igraph::layout_in_circle(g)
igraph::layout_with_kk(g)
})
}
set.seed(seed)
glayout <- igraph::merge_coords(ls_ig, layouts)
node.xcoord <- glayout[, 1]
node.ycoord <- glayout[, 2]
if (max(node.xcoord) != min(node.xcoord)) {
node.xcoord <- (node.xcoord - min(node.xcoord))/(max(node.xcoord) -
min(node.xcoord)) * 2 - 1
}
if (max(node.ycoord) != min(node.ycoord)) {
node.ycoord <- (node.ycoord - min(node.ycoord))/(max(node.ycoord) -
min(node.ycoord)) * 2 - 1
}
glayout <- cbind(node.xcoord, node.ycoord)
g_tmp <- igraph::disjoint_union(ls_ig)
ind <- match(V(ig)$name, V(g_tmp)$name)
glayout <- glayout[ind, ]
V(ig)$xcoord <- glayout[, 1]
V(ig)$ycoord <- glayout[, 2]
m <- df_nodes$node.attr
names(m) <- df_nodes$name
df_edges <- igraph::get.data.frame(ig, what = "edges")
ind <- match(df_edges$from, names(m))
m1 <- m[ind]
ind <- match(df_edges$to, names(m))
m2 <- m[ind]
res <- m1 != m2
if (!is.null(names(m1))) {
names(res) <- paste(names(m1), names(m2), sep = "|")
}
if (length(edge.color.alternative) == 2) {
edge.color.alternative <- edge.color.alternative
}
else {
edge.color.alternative <- c("grey70", "grey95")
}
edge.color <- edge.color.alternative[res + 1]
E(ig)$color <- edge.color
}
invisible(ig)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
