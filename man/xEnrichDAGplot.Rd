\name{xEnrichDAGplot}
\alias{xEnrichDAGplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichDAGplot(eTerm, top_num = 10, ig = NULL, displayBy = c("fc",
"adjp", "fdr", "zscore", "pvalue"), path.mode = c("all_paths",
"shortest_paths", "all_shortest_paths"), height = 7, width = 7, margin
= rep(0.1, 4), colormap = c("yr", "bwr", "jet", "gbr", "wyr", "br",
"rainbow", "wb", "lightyellow-orange"), ncolors = 40, zlim = NULL,
colorbar = T, colorbar.fraction = 0.1, newpage = T, layout.orientation
= c("top_bottom", "left_right", "bottom_top", "right_left"), node.info
= c("none", "term_id", "term_name", "both", "full_term_name"),
wrap.width = NULL, graph.node.attrs = NULL, graph.edge.attrs = NULL,
node.attrs = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{eTerm}{
%% ~~Describe \code{eTerm} here~~
}
\item{top_num}{
%% ~~Describe \code{top_num} here~~
}
\item{ig}{
%% ~~Describe \code{ig} here~~
}
\item{displayBy}{
%% ~~Describe \code{displayBy} here~~
}
\item{path.mode}{
%% ~~Describe \code{path.mode} here~~
}
\item{height}{
%% ~~Describe \code{height} here~~
}
\item{width}{
%% ~~Describe \code{width} here~~
}
\item{margin}{
%% ~~Describe \code{margin} here~~
}
\item{colormap}{
%% ~~Describe \code{colormap} here~~
}
\item{ncolors}{
%% ~~Describe \code{ncolors} here~~
}
\item{zlim}{
%% ~~Describe \code{zlim} here~~
}
\item{colorbar}{
%% ~~Describe \code{colorbar} here~~
}
\item{colorbar.fraction}{
%% ~~Describe \code{colorbar.fraction} here~~
}
\item{newpage}{
%% ~~Describe \code{newpage} here~~
}
\item{layout.orientation}{
%% ~~Describe \code{layout.orientation} here~~
}
\item{node.info}{
%% ~~Describe \code{node.info} here~~
}
\item{wrap.width}{
%% ~~Describe \code{wrap.width} here~~
}
\item{graph.node.attrs}{
%% ~~Describe \code{graph.node.attrs} here~~
}
\item{graph.edge.attrs}{
%% ~~Describe \code{graph.edge.attrs} here~~
}
\item{node.attrs}{
%% ~~Describe \code{node.attrs} here~~
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
function (eTerm, top_num = 10, ig = NULL, displayBy = c("fc",
"adjp", "fdr", "zscore", "pvalue"), path.mode = c("all_paths",
"shortest_paths", "all_shortest_paths"), height = 7, width = 7,
margin = rep(0.1, 4), colormap = c("yr", "bwr", "jet", "gbr",
"wyr", "br", "rainbow", "wb", "lightyellow-orange"),
ncolors = 40, zlim = NULL, colorbar = T, colorbar.fraction = 0.1,
newpage = T, layout.orientation = c("top_bottom", "left_right",
"bottom_top", "right_left"), node.info = c("none", "term_id",
"term_name", "both", "full_term_name"), wrap.width = NULL,
graph.node.attrs = NULL, graph.edge.attrs = NULL, node.attrs = NULL)
{
displayBy <- match.arg(displayBy)
path.mode <- match.arg(path.mode)
layout.orientation <- match.arg(layout.orientation)
node.info <- match.arg(node.info)
if (is.null(eTerm)) {
warnings("There is no enrichment in the 'eTerm' object.\n")
return(NULL)
}
if (class(eTerm)[1] == "eTerm") {
df <- xEnrichViewer(eTerm, top_num = "all")
if (top_num == "auto") {
top_num <- sum(df$adjp < 0.05)
if (top_num <= 1) {
top_num <- sum(df$adjp < 0.1)
}
if (top_num <= 1) {
top_num <- 10
}
}
df <- xEnrichViewer(eTerm, top_num = top_num, sortBy = "adjp")
nodes_query <- rownames(df)
if (class(ig) == "igraph") {
nodes_query_tmp <- intersect(nodes_query, V(ig)$name)
if (length(nodes_query_tmp) > 0) {
nodes_query <- nodes_query_tmp
}
}
g <- eTerm$g
subg <- dnet::dDAGinduce(g, nodes_query, path.mode = path.mode)
if (displayBy == "adjp" | displayBy == "fdr") {
data <- -1 * log10(df$adjp)
}
else if (displayBy == "fc") {
data <- df$fc
}
else if (displayBy == "pvalue") {
data <- -1 * log10(df$pvalue)
}
else if (displayBy == "zscore") {
data <- df$zscore
}
names(data) <- rownames(df)
data <- data[V(subg)$name]
}
else {
stop("Cannot find an 'eTerm' object.\n")
}
if (is.null(graph.edge.attrs)) {
}
if (is.null(node.attrs)) {
hightlighted.shape <- rep("box", length(nodes_query))
names(hightlighted.shape) <- V(subg)[nodes_query]$name
node.attrs <- list(shape = hightlighted.shape)
}
agDAG <- visDAG(g = subg, data = data, height = height, width = width,
margin = margin, colormap = colormap, ncolors = ncolors,
zlim = zlim, colorbar = colorbar, colorbar.fraction =
colorbar.fraction,
newpage = newpage, layout.orientation = layout.orientation,
node.info = node.info, numChar = wrap.width, graph.node.attrs =
graph.node.attrs,
graph.edge.attrs = graph.edge.attrs, node.attrs = node.attrs)
invisible(agDAG)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
