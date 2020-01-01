\name{xSocialiserDAGplot}
\alias{xSocialiserDAGplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSocialiserDAGplot(g, query, displayBy = c("IC", "none"), path.mode =
c("all_paths", "shortest_paths", "all_shortest_paths"), height = 7,
width = 7, margin = rep(0.1, 4), colormap = c("yr", "bwr", "jet",
"gbr", "wyr", "br", "rainbow", "wb", "lightyellow-orange"), ncolors =
40, zlim = NULL, colorbar = T, colorbar.fraction = 0.1, newpage = T,
layout.orientation = c("top_bottom", "left_right", "bottom_top",
"right_left"), node.info = c("none", "term_id", "term_name", "both",
"full_term_name"), wrap.width = NULL, graph.node.attrs = NULL,
graph.edge.attrs = NULL, node.attrs = NULL, output.format =
c("Ragraph", "igraph"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{g}{
%% ~~Describe \code{g} here~~
}
\item{query}{
%% ~~Describe \code{query} here~~
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
\item{output.format}{
%% ~~Describe \code{output.format} here~~
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
function (g, query, displayBy = c("IC", "none"), path.mode =
c("all_paths",
"shortest_paths", "all_shortest_paths"), height = 7, width = 7,
margin = rep(0.1, 4), colormap = c("yr", "bwr", "jet", "gbr",
"wyr", "br", "rainbow", "wb", "lightyellow-orange"),
ncolors = 40, zlim = NULL, colorbar = T, colorbar.fraction = 0.1,
newpage = T, layout.orientation = c("top_bottom", "left_right",
"bottom_top", "right_left"), node.info = c("none", "term_id",
"term_name", "both", "full_term_name"), wrap.width = NULL,
graph.node.attrs = NULL, graph.edge.attrs = NULL, node.attrs = NULL,
output.format = c("Ragraph", "igraph"))
{
displayBy <- match.arg(displayBy)
path.mode <- match.arg(path.mode)
layout.orientation <- match.arg(layout.orientation)
node.info <- match.arg(node.info)
output.format <- match.arg(output.format)
if (is.logical(g)) {
stop("There is no similarity in the 'igraph' object.\n")
}
if (class(g) != "igraph") {
stop("The function must apply to the 'igraph' object.\n")
}
if (is.null(g$dag)) {
dag <- g
}
else {
dag <- g$dag
}
if (is.null(V(dag)$anno) | is.null(V(dag)$IC)) {
stop("The function requires that input graph has already contained
annotation data and also information content (IC).\n")
}
flag <- sapply(V(dag)$anno, function(x) {
ind <- match(query, x)
if (is.na(ind)) {
0
}
else {
if (is.null(names(x[ind]))) {
1
}
else if (names(x[ind]) == "o") {
1
}
else if (names(x[ind]) == "i") {
2
}
}
})
terms <- V(dag)$name[flag > 0]
terms_origin <- V(dag)$name[flag == 1]
if (length(terms) == 0) {
warnings("No terms are found to annotate the entity in query!")
return(NULL)
}
subg <- dnet::dDAGinduce(g = dag, nodes_query = terms, path.mode =
path.mode)
inherited <- rep(1, length(V(subg)$name))
names(inherited) <- V(subg)$name
if (length(terms_origin) > 0) {
inherited[terms_origin] <- 0
}
V(subg)$inherited <- inherited
if (displayBy == "IC") {
data <- V(subg)$IC
names(data) <- V(subg)$name
}
else {
data <- NULL
}
if (output.format == "igraph") {
invisible(subg)
}
else if (output.format == "Ragraph") {
if (is.null(graph.edge.attrs)) {
}
if (is.null(node.attrs)) {
if (length(terms_origin) > 0) {
hightlighted.shape <- rep("box", length(terms_origin))
names(hightlighted.shape) <- V(subg)[terms_origin]$name
node.attrs <- list(shape = hightlighted.shape)
}
}
agDAG <- visDAG(g = subg, data = data, height = height,
width = width, margin = margin, colormap = colormap,
ncolors = ncolors, zlim = zlim, colorbar = colorbar,
colorbar.fraction = colorbar.fraction, newpage = newpage,
layout.orientation = layout.orientation, node.info = node.info,
numChar = wrap.width, graph.node.attrs = graph.node.attrs,
graph.edge.attrs = graph.edge.attrs, node.attrs = node.attrs)
invisible(agDAG)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
