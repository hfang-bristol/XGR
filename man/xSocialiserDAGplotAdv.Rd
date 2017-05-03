\name{xSocialiserDAGplotAdv}
\alias{xSocialiserDAGplotAdv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSocialiserDAGplotAdv(g, query1, query2, displayBy = c("IC", "none"),
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
height = 7, width = 7, margin = rep(0.1, 4), colormap = c("wyr", "bwr",
"jet", "gbr", "yr", "br", "rainbow", "wb", "lightyellow-orange"),
ncolors = 40, zlim = NULL, colorbar = T, colorbar.fraction = 0.1,
newpage = T, layout.orientation = c("top_bottom", "left_right",
"bottom_top", "right_left"), node.info = c("term_name", "term_id"),
wrap.width = NULL, graph.node.attrs = NULL, graph.edge.attrs = NULL,
node.attrs = NULL, output.format = c("Ragraph", "igraph"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{g}{
%% ~~Describe \code{g} here~~
}
\item{query1}{
%% ~~Describe \code{query1} here~~
}
\item{query2}{
%% ~~Describe \code{query2} here~~
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
function (g, query1, query2, displayBy = c("IC", "none"), path.mode =
c("all_paths",
"shortest_paths", "all_shortest_paths"), height = 7, width = 7,
margin = rep(0.1, 4), colormap = c("wyr", "bwr", "jet", "gbr",
"yr", "br", "rainbow", "wb", "lightyellow-orange"), ncolors = 40,
zlim = NULL, colorbar = T, colorbar.fraction = 0.1, newpage = T,
layout.orientation = c("top_bottom", "left_right", "bottom_top",
"right_left"), node.info = c("term_name", "term_id"),
wrap.width = NULL, graph.node.attrs = NULL, graph.edge.attrs = NULL,
node.attrs = NULL, output.format = c("Ragraph", "igraph"))
{
displayBy <- match.arg(displayBy)
path.mode <- match.arg(path.mode)
layout.orientation <- match.arg(layout.orientation)
node.info <- match.arg(node.info)
output.format <- match.arg(output.format)
ig_SNP1 <- xSocialiserDAGplot(g, query = query1, output.format =
"igraph")
ig_SNP2 <- xSocialiserDAGplot(g, query = query2, output.format =
"igraph")
if (is.null(ig_SNP1) | is.null(ig_SNP2)) {
warnings("No terms are found to annotate one or two entities in
query!")
return(NULL)
}
list_igraph <- list(ig_SNP1, ig_SNP2)
names(list_igraph) <- c(query1, query2)
ls_edges <- lapply(list_igraph, function(x) {
df_edge <- igraph::get.data.frame(x, what = "edges")
})
relations <- do.call(rbind, ls_edges)
relations <- relations[!duplicated(relations), ]
ls_nodes <- lapply(list_igraph, function(x) {
df_nodes <- igraph::get.data.frame(x, what = "vertices")
})
nodes <- do.call(rbind, ls_nodes)
nodes <- nodes[, -7]
nodes <- nodes[!duplicated(nodes), ]
ig <- igraph::graph.data.frame(d = relations, directed = T,
vertices = nodes)
ls_res <- lapply(list_igraph, function(x) {
ind <- match(V(ig)$name, V(x)$name)
res <- V(x)$inherited[ind] * (-1) + 2
res[is.na(res)] <- 0
res
})
df_res <- do.call(cbind, ls_res)
code <- apply(df_res, 1, paste, collapse = "-")
V(ig)$code <- code
subg <- ig
V(subg)$name <- paste(V(subg)$code, V(subg)$name, sep = "\\\n")
V(subg)$term_name <- paste(V(subg)$code, V(subg)$term_name,
sep = "\\\n")
if (displayBy == "IC") {
data <- V(subg)$IC
names(data) <- V(subg)$name
}
else {
data <- NULL
}
if (output.format == "igraph") {
invisible(ig)
}
else if (output.format == "Ragraph") {
if (is.null(graph.edge.attrs)) {
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
