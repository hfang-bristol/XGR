\name{xGGraph}
\alias{xGGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGGraph(ig, layout = "partition", circular = T, leave = T,
node.label.size = 2, node.label.direction = c("none", "leftright",
"topbottom"), node.label.color = "steelblue", node.label.alpha = 0.7,
node.label.wrap = NULL, node.label.offset = 0.5, node.size = 2,
limit.expansion = NULL, edge = c("diagonal", "link", "arc", "fan",
"elbow"), edge.color = "grey", edge.alpha = 0.5, edge.width = 0.5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ig}{
%% ~~Describe \code{ig} here~~
}
\item{layout}{
%% ~~Describe \code{layout} here~~
}
\item{circular}{
%% ~~Describe \code{circular} here~~
}
\item{leave}{
%% ~~Describe \code{leave} here~~
}
\item{node.label.size}{
%% ~~Describe \code{node.label.size} here~~
}
\item{node.label.direction}{
%% ~~Describe \code{node.label.direction} here~~
}
\item{node.label.color}{
%% ~~Describe \code{node.label.color} here~~
}
\item{node.label.alpha}{
%% ~~Describe \code{node.label.alpha} here~~
}
\item{node.label.wrap}{
%% ~~Describe \code{node.label.wrap} here~~
}
\item{node.label.offset}{
%% ~~Describe \code{node.label.offset} here~~
}
\item{node.size}{
%% ~~Describe \code{node.size} here~~
}
\item{limit.expansion}{
%% ~~Describe \code{limit.expansion} here~~
}
\item{edge}{
%% ~~Describe \code{edge} here~~
}
\item{edge.color}{
%% ~~Describe \code{edge.color} here~~
}
\item{edge.alpha}{
%% ~~Describe \code{edge.alpha} here~~
}
\item{edge.width}{
%% ~~Describe \code{edge.width} here~~
}
\item{\dots}{
%% ~~Describe \code{\dots} here~~
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
function (ig, layout = "partition", circular = T, leave = T,
node.label.size = 2, node.label.direction = c("none", "leftright",
"topbottom"), node.label.color = "steelblue", node.label.alpha = 0.7,
node.label.wrap = NULL, node.label.offset = 0.5, node.size = 2,
limit.expansion = NULL, edge = c("diagonal", "link", "arc",
"fan", "elbow"), edge.color = "grey", edge.alpha = 0.5,
edge.width = 0.5, ...)
{
node.label.direction <- match.arg(node.label.direction)
edge <- match.arg(edge)
if (0) {
if (any(duplicated(tree$node.label)))
tree$node.label <- NULL
ig <- as.igraph(tree, directed = T, use.labels = T)
}
if (class(ig) != "igraph") {
warnings("The function must apply to the 'igraph' object.\n")
return(NULL)
}
else {
if (!all(c("name") \%in\% igraph::vertex_attr_names(ig))) {
warnings("The igraph object must have vertex attribute 'name'.\n")
return(NULL)
}
}
if (!("label" \%in\% igraph::vertex_attr_names(ig))) {
V(ig)$label <- V(ig)$name
}
if (!is.null(node.label.wrap)) {
width <- as.integer(node.label.wrap)
res_list <- lapply(V(ig)$label, function(x) {
x <- gsub("_", " ", x)
y <- strwrap(x, width = width)
if (length(y) > 1) {
paste(y, collapse = "\n")
}
else {
y
}
})
V(ig)$label <- unlist(res_list)
}
x <- y <- leaf <- label <- name <- NULL
if (length(suppressWarnings(tryCatch(gp <- ggraph::ggraph(ig,
layout = layout, circular = circular), error = function(e) e,
warning = function(w) w))) == 2) {
warning("The layout does not support circularity!")
circular <- F
gp <- ggraph::ggraph(ig, layout = layout, circular = circular)
}
if (leave & !is.null(gp$data$leaf)) {
gp <- gp + ggraph::geom_node_point(aes(filter = leaf),
size = node.size, color = edge.color, alpha = edge.alpha)
}
else {
gp <- gp + ggraph::geom_node_point(size = node.size,
color = edge.color, alpha = edge.alpha)
}
if (edge == "diagonal") {
gp <- gp + ggraph::geom_edge_diagonal(color = edge.color,
alpha = edge.alpha, width = edge.width)
}
else if (edge == "link") {
gp <- gp + ggraph::geom_edge_link(color = edge.color,
alpha = edge.alpha, width = edge.width)
}
else if (edge == "arc") {
gp <- gp + ggraph::geom_edge_arc(color = edge.color,
alpha = edge.alpha, width = edge.width)
}
else if (edge == "fan") {
gp <- gp + ggraph::geom_edge_fan(color = edge.color,
alpha = edge.alpha, width = edge.width)
}
else if (edge == "elbow") {
gp <- gp + ggraph::geom_edge_elbow(color = edge.color,
alpha = edge.alpha, width = edge.width)
}
if (node.label.size > 0) {
if (leave & !is.null(gp$data$leaf)) {
df <- subset(gp$data, leaf == T)
}
else {
df <- gp$data
}
if (node.label.direction == "none") {
gp <- gp + ggrepel::geom_text_repel(data = df, aes(x = x,
y = y, label = label), color = node.label.color,
size = node.label.size, alpha = node.label.alpha,
show.legend = F, segment.alpha = 0.5, segment.color = "grey50",
segment.size = 0.2, arrow = arrow(length = unit(0.01,
"npc")), ...)
}
else if (node.label.direction == "leftright") {
offset <- (range(gp$data$x)[2] - range(gp$data$x)[1]) *
node.label.offset
root <- subset(gp$data, name == dnet::dDAGroot(ig))
df1 <- subset(df, x < root$x)
df1$nudge_x <- -1 * (df1$x - min(gp$data$x)) - offset
gp <- gp + ggrepel::geom_text_repel(data = df1, aes(x = x,
y = y, label = label), color = node.label.color,
size = node.label.size, alpha = node.label.alpha,
show.legend = F, segment.alpha = 0.5, segment.color = "grey50",
segment.size = 0.2, arrow = arrow(length = unit(0.01,
"npc")), direction = "y", hjust = 0, nudge_x = df1$nudge_x)
df2 <- subset(df, x >= root$x)
df2$nudge_x <- (max(gp$data$x) - df2$x) + offset
gp <- gp + ggrepel::geom_text_repel(data = df2, aes(x = x,
y = y, label = label), color = node.label.color,
size = node.label.size, alpha = node.label.alpha,
show.legend = F, segment.alpha = 0.5, segment.color = "grey50",
segment.size = 0.2, arrow = arrow(length = unit(0.01,
"npc")), direction = "y", hjust = 1, nudge_x = df2$nudge_x)
if (is.null(limit.expansion)) {
gp <- gp + expand_limits(x = range(gp$data$x) *
(1 + 2 * node.label.offset), y = range(gp$data$y))
}
else {
gp <- gp + expand_limits(x = c(-limit.expansion,
limit.expansion), y = c(-limit.expansion, limit.expansion))
}
}
else if (node.label.direction == "topbottom") {
offset <- (range(gp$data$y)[2] - range(gp$data$y)[1]) *
node.label.offset
root <- subset(gp$data, name == dnet::dDAGroot(ig))
df1 <- subset(df, y < root$y)
df1$nudge_y <- -1 * (df1$y - min(gp$data$y)) - offset
gp <- gp + ggrepel::geom_text_repel(data = df1, aes(x = x,
y = y, label = label), color = node.label.color,
size = node.label.size, alpha = node.label.alpha,
show.legend = F, segment.alpha = 0.5, segment.color = "grey50",
segment.size = 0.2, arrow = arrow(length = unit(0.01,
"npc")), direction = "x", hjust = 0, nudge_y = df1$nudge_y,
angle = 90)
df2 <- subset(df, y >= root$y)
df2$nudge_y <- (max(gp$data$y) - df2$y) + offset
gp <- gp + ggrepel::geom_text_repel(data = df2, aes(x = x,
y = y, label = label), color = node.label.color,
size = node.label.size, alpha = node.label.alpha,
show.legend = F, segment.alpha = 0.5, segment.color = "grey50",
segment.size = 0.2, arrow = arrow(length = unit(0.01,
"npc")), direction = "x", hjust = 1, nudge_y = df2$nudge_y,
angle = 90)
if (is.null(limit.expansion)) {
gp <- gp + expand_limits(x = range(gp$data$x),
y = range(gp$data$y) * (1 + 2 * node.label.offset))
}
else {
gp <- gp + expand_limits(x = c(-limit.expansion,
limit.expansion), y = c(-limit.expansion, limit.expansion))
}
}
}
gp <- gp + ggraph::theme_graph(base_family = "sans", plot_margin =
margin(0,
0, 0, 0))
if (0) {
tipid <- NULL
}
gp$ig <- ig
invisible(gp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
