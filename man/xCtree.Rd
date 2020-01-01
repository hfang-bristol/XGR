\name{xCtree}
\alias{xCtree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xCtree(ig, leave.label.orientation = c("outwards", "inwards"),
leave.label.size = 2, leave.label.color = "steelblue",
leave.label.alpha = 0.7, leave.label.wrap = NULL, leave.label.expansion
= NULL, leave.size = 0, limit.expansion = 1.1, edge.color = "grey",
edge.alpha = 0.5, edge.width = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ig}{
%% ~~Describe \code{ig} here~~
}
\item{leave.label.orientation}{
%% ~~Describe \code{leave.label.orientation} here~~
}
\item{leave.label.size}{
%% ~~Describe \code{leave.label.size} here~~
}
\item{leave.label.color}{
%% ~~Describe \code{leave.label.color} here~~
}
\item{leave.label.alpha}{
%% ~~Describe \code{leave.label.alpha} here~~
}
\item{leave.label.wrap}{
%% ~~Describe \code{leave.label.wrap} here~~
}
\item{leave.label.expansion}{
%% ~~Describe \code{leave.label.expansion} here~~
}
\item{leave.size}{
%% ~~Describe \code{leave.size} here~~
}
\item{limit.expansion}{
%% ~~Describe \code{limit.expansion} here~~
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
function (ig, leave.label.orientation = c("outwards", "inwards"),
leave.label.size = 2, leave.label.color = "steelblue",
leave.label.alpha = 0.7,
leave.label.wrap = NULL, leave.label.expansion = NULL, leave.size = 0,
limit.expansion = 1.1, edge.color = "grey", edge.alpha = 0.5,
edge.width = 0.5)
{
leave.label.orientation <- match.arg(leave.label.orientation)
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
if (!is.null(leave.label.wrap)) {
width <- as.integer(leave.label.wrap)
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
if (!("tipid" \%in\% igraph::vertex_attr_names(ig))) {
ind <- match(V(ig)$name, dnet::dDAGtip(ig))
V(ig)$tipid <- NA
V(ig)$tipid[!is.na(ind)] <- 1:sum(!is.na(ind))
}
V(ig)$angle <- 90 - 360 * V(ig)$tipid/sum(!is.na(V(ig)$tipid))
if (leave.label.orientation == "outwards") {
V(ig)$hjust <- ifelse(V(ig)$angle < -90, 1, 0)
if (is.null(leave.label.expansion)) {
leave.label.expansion <- 1.05
}
}
else if (leave.label.orientation == "inwards") {
V(ig)$hjust <- ifelse(V(ig)$angle < -90, 0, 1)
if (is.null(leave.label.expansion)) {
leave.label.expansion <- 0.98
}
}
V(ig)$angle <- ifelse(V(ig)$angle < -90, V(ig)$angle + 180,
V(ig)$angle)
x <- y <- leaf <- label <- angle <- hjust <- NULL
gp <- ggraph::ggraph(ig, layout = "dendrogram", circular = TRUE)
gp <- gp + ggraph::geom_edge_diagonal(color = edge.color,
alpha = edge.alpha, width = edge.width)
gp <- gp + ggraph::geom_node_point(aes(filter = leaf), size =
leave.size,
color = edge.color, alpha = edge.alpha)
gp <- gp + ggraph::geom_node_text(aes(x = x * leave.label.expansion,
y = y * leave.label.expansion, filter = leaf, label = label,
angle = angle, hjust = hjust), show.legend = F, color =
leave.label.color,
size = leave.label.size, alpha = leave.label.alpha, fontface = "bold")
+
expand_limits(x = c(-limit.expansion, limit.expansion),
y = c(-limit.expansion, limit.expansion))
gp <- gp + coord_fixed() + theme(legend.position = "bottom") +
ggraph::theme_graph(base_family = "sans", plot_margin = margin(0,
0, 0, 0))
if (0) {
tipid <- NULL
}
if (1) {
ggraph.orig_index <- NULL
data <- gp$data \%>\% dplyr::arrange(ggraph.orig_index)
V(ig)$xcoord <- data[, "x"]
V(ig)$ycoord <- data[, "y"]
}
gp$ig <- ig
invisible(gp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
