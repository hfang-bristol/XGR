\name{xGraphML}
\alias{xGraphML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGraphML(g, node.label = NULL, label.wrap.width = NULL, node.label.size
= 12, node.label.color = "#000000", node.tooltip = NULL, node.link =
NULL, node.xcoord = "xcoord", node.ycoord = "ycoord", node.color.na =
"#dddddd", node.color = NULL, colormap = "grey-orange-darkred", ncolors
= 64, nlegend = 11, legend.label.size = 10, legend.interval = 0.05,
zlim = NULL, node.size = 30, node.coord.scale = 300, edge.color =
"#00000033", edge.width = 1, filename = "xGraphML", verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{g}{
%% ~~Describe \code{g} here~~
}
\item{node.label}{
%% ~~Describe \code{node.label} here~~
}
\item{label.wrap.width}{
%% ~~Describe \code{label.wrap.width} here~~
}
\item{node.label.size}{
%% ~~Describe \code{node.label.size} here~~
}
\item{node.label.color}{
%% ~~Describe \code{node.label.color} here~~
}
\item{node.tooltip}{
%% ~~Describe \code{node.tooltip} here~~
}
\item{node.link}{
%% ~~Describe \code{node.link} here~~
}
\item{node.xcoord}{
%% ~~Describe \code{node.xcoord} here~~
}
\item{node.ycoord}{
%% ~~Describe \code{node.ycoord} here~~
}
\item{node.color.na}{
%% ~~Describe \code{node.color.na} here~~
}
\item{node.color}{
%% ~~Describe \code{node.color} here~~
}
\item{colormap}{
%% ~~Describe \code{colormap} here~~
}
\item{ncolors}{
%% ~~Describe \code{ncolors} here~~
}
\item{nlegend}{
%% ~~Describe \code{nlegend} here~~
}
\item{legend.label.size}{
%% ~~Describe \code{legend.label.size} here~~
}
\item{legend.interval}{
%% ~~Describe \code{legend.interval} here~~
}
\item{zlim}{
%% ~~Describe \code{zlim} here~~
}
\item{node.size}{
%% ~~Describe \code{node.size} here~~
}
\item{node.coord.scale}{
%% ~~Describe \code{node.coord.scale} here~~
}
\item{edge.color}{
%% ~~Describe \code{edge.color} here~~
}
\item{edge.width}{
%% ~~Describe \code{edge.width} here~~
}
\item{filename}{
%% ~~Describe \code{filename} here~~
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
function (g, node.label = NULL, label.wrap.width = NULL,
node.label.size = 12,
node.label.color = "#000000", node.tooltip = NULL, node.link = NULL,
node.xcoord = "xcoord", node.ycoord = "ycoord", node.color.na =
"#dddddd",
node.color = NULL, colormap = "grey-orange-darkred", ncolors = 64,
nlegend = 11, legend.label.size = 10, legend.interval = 0.05,
zlim = NULL, node.size = 30, node.coord.scale = 300, edge.color =
"#00000033",
edge.width = 1, filename = "xGraphML", verbose = T)
{
if (class(g) != "igraph") {
stop("The function must apply to either 'igraph' or 'graphNEL'
object.\n")
}
else {
ig <- g
}
df_legends <- NULL
if (!is.null(node.color)) {
pattern <- igraph::vertex_attr(ig, node.color)
if (!is.null(pattern)) {
pattern <- as.numeric(pattern)
pattern_nona <- pattern[!is.na(pattern)]
pattern_nona <- as.numeric(pattern_nona)
if (is.null(zlim)) {
vmin <- floor(stats::quantile(pattern_nona, 0.05))
vmax <- ceiling(stats::quantile(pattern_nona,
0.95))
if (vmin < 0 & vmax > 0) {
vsym <- abs(min(vmin, vmax))
vmin <- -1 * vsym
vmax <- vsym
}
zlim <- c(vmin, vmax)
}
vec2color <- function(vec, colormap = colormap, ncolors = ncolors,
zlim = zlim) {
palette.name <- xColormap(colormap = colormap)
colors <- palette.name(ncolors)
scale <- length(colors)/(max(zlim) - min(zlim))
sapply(1:length(vec), function(x) {
if (is.na(vec[x])) {
"#dddddd"
}
else {
ind <- floor(1 + (vec[x] - min(zlim)) * scale)
colors[max(1, min(ncolors, ind))]
}
})
}
node.color <- vec2color(pattern, colormap = colormap,
ncolors = ncolors, zlim = zlim)
colors <- xColormap(colormap = colormap)(ncolors)
legend_colors <- colors[round(seq(1, ncolors, length.out = nlegend))]
df_legends <- data.frame(name = paste0("l", 1:nlegend),
colors = legend_colors, labels = signif(seq(min(zlim),
max(zlim), length.out = nlegend), digits = 2),
x = rep(-1.1, nlegend), y = 1 - cumsum(c(0, rep(legend.interval,
nlegend - 1))), stringsAsFactors = FALSE)
}
else {
warning("The input 'pattern' is ignored. Please check the help for
enabling your input")
node.color <- rep("#BFFFBF", vcount(ig))
}
}
else {
node.color <- rep("#BFFFBF", vcount(ig))
}
output.head <- "<?xml version=\"1.0\" encoding=\"UTF-8\"
standalone=\"no\"?>\n<graphml
xmlns=\"http://graphml.graphdrawing.org/xmlns\"
xmlns:java=\"http://www.yworks.com/xml/yfiles-common/1.0/java\"
xmlns:sys=\"http://www.yworks.com/xml/yfiles-common/markup/primitives/2.0\"
xmlns:x=\"http://www.yworks.com/xml/yfiles-common/markup/2.0\"
xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"
xmlns:y=\"http://www.yworks.com/xml/graphml\"
xmlns:yed=\"http://www.yworks.com/xml/yed/3\"
xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns
http://www.yworks.com/xml/schema/graphml/1.1/ygraphml.xsd\">\n
<!--Created by yEd 3.17-->\n <key for=\"node\" id=\"d1\"
attr.name=\"url\" attr.type=\"string\"/>\n <key for=\"node\" id=\"d2\"
attr.name=\"description\" attr.type=\"string\"/>\n <key for=\"node\"
id=\"d3\" yfiles.type=\"nodegraphics\"/>\n <key for=\"edge\" id=\"d4\"
attr.name=\"description\" attr.type=\"string\"/>\n <key for=\"edge\"
id=\"d5\" yfiles.type=\"edgegraphics\"/>\n <key for=\"graphml\"
id=\"d6\" yfiles.type=\"resources\"/>\n <graph edgedefault=\"directed\"
id=\"G\">"
nnode <- igraph::vcount(ig)
if (length(node.label) != nnode) {
if (!is.null(node.label)) {
node.label <- igraph::vertex_attr(ig, node.label)
}
if (is.null(node.label)) {
node.label <- rep(NA, nnode)
}
}
node.label <- unlist(lapply(node.label, function(x) gsub("/",
"-", x)))
node.label <- unlist(lapply(node.label, function(x) gsub("&",
"-", x)))
if (!is.null(label.wrap.width)) {
width <- as.integer(label.wrap.width)
res_list <- lapply(node.label, function(x) {
if (!is.na(x)) {
x <- gsub("_", " ", x)
y <- strwrap(x, width = width)
if (length(y) == 2) {
paste(y, collapse = "\n")
}
else if (length(y) > 2) {
paste(y, collapse = "\n")
}
else {
y
}
}
else {
x
}
})
node.label <- unlist(res_list)
}
if (length(node.tooltip) != nnode) {
if (!is.null(node.tooltip)) {
node.tooltip <- igraph::vertex_attr(ig, node.tooltip)
}
if (is.null(node.tooltip)) {
node.tooltip <- igraph::vertex_attr(ig, "name")
}
}
node.tooltip <- unlist(lapply(node.tooltip, function(x) gsub("/",
"-", x)))
node.tooltip <- unlist(lapply(node.tooltip, function(x) gsub("&",
"-", x)))
if (length(node.link) != nnode) {
if (!is.null(node.link)) {
node.link <- igraph::vertex_attr(ig, node.link)
}
if (is.null(node.link)) {
node.link <- igraph::vertex_attr(ig, "name")
}
}
if (length(node.size) != nnode) {
if (!is.null(node.size)) {
tmp.node.size <- igraph::vertex_attr(ig, node.size)
}
else {
tmp.node.size <- rep(30, nnode)
}
if (is.null(tmp.node.size)) {
node.size <- rep(node.size, nnode)
}
else {
node.size <- tmp.node.size
}
}
V(ig)$label <- V(ig)$name
V(ig)$name <- paste0("n", 1:vcount(ig))
if (length(node.xcoord) != nnode | length(node.ycoord) !=
nnode) {
if (!is.null(node.xcoord)) {
node.xcoord <- igraph::vertex_attr(ig, node.xcoord)
}
if (!is.null(node.ycoord)) {
node.ycoord <- igraph::vertex_attr(ig, node.ycoord)
}
if (is.null(node.xcoord) | is.null(node.ycoord)) {
glayout <- igraph::layout_as_tree(ig, root = dnet::dDAGroot(ig),
circular = TRUE, flip.y = TRUE)
glayout <- glayout[, c(2:1)]
node.xcoord <- glayout[, 1]
node.ycoord <- glayout[, 2]
}
}
if (max(node.xcoord) != min(node.xcoord)) {
node.xcoord <- (node.xcoord - min(node.xcoord))/(max(node.xcoord) -
min(node.xcoord)) * 2 - 1
}
if (max(node.ycoord) != min(node.ycoord)) {
node.ycoord <- (node.ycoord - min(node.ycoord))/(max(node.ycoord) -
min(node.ycoord)) * 2 - 1
}
node.ycoord <- -node.ycoord
node.xcoord <- node.xcoord * node.coord.scale
node.ycoord <- node.ycoord * node.coord.scale
df_nodes <- igraph::get.data.frame(ig, what = "vertices")
df_nodes$node.label <- node.label
df_nodes$node.tooltip <- node.tooltip
df_nodes$node.link <- node.link
df_nodes$node.color <- paste0(node.color, "cc")
df_nodes$node.size <- node.size
ls_nodes <- lapply(1:nrow(df_nodes), function(i) {
k <- 0
vec <- vector()
k <- k + 1
vec[k] <- paste0("<node id=\"", df_nodes$name[i], "\">")
if (!is.null(node.link)) {
k <- k + 1
vec[k] <- paste0("<data key=\"d1\"><![CDATA[", df_nodes$node.link[i],
"]]></data>")
}
k <- k + 1
vec[k] <- paste0("<data key=\"d2\"><![CDATA[",
df_nodes$node.tooltip[i],
"]]></data>")
k <- k + 1
vec[k] <- paste0("<data key=\"d3\">")
k <- k + 1
vec[k] <- paste0("<y:ShapeNode>")
k <- k + 1
vec[k] <- paste0("<y:Geometry height=\"", df_nodes$node.size[i],
"\" width=\"", df_nodes$node.size[i], "\" x=\"",
node.xcoord[i], "\" y=\"", node.ycoord[i], "\"/>")
k <- k + 1
vec[k] <- paste0("<y:Fill color=\"", df_nodes$node.color[i],
"\" transparent=\"false\"/>")
k <- k + 1
vec[k] <- paste0("<y:BorderStyle color=\"#dddddd\" raised=\"false\"
type=\"line\" width=\"1\"/>")
k <- k + 1
if (!is.na(df_nodes$node.label[i])) {
vec[k] <- paste0("<y:NodeLabel alignment=\"center\"
autoSizePolicy=\"content\" borderDistance=\"0.0\" fontFamily=\"Arial\"
fontSize=\"",
node.label.size, "\" fontStyle=\"italic\" hasBackgroundColor=\"false\"
hasLineColor=\"false\" height=\"30\" horizontalTextPosition=\"center\"
iconTextGap=\"4\" modelName=\"sides\" modelPosition=\"n\"
textColor=\"",
node.label.color, "\" verticalTextPosition=\"bottom\" visible=\"true\"
width=\"30\" x=\"0\" y=\"0\">",
df_nodes$node.label[i], "</y:NodeLabel>")
}
else {
vec[k] <- paste0("<y:NodeLabel alignment=\"center\"
autoSizePolicy=\"content\" borderDistance=\"0.0\" fontFamily=\"Arial\"
fontSize=\"",
node.label.size, "\" fontStyle=\"italic\" hasBackgroundColor=\"false\"
hasLineColor=\"false\" height=\"30\" horizontalTextPosition=\"center\"
iconTextGap=\"4\" modelName=\"sides\" modelPosition=\"n\"
textColor=\"",
node.label.color, "\" verticalTextPosition=\"bottom\" visible=\"false\"
width=\"30\" x=\"0\" y=\"0\">",
df_nodes$node.tooltip[i], "</y:NodeLabel>")
}
k <- k + 1
vec[k] <- paste0("<y:Shape type=\"ellipse\"/>")
k <- k + 1
vec[k] <- paste0("</y:ShapeNode>")
k <- k + 1
vec[k] <- paste0("</data>")
k <- k + 1
vec[k] <- paste0("</node>")
paste(vec, collapse = "\n")
})
vec_nodes <- unlist(ls_nodes)
output.nodes <- paste(vec_nodes, collapse = "\n")
if (!is.null(df_legends)) {
ls_legends <- lapply(1:nrow(df_legends), function(i) {
k <- 0
vec <- vector()
k <- k + 1
vec[k] <- paste0("<node id=\"", df_legends$name[i],
"\">")
k <- k + 1
vec[k] <- paste0("<data key=\"d2\"><![CDATA[", df_legends$labels[i],
"]]></data>")
k <- k + 1
vec[k] <- paste0("<data key=\"d3\">")
k <- k + 1
vec[k] <- paste0("<y:ShapeNode>")
k <- k + 1
vec[k] <- paste0("<y:Geometry height=\"", legend.interval *
node.coord.scale, "\" width=\"", legend.interval *
node.coord.scale, "\" x=\"", df_legends$x[i] *
node.coord.scale, "\" y=\"", df_legends$y[i] *
node.coord.scale, "\"/>")
k <- k + 1
vec[k] <- paste0("<y:Fill color=\"", df_legends$colors[i],
"\" transparent=\"false\"/>")
k <- k + 1
vec[k] <- paste0("<y:BorderStyle color=\"#dddddd\" raised=\"false\"
type=\"line\" width=\"1\"/>")
k <- k + 1
vec[k] <- paste0("<y:NodeLabel alignment=\"center\"
autoSizePolicy=\"content\" borderDistance=\"0.0\" fontFamily=\"Arial\"
fontSize=\"",
legend.label.size, "\" fontStyle=\"plain\" hasBackgroundColor=\"false\"
hasLineColor=\"false\" height=\"30\" horizontalTextPosition=\"center\"
iconTextGap=\"4\" modelName=\"sides\" modelPosition=\"w\"
textColor=\"#000000\" verticalTextPosition=\"bottom\" visible=\"true\"
width=\"30\" x=\"0\" y=\"0\">",
df_legends$labels[i], "</y:NodeLabel>")
k <- k + 1
vec[k] <- paste0("<y:Shape type=\"ellipse\"/>")
k <- k + 1
vec[k] <- paste0("</y:ShapeNode>")
k <- k + 1
vec[k] <- paste0("</data>")
k <- k + 1
vec[k] <- paste0("</node>")
paste(vec, collapse = "\n")
})
vec_legends <- unlist(ls_legends)
output.legends <- paste(vec_legends, collapse = "\n")
}
else {
output.legends <- NULL
}
df_edges <- igraph::get.data.frame(ig, what = "edges")
ls_edges <- lapply(1:nrow(df_edges), function(i) {
source <- df_edges$from[i]
target <- df_edges$to[i]
k <- 0
vec <- vector()
k <- k + 1
vec[k] <- paste0("<edge id=\"", "e", i, "\" source=\"",
source, "\" target=\"", target, "\">")
k <- k + 1
vec[k] <- paste0("<data key=\"d5\">")
k <- k + 1
vec[k] <- paste0("<y:GenericEdge configuration=\"DEFAULT\">")
if (0) {
k <- k + 1
vec[k] <- paste0("<y:Path sx=\"0\" sy=\"0\" tx=\"0\" ty=\"0\"/>")
}
k <- k + 1
vec[k] <- paste0("<y:LineStyle color=\"", edge.color,
"\" type=\"line\" width=\"", edge.width, "\"/>")
k <- k + 1
if (igraph::is_directed(ig)) {
vec[k] <- paste0("<y:Arrows source=\"none\" target=\"standard\"/>")
}
else {
vec[k] <- paste0("<y:Arrows source=\"none\" target=\"none\"/>")
}
k <- k + 1
vec[k] <- paste0("</y:GenericEdge>")
k <- k + 1
vec[k] <- paste0("</data>")
k <- k + 1
vec[k] <- paste0("</edge>")
paste(vec, collapse = "\n")
})
vec_edges <- unlist(ls_edges)
output.edges <- paste(vec_edges, collapse = "\n")
output.tail <- "</graph>\n <data key=\"d6\">\n <y:Resources/>\n
</data>\n</graphml>"
if (is.null(output.legends)) {
output <- paste0(output.head, "\n", output.nodes, "\n",
output.edges, "\n", output.tail, "\n")
}
else {
output <- paste0(output.head, "\n", output.nodes, "\n",
output.legends, "\n", output.edges, "\n", output.tail,
"\n")
}
if (!is.null(filename)) {
filename <- gsub(".graphml$", "", filename)
outputfile <- paste0(filename, ".graphml")
fileConn <- base::file(outputfile)
base::writeLines(output, fileConn)
base::close(fileConn)
if (verbose) {
message(sprintf("Congratulations! A file '\%s' (in the directory \%s)
has been created!",
outputfile, getwd()), appendLF = T)
}
}
invisible(output)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
