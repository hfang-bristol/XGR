\name{xGGnetwork}
\alias{xGGnetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGGnetwork(g, node.label = NULL, label.wrap.width = NULL,
label.wrap.lineheight = 0.8, node.label.size = NULL,
node.label.fontface = "plain", node.label.color = "darkblue",
node.label.alpha = 0.8, node.label.padding = 1, node.label.arrow =
0.01, node.label.force = 1, node.shape = 19, node.shape.title = NULL,
node.xcoord = NULL, node.ycoord = NULL, node.color = NULL,
node.color.title = NULL, colormap = "grey-orange-darkred", ncolors =
64, zlim = NULL, na.color = "grey80", node.color.alpha = 1, node.size =
NULL, node.size.title = NULL, node.size.range = c(1, 4), slim = NULL,
title = "", edge.size = 0.5, edge.color = "black", edge.color.alpha =
0.5, edge.curve = 0.1, edge.arrow = 2, edge.arrow.gap = 0.02, ncolumns
= NULL)
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
\item{label.wrap.lineheight}{
%% ~~Describe \code{label.wrap.lineheight} here~~
}
\item{node.label.size}{
%% ~~Describe \code{node.label.size} here~~
}
\item{node.label.fontface}{
%% ~~Describe \code{node.label.fontface} here~~
}
\item{node.label.color}{
%% ~~Describe \code{node.label.color} here~~
}
\item{node.label.alpha}{
%% ~~Describe \code{node.label.alpha} here~~
}
\item{node.label.padding}{
%% ~~Describe \code{node.label.padding} here~~
}
\item{node.label.arrow}{
%% ~~Describe \code{node.label.arrow} here~~
}
\item{node.label.force}{
%% ~~Describe \code{node.label.force} here~~
}
\item{node.shape}{
%% ~~Describe \code{node.shape} here~~
}
\item{node.shape.title}{
%% ~~Describe \code{node.shape.title} here~~
}
\item{node.xcoord}{
%% ~~Describe \code{node.xcoord} here~~
}
\item{node.ycoord}{
%% ~~Describe \code{node.ycoord} here~~
}
\item{node.color}{
%% ~~Describe \code{node.color} here~~
}
\item{node.color.title}{
%% ~~Describe \code{node.color.title} here~~
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
\item{na.color}{
%% ~~Describe \code{na.color} here~~
}
\item{node.color.alpha}{
%% ~~Describe \code{node.color.alpha} here~~
}
\item{node.size}{
%% ~~Describe \code{node.size} here~~
}
\item{node.size.title}{
%% ~~Describe \code{node.size.title} here~~
}
\item{node.size.range}{
%% ~~Describe \code{node.size.range} here~~
}
\item{slim}{
%% ~~Describe \code{slim} here~~
}
\item{title}{
%% ~~Describe \code{title} here~~
}
\item{edge.size}{
%% ~~Describe \code{edge.size} here~~
}
\item{edge.color}{
%% ~~Describe \code{edge.color} here~~
}
\item{edge.color.alpha}{
%% ~~Describe \code{edge.color.alpha} here~~
}
\item{edge.curve}{
%% ~~Describe \code{edge.curve} here~~
}
\item{edge.arrow}{
%% ~~Describe \code{edge.arrow} here~~
}
\item{edge.arrow.gap}{
%% ~~Describe \code{edge.arrow.gap} here~~
}
\item{ncolumns}{
%% ~~Describe \code{ncolumns} here~~
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
label.wrap.lineheight = 0.8,
node.label.size = NULL, node.label.fontface = "plain", node.label.color
= "darkblue",
node.label.alpha = 0.8, node.label.padding = 1, node.label.arrow =
0.01,
node.label.force = 1, node.shape = 19, node.shape.title = NULL,
node.xcoord = NULL, node.ycoord = NULL, node.color = NULL,
node.color.title = NULL, colormap = "grey-orange-darkred",
ncolors = 64, zlim = NULL, na.color = "grey80", node.color.alpha = 1,
node.size = NULL, node.size.title = NULL, node.size.range = c(1,
4), slim = NULL, title = "", edge.size = 0.5, edge.color = "black",
edge.color.alpha = 0.5, edge.curve = 0.1, edge.arrow = 2,
edge.arrow.gap = 0.02, ncolumns = NULL)
{
if (any(class(g) \%in\% c("igraph"))) {
ls_ig <- list(g)
}
else if (class(g) == "list") {
ls_ig <- base::Filter(base::Negate(is.null), g)
if (length(ls_ig) == 0) {
return(NULL)
}
}
else {
stop("The function must apply to 'list' of 'igraph' objects or a
'igraph' object.\n")
}
ls_names <- names(ls_ig)
if (is.null(ls_names)) {
ls_names <- paste("IG", 1:length(ls_ig), sep = "_")
names(ls_ig) <- ls_names
}
ls_df <- lapply(1:length(ls_ig), function(i) {
ig <- ls_ig[[i]]
if (igraph::vcount(ig) == 0) {
return(NULL)
}
if (is.null(V(ig)$name)) {
V(ig)$name <- 1:vcount(ig)
}
node_attrs <- igraph::vertex_attr_names(ig)
for (k in 1:length(node_attrs)) {
if (class(igraph::vertex_attr(ig, node_attrs[k])) ==
"list") {
ig <- ig \%>\% igraph::delete_vertex_attr(node_attrs[k])
}
}
nnode <- igraph::vcount(ig)
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
if (all(is.na(glayout))) {
glayout <- igraph::layout_with_kk(ig)
}
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
if (length(node.label) != nnode) {
if (!is.null(node.label)) {
node.label <- igraph::vertex_attr(ig, node.label)
}
if (is.null(node.label)) {
node.label <- rep("", nnode)
}
}
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
paste0(paste(y[1:2], collapse = "\n"), "...")
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
V(ig)$n.label <- node.label
if (length(node.label.size) != nnode) {
if (!is.null(node.label.size)) {
tmp.node.label.size <- igraph::vertex_attr(ig,
node.label.size)
}
else {
tmp.node.label.size <- rep(0, nnode)
}
if (is.null(tmp.node.label.size)) {
node.label.size <- rep(node.label.size, nnode)
}
else {
node.label.size <- tmp.node.label.size
}
}
V(ig)$n.label.size <- node.label.size
if (length(node.label.fontface) != nnode) {
if (!is.null(node.label.fontface)) {
tmp.node.label.fontface <- igraph::vertex_attr(ig,
node.label.fontface)
}
else {
tmp.node.label.fontface <- rep(0, nnode)
}
if (is.null(tmp.node.label.fontface)) {
node.label.fontface <- rep(node.label.fontface,
nnode)
}
else {
node.label.fontface <- tmp.node.label.fontface
}
}
V(ig)$n.label.fontface <- node.label.fontface
if (length(node.label.color) != nnode) {
if (!is.null(node.label.color)) {
tmp.node.label.color <- igraph::vertex_attr(ig,
node.label.color)
}
else {
tmp.node.label.color <- rep(0, nnode)
}
if (is.null(tmp.node.label.color)) {
node.label.color <- rep(node.label.color, nnode)
}
else {
node.label.color <- tmp.node.label.color
}
}
V(ig)$n.label.color <- node.label.color
if (length(node.color) != nnode) {
if (!is.null(node.color)) {
node.color <- igraph::vertex_attr(ig, node.color)
}
if (is.null(node.color)) {
node.color <- rep(0, nnode)
}
}
if (is.null(zlim)) {
zlim <- c(min(node.color), max(node.color))
}
node.color[node.color <= zlim[1]] <- zlim[1]
node.color[node.color >= zlim[2]] <- zlim[2]
V(ig)$n.color <- node.color
if (length(node.size) != nnode) {
if (!is.null(node.size)) {
tmp.node.size <- igraph::vertex_attr(ig, node.size)
}
else {
tmp.node.size <- rep(0, nnode)
}
if (is.null(tmp.node.size)) {
node.size <- rep(node.size, nnode)
}
else {
node.size <- tmp.node.size
}
}
if (is.null(slim)) {
slim <- c(min(node.size), max(node.size))
}
node.size[node.size <= slim[1]] <- slim[1]
node.size[node.size >= slim[2]] <- slim[2]
V(ig)$n.size <- node.size
if (length(node.shape) != nnode) {
if (!is.null(node.shape)) {
tmp.node.shape <- igraph::vertex_attr(ig, node.shape)
}
else {
tmp.node.shape <- rep(19, nnode)
}
if (is.null(tmp.node.shape)) {
node.shape <- rep(node.shape, nnode)
}
else {
node.shape <- tmp.node.shape
}
}
V(ig)$n.shape <- node.shape
nedge <- igraph::ecount(ig)
if (length(edge.color) != nedge) {
if (!is.null(edge.color)) {
tmp.edge.color <- igraph::edge_attr(ig, edge.color)
}
else {
tmp.edge.color <- rep("black", nedge)
}
if (is.null(tmp.edge.color)) {
edge.color <- rep(edge.color, nedge)
}
else {
edge.color <- tmp.edge.color
}
}
E(ig)$e.color <- edge.color
if (length(edge.color.alpha) != nedge) {
if (!is.null(edge.color.alpha)) {
tmp.edge.color.alpha <- igraph::edge_attr(ig,
edge.color.alpha)
}
else {
tmp.edge.color.alpha <- rep(0.5, nedge)
}
if (is.null(tmp.edge.color.alpha)) {
edge.color.alpha <- rep(edge.color.alpha, nedge)
}
else {
edge.color.alpha <- tmp.edge.color.alpha
}
}
E(ig)$e.color.alpha <- edge.color.alpha
if (length(edge.size) != nedge) {
if (!is.null(edge.size)) {
tmp.edge.size <- igraph::edge_attr(ig, edge.size)
}
else {
tmp.edge.size <- rep(0.5, nedge)
}
if (is.null(tmp.edge.size)) {
edge.size <- rep(edge.size, nedge)
}
else {
edge.size <- tmp.edge.size
}
}
E(ig)$e.size <- edge.size
gnet <- ggnetwork::ggnetwork(ig, layout = cbind(node.xcoord,
node.ycoord), arrow.gap = edge.arrow.gap, cell.jitter = 0.75)
data.frame(gnet, group = rep(names(ls_ig)[i], nrow(gnet)),
stringsAsFactors = F)
})
df <- do.call(rbind, ls_df)
i <- sapply(df, is.factor)
df[i] <- lapply(df[i], as.character)
df$group <- factor(df$group, levels = names(ls_ig))
df$n.color <- as.numeric(df$n.color)
df$n.size <- as.numeric(df$n.size)
df$n.shape <- factor(df$n.shape, levels = sort(unique(df$n.shape)))
n.color <- n.size <- n.shape <- n.label <- n.label.size <-
n.label.color <- NULL
x <- y <- xend <- yend <- NULL
gp <- ggplot(df, aes(x = x, y = y, xend = xend, yend = yend))
e.color <- subset(df, !is.na(na.y))$e.color
e.color.alpha <- subset(df, !is.na(na.y))$e.color.alpha
e.size <- subset(df, !is.na(na.y))$e.size
if (igraph::is_directed(ls_ig[[1]])) {
gp <- gp + ggnetwork::geom_edges(color = e.color, size = e.size,
alpha = e.color.alpha, curvature = edge.curve, arrow = arrow(length =
unit(edge.arrow,
"pt"), type = "closed"), show.legend = FALSE)
}
else {
gp <- gp + ggnetwork::geom_edges(color = e.color, size = e.size,
alpha = e.color.alpha, curvature = edge.curve, show.legend = FALSE)
}
if (length(unique(df$n.shape)) == 1) {
if (!is.numeric(node.shape)) {
node.shape <- 19
}
gp <- gp + ggnetwork::geom_nodes(aes(color = n.color,
size = n.size), shape = node.shape, alpha = node.color.alpha)
}
else {
gp <- gp + ggnetwork::geom_nodes(aes(color = n.color,
size = n.size, shape = n.shape), alpha = node.color.alpha)
gp <- gp + scale_shape(guide = guide_legend(node.shape.title,
title.position = "top", ncol = 1))
}
if (is.null(node.color.title)) {
node.color.title <- "Node color"
}
if (is.null(zlim)) {
zlim <- range(df$n.color[!is.na(df$n.color)])
}
if (zlim[1] != zlim[2]) {
gp <- gp + scale_colour_gradientn(colors =
xColormap(colormap)(ncolors),
limits = zlim, guide = guide_colorbar(title = node.color.title,
title.position = "top", barwidth = 0.5, nbin = 64,
draw.ulim = FALSE, draw.llim = FALSE), na.value = na.color)
}
else {
gp <- gp + scale_colour_gradientn(colors =
xColormap(colormap)(ncolors),
guide = guide_colorbar(title = node.color.title,
title.position = "top", barwidth = 0.5))
}
if (is.null(node.size.title)) {
node.size.title <- "Node size"
}
if (is.null(slim)) {
slim <- range(df$n.size)
}
if (slim[1] != slim[2]) {
gp <- gp + scale_size_continuous(limits = slim, range =
node.size.range,
guide = guide_legend(node.size.title, title.position = "top",
ncol = 1))
}
else {
gp <- gp + scale_size_continuous(limits = slim, range =
node.size.range,
guide = guide_legend(node.size.title, title.position = "top",
ncol = 1))
}
gp <- gp + ggnetwork::theme_blank()
gp <- gp + theme(text = element_text(family = "sans")) +
labs(title = title) + theme(plot.title = element_text(hjust = 0.5,
size = 10, face = "bold"), plot.margin = unit(rep(0,
4), rep("lines", 4)))
if ((zlim[1] != zlim[2]) & (slim[1] != slim[2])) {
gp <- gp + theme(legend.position = "right")
}
else {
if (slim[1] == slim[2]) {
gp <- gp + guides(size = "none")
}
if (zlim[1] == zlim[2]) {
gp <- gp + guides(color = "none")
}
if (length(unique(df$n.shape)) == 1) {
gp <- gp + guides(shape = "none")
}
}
gp <- gp + theme(legend.title = element_text(size = 8, face = "bold"),
legend.text = element_text(size = 6))
if (length(ls_ig) > 1) {
if (is.null(ncolumns)) {
ncolumns <- ceiling(sqrt(length(ls_ig)))
}
group <- NULL
gp <- gp + facet_wrap(~group, ncol = ncolumns)
gp <- gp + theme(strip.background = element_rect(fill = "transparent",
color = "transparent"), strip.text = element_text(size = 12,
face = "bold"), strip.placement = "inside", panel.spacing = unit(0,
"lines"))
}
if (sum(df$n.label == "") != nrow(df)) {
StatNodes <- ggplot2::ggproto("StatNodes", ggplot2::Stat,
compute_layer = function(data, scales, params) {
if (all(c("xend", "yend") \%in\% names(data))) {
unique(subset(data, select = c(-xend, -yend)))
}
else {
unique(data)
}
})
my_geom_nodetext_repel <- function(mapping = NULL, data = NULL,
parse = FALSE, ..., box.padding = unit(0.25, "lines"),
point.padding = unit(1e-06, "lines"), segment.size = 0.5,
arrow = NULL, force = 1, max.iter = 2000, nudge_x = 0,
nudge_y = 0, na.rm = FALSE, show.legend = F, inherit.aes = TRUE) {
ggplot2::layer(data = data, mapping = mapping, stat = StatNodes,
geom = ggrepel::GeomTextRepel, position = "identity",
show.legend = show.legend, inherit.aes = inherit.aes,
params = list(parse = parse, na.rm = na.rm, box.padding = box.padding,
point.padding = point.padding, segment.size = segment.size,
arrow = arrow, force = force, max.iter = max.iter,
nudge_x = nudge_x, nudge_y = nudge_y, ...))
}
n.label.size <- subset(df, is.na(na.y))$n.label.size
n.label.fontface <- subset(df, is.na(na.y))$n.label.fontface
n.label.color <- subset(df, is.na(na.y))$n.label.color
gp <- gp + my_geom_nodetext_repel(aes(label = n.label),
lineheight = label.wrap.lineheight, size = n.label.size,
color = n.label.color, fontface = n.label.fontface,
alpha = node.label.alpha, box.padding = unit(0.5,
"lines"), point.padding = unit(node.label.padding,
"lines"), segment.alpha = 0.2, segment.size = 0.2,
arrow = arrow(length = unit(node.label.arrow, "npc")),
force = node.label.force)
}
if (1) {
df <- gp$data
na.y <- NULL
df_sub <- subset(df, is.na(na.y))
ind <- match(colnames(df_sub), c("na.x", "na.y", "e.color",
"edge.color", "e.size"))
df_sub <- df_sub[, is.na(ind)]
df_sub <- df_sub[!duplicated(df_sub), ]
gp$data_nodes <- df_sub
df_sub <- subset(df, !is.na(na.y))
ind <- match(colnames(df_sub), c("x", "y", "xend", "yend",
"e.color", "edge.color", "e.size"))
df_sub <- df_sub[, !is.na(ind)]
df_sub <- df_sub[!duplicated(df_sub), ]
gp$data_edges <- df_sub
}
invisible(gp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
