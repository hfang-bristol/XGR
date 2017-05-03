\name{xEnrichNetplot}
\alias{xEnrichNetplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichNetplot(eTerm, top_num = 10, displayBy = c("fc", "adjp", "fdr",
"zscore", "pvalue"), path.mode = c("all_paths", "shortest_paths",
"all_shortest_paths"), node.info = c("none", "term_id", "term_name",
"both", "full_term_name"), wrap.width = 15, colormap = c("yr", "jet",
"gbr", "wyr", "br", "bwr", "rainbow", "wb"), ncolors = 40, zlim = NULL,
colorbar = T, newpage = T, glayout = layout_as_tree, vertex.frame.color
= NA, vertex.size = NULL, vertex.color = NULL, vertex.shape = NULL,
vertex.label = NULL, vertex.label.cex = NULL, vertex.label.dist = 0.3,
vertex.label.color = "blue", edge.arrow.size = 0.3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{eTerm}{
%% ~~Describe \code{eTerm} here~~
}
\item{top_num}{
%% ~~Describe \code{top_num} here~~
}
\item{displayBy}{
%% ~~Describe \code{displayBy} here~~
}
\item{path.mode}{
%% ~~Describe \code{path.mode} here~~
}
\item{node.info}{
%% ~~Describe \code{node.info} here~~
}
\item{wrap.width}{
%% ~~Describe \code{wrap.width} here~~
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
\item{newpage}{
%% ~~Describe \code{newpage} here~~
}
\item{glayout}{
%% ~~Describe \code{glayout} here~~
}
\item{vertex.frame.color}{
%% ~~Describe \code{vertex.frame.color} here~~
}
\item{vertex.size}{
%% ~~Describe \code{vertex.size} here~~
}
\item{vertex.color}{
%% ~~Describe \code{vertex.color} here~~
}
\item{vertex.shape}{
%% ~~Describe \code{vertex.shape} here~~
}
\item{vertex.label}{
%% ~~Describe \code{vertex.label} here~~
}
\item{vertex.label.cex}{
%% ~~Describe \code{vertex.label.cex} here~~
}
\item{vertex.label.dist}{
%% ~~Describe \code{vertex.label.dist} here~~
}
\item{vertex.label.color}{
%% ~~Describe \code{vertex.label.color} here~~
}
\item{edge.arrow.size}{
%% ~~Describe \code{edge.arrow.size} here~~
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
function (eTerm, top_num = 10, displayBy = c("fc", "adjp", "fdr",
"zscore", "pvalue"), path.mode = c("all_paths", "shortest_paths",
"all_shortest_paths"), node.info = c("none", "term_id", "term_name",
"both", "full_term_name"), wrap.width = 15, colormap = c("yr",
"jet", "gbr", "wyr", "br", "bwr", "rainbow", "wb"), ncolors = 40,
zlim = NULL, colorbar = T, newpage = T, glayout = layout_as_tree,
vertex.frame.color = NA, vertex.size = NULL, vertex.color = NULL,
vertex.shape = NULL, vertex.label = NULL, vertex.label.cex = NULL,
vertex.label.dist = 0.3, vertex.label.color = "blue", edge.arrow.size =
0.3,
...)
{
displayBy <- match.arg(displayBy)
path.mode <- match.arg(path.mode)
node.info <- match.arg(node.info)
if (is.logical(eTerm)) {
stop("There is no enrichment in the 'eTerm' object.\n")
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
g <- eTerm$g
nodes_query <- rownames(df)
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
ind <- match(V(subg)$name, names(data))
data <- data[ind]
names(data) <- V(subg)$name
}
else if (any(class(eTerm) \%in\% c("gg", "ggplot"))) {
bp <- eTerm
if (!is.null(bp$g)) {
bp_ig <- bp$g
bp_data <- bp$data
nodes_query <- unique(bp_data$id)
subg <- dDAGinduce(g = bp_ig, nodes_query = nodes_query,
path.mode = path.mode)
data <- bp_data$nSig
if (is.null(zlim)) {
zlim <- c(0, max(data))
}
}
else {
stop("No 'igraph' object can not be found in the input 'ggplot'
object.\n")
}
}
else {
stop("Cannot find either an 'eTerm' object or an 'ggplot' object.\n")
}
if (is.null(vertex.label)) {
getTermInfo <- function(g, vids, numChar = 15, mulLines = F) {
fullNames <- V(g)[vids]$term_name
names(fullNames) <- V(g)[vids]$name
if (mulLines == F) {
shortNames <- paste(substr(fullNames, 1, numChar),
ifelse(nchar(fullNames) > numChar, "...", ""),
sep = "")
}
else {
shortNames <- sapply(fullNames, function(x) {
return(paste(strwrap(x, numChar), sep = "",
collapse = "\n"))
})
}
names(shortNames) <- names(fullNames)
return(shortNames)
}
termNames <- getTermInfo(subg, vids = V(subg)$term_id,
numChar = 15, mulLines = F)
vertex.label <- switch(node.info, none = NULL, term_id =
V(subg)$term_id,
term_name = V(subg)$term_name, both = paste(V(subg)$term_id,
termNames, sep = "\n"), full_term_name = getTermInfo(subg,
vids = V(subg)$term_id, numChar = wrap.width,
mulLines = T))
}
pattern <- data
if (!is.null(pattern)) {
flag <- 0
if (!is.null(names(pattern))) {
pattern <- pattern[V(subg)$name]
}
if (length(pattern) == vcount(subg)) {
flag <- 1
}
if (flag == 1) {
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
palette.name <- supraHex::visColormap(colormap = colormap)
colors <- palette.name(ncolors)
scale <- length(colors)/(max(zlim) - min(zlim))
sapply(1:length(vec), function(x) {
if (is.na(vec[x])) {
"transparent"
}
else {
ind <- floor(1 + (vec[x] - min(zlim)) * scale)
colors[max(1, min(ncolors, ind))]
}
})
}
vertex.color <- vec2color(pattern, colormap = colormap,
ncolors = ncolors, zlim = zlim)
vertex.frame.color <- vec2color(pattern, colormap = colormap,
ncolors = ncolors, zlim = zlim)
vertex.frame.color[vertex.frame.color == "transparent"] <- "grey"
}
else {
warning("The input 'pattern' is ignored. Please check the help for
enabling your input")
pattern <- NULL
if (is.null(vertex.color)) {
vertex.color <- "SkyBlue2"
}
}
}
else {
if (is.null(vertex.color)) {
vertex.color <- "SkyBlue2"
}
}
V(subg)$enrichment <- pattern
par_old <- graphics::par()
dnet::visNet(g = subg, pattern = pattern, colormap = colormap,
ncolors = ncolors, zlim = zlim, colorbar = colorbar,
newpage = newpage, glayout = glayout, vertex.frame.color =
vertex.frame.color,
vertex.size = vertex.size, vertex.color = vertex.color,
vertex.shape = vertex.shape, vertex.label = vertex.label,
vertex.label.cex = vertex.label.cex, vertex.label.dist =
vertex.label.dist,
vertex.label.color = vertex.label.color, edge.arrow.size =
edge.arrow.size,
...)
suppressWarnings(graphics::par(par_old))
invisible(subg)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
