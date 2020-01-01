\name{xEnrichDAGplotAdv}
\alias{xEnrichDAGplotAdv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichDAGplotAdv(ggplot, displayBy = c("nSig", "none"), path.mode =
c("all_paths", "shortest_paths", "all_shortest_paths"), height = 7,
width = 7, margin = rep(0.1, 4), colormap = c("white-lightcyan-cyan",
"yr", "bwr", "jet", "gbr", "wyr", "br", "rainbow", "wb",
"lightyellow-orange"), ncolors = 40, zlim = NULL, colorbar = T,
colorbar.fraction = 0.1, newpage = T, layout.orientation =
c("left_right", "top_bottom", "bottom_top", "right_left"), node.info =
c("term_name", "term_id", "none"), wrap.width = NULL, graph.node.attrs
= NULL, graph.edge.attrs = NULL, node.attrs = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ggplot}{
%% ~~Describe \code{ggplot} here~~
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
function (ggplot, displayBy = c("nSig", "none"), path.mode =
c("all_paths",
"shortest_paths", "all_shortest_paths"), height = 7, width = 7,
margin = rep(0.1, 4), colormap = c("white-lightcyan-cyan",
"yr", "bwr", "jet", "gbr", "wyr", "br", "rainbow", "wb",
"lightyellow-orange"), ncolors = 40, zlim = NULL, colorbar = T,
colorbar.fraction = 0.1, newpage = T, layout.orientation =
c("left_right",
"top_bottom", "bottom_top", "right_left"), node.info = c("term_name",
"term_id", "none"), wrap.width = NULL, graph.node.attrs = NULL,
graph.edge.attrs = NULL, node.attrs = NULL)
{
displayBy <- match.arg(displayBy)
path.mode <- match.arg(path.mode)
layout.orientation <- match.arg(layout.orientation)
node.info <- match.arg(node.info)
if (any(class(ggplot) \%in\% c("gg", "ggplot"))) {
bp <- ggplot
if (!is.null(bp$g)) {
bp_ig <- bp$g
bp_data <- bp$data
nodes_query <- unique(bp_data$id)
subg <- dDAGinduce(g = bp_ig, nodes_query = nodes_query,
path.mode = path.mode)
}
else {
stop("No 'igraph' object can not be found in the input 'ggplot'
object.\n")
}
}
else {
stop("Cannot find either an 'eTerm' object or an 'ggplot' object.\n")
}
code <- bp_data$code
names(code) <- bp_data$id
ind <- match(V(subg)$name, names(code))
code_in_subg <- code[ind]
names(code_in_subg) <- V(subg)$name
V(subg)$code <- code_in_subg
nSig <- bp_data$nSig
names(nSig) <- bp_data$id
ind <- match(V(subg)$name, names(nSig))
nSig_in_subg <- nSig[ind]
names(nSig_in_subg) <- V(subg)$name
V(subg)$nSig <- nSig_in_subg
for (i in 1:length(V(subg)$code)) {
if (!is.na(V(subg)$code[i])) {
V(subg)$name[i] <- paste(V(subg)$code[i], V(subg)$name[i],
sep = "\\\n")
V(subg)$term_name[i] <- paste(V(subg)$code[i], V(subg)$term_name[i],
sep = "\\\n")
}
}
if (displayBy == "nSig") {
data <- V(subg)$nSig
names(data) <- V(subg)$name
data <- data[!is.na(data)]
if (is.null(zlim)) {
zlim <- c(0, max(data))
}
}
else {
data <- NULL
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
