\name{xVisNet}
\alias{xVisNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xVisNet(g, pattern = NULL, colormap = c("yr", "jet", "gbr", "wyr",
"br", "bwr", "rainbow", "wb"), ncolors = 40, zlim = NULL, colorbar =
TRUE, newpage = TRUE, signature = TRUE, glayout = layout_with_kk,
vertex.frame.color = NA, vertex.size = NULL, vertex.color = NULL,
vertex.shape = NULL, vertex.label = NULL, vertex.label.cex = NULL,
vertex.label.dist = 0.3, vertex.label.color = "blue",
vertex.label.family = "sans", edge.arrow.size = 0.3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{g}{
%% ~~Describe \code{g} here~~
}
\item{pattern}{
%% ~~Describe \code{pattern} here~~
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
\item{signature}{
%% ~~Describe \code{signature} here~~
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
\item{vertex.label.family}{
%% ~~Describe \code{vertex.label.family} here~~
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
function (g, pattern = NULL, colormap = c("yr", "jet", "gbr",
"wyr", "br", "bwr", "rainbow", "wb"), ncolors = 40, zlim = NULL,
colorbar = TRUE, newpage = TRUE, signature = TRUE, glayout =
layout_with_kk,
vertex.frame.color = NA, vertex.size = NULL, vertex.color = NULL,
vertex.shape = NULL, vertex.label = NULL, vertex.label.cex = NULL,
vertex.label.dist = 0.3, vertex.label.color = "blue",
vertex.label.family = "sans",
edge.arrow.size = 0.3, ...)
{
if (class(g) != "igraph") {
stop("The function must apply to either 'igraph' or 'graphNEL'
object.\n")
}
else {
ig <- g
}
if (igraph::vcount(ig) == 0) {
return(NULL)
}
if (is_directed(ig)) {
}
if (!newpage) {
if (dev.cur() > 1) {
}
}
if (!is.null(pattern)) {
if (class(pattern) == "character") {
tmp <- as.numeric(pattern)
if (!is.null(names(pattern))) {
names(tmp) <- names(pattern)
}
pattern <- tmp
}
pattern[is.infinite(pattern)] <- max(pattern[!is.infinite(pattern)])
if (is.null(zlim)) {
if (max(pattern) < 1 & max(pattern) > 0) {
vmax <- ceiling(stats::quantile(pattern, 0.75) *
10)/10
}
else {
vmax <- ceiling(stats::quantile(pattern, 0.75))
}
vmin <- floor(min(pattern))
zlim <- c(vmin, vmax)
}
}
suppressWarnings(par_old <- graphics::par(no.readonly = TRUE))
dnet::visNet(g = ig, pattern = pattern, colormap = colormap,
ncolors = ncolors, zlim = zlim, colorbar = colorbar,
newpage = newpage, glayout = glayout, vertex.frame.color =
vertex.frame.color,
vertex.size = vertex.size, vertex.color = vertex.color,
vertex.shape = vertex.shape, vertex.label = vertex.label,
vertex.label.cex = vertex.label.cex, vertex.label.dist =
vertex.label.dist,
vertex.label.color = vertex.label.color, vertex.label.family =
vertex.label.family,
edge.arrow.size = edge.arrow.size, ...)
suppressWarnings(graphics::par(par_old))
if (signature) {
caption <- paste("Created by xVisNet from XGR version",
utils::packageVersion("XGR"))
graphics::mtext(caption, side = 1, line = 2, adj = 1,
cex = 0.66, font = 3, family = vertex.label.family)
}
invisible()
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
