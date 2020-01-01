\name{xHeatmapAdv}
\alias{xHeatmapAdv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xHeatmapAdv(data.main, data.meta, reorder = c("none", "row", "col",
"both"), colormap = "spectral", ncolors = 64, zlim = NULL, barwidth =
0.3, barheight = 4, nbin = 64, legend.title = "Main", x.rotate = 60,
x.text.size = 6, x.text.hjust = 0.5, y.text.size = 6, legend.text.size
= 5, legend.title.size = 6, shape = 19, size = 2, plot.margin =
unit(c(5.5, 5.5, 5.5, 5.5), "pt"), font.family = "sans", na.color =
"grey80", data.label = NULL, label.size = 1, label.color = "black",
meta.colormap = "spectral", meta.x.rotate = 75, meta.shape.continuous =
15, meta.shape.discrete = 95, meta.size = 2, meta.location = c("right",
"left"), meta.width = 0.5, gap.width = 0.5, legend.width = NULL,
legend.direction = c("vertical", "horizontal"), legend.nrow = NULL,
verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data.main}{
%% ~~Describe \code{data.main} here~~
}
\item{data.meta}{
%% ~~Describe \code{data.meta} here~~
}
\item{reorder}{
%% ~~Describe \code{reorder} here~~
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
\item{barwidth}{
%% ~~Describe \code{barwidth} here~~
}
\item{barheight}{
%% ~~Describe \code{barheight} here~~
}
\item{nbin}{
%% ~~Describe \code{nbin} here~~
}
\item{legend.title}{
%% ~~Describe \code{legend.title} here~~
}
\item{x.rotate}{
%% ~~Describe \code{x.rotate} here~~
}
\item{x.text.size}{
%% ~~Describe \code{x.text.size} here~~
}
\item{x.text.hjust}{
%% ~~Describe \code{x.text.hjust} here~~
}
\item{y.text.size}{
%% ~~Describe \code{y.text.size} here~~
}
\item{legend.text.size}{
%% ~~Describe \code{legend.text.size} here~~
}
\item{legend.title.size}{
%% ~~Describe \code{legend.title.size} here~~
}
\item{shape}{
%% ~~Describe \code{shape} here~~
}
\item{size}{
%% ~~Describe \code{size} here~~
}
\item{plot.margin}{
%% ~~Describe \code{plot.margin} here~~
}
\item{font.family}{
%% ~~Describe \code{font.family} here~~
}
\item{na.color}{
%% ~~Describe \code{na.color} here~~
}
\item{data.label}{
%% ~~Describe \code{data.label} here~~
}
\item{label.size}{
%% ~~Describe \code{label.size} here~~
}
\item{label.color}{
%% ~~Describe \code{label.color} here~~
}
\item{meta.colormap}{
%% ~~Describe \code{meta.colormap} here~~
}
\item{meta.x.rotate}{
%% ~~Describe \code{meta.x.rotate} here~~
}
\item{meta.shape.continuous}{
%% ~~Describe \code{meta.shape.continuous} here~~
}
\item{meta.shape.discrete}{
%% ~~Describe \code{meta.shape.discrete} here~~
}
\item{meta.size}{
%% ~~Describe \code{meta.size} here~~
}
\item{meta.location}{
%% ~~Describe \code{meta.location} here~~
}
\item{meta.width}{
%% ~~Describe \code{meta.width} here~~
}
\item{gap.width}{
%% ~~Describe \code{gap.width} here~~
}
\item{legend.width}{
%% ~~Describe \code{legend.width} here~~
}
\item{legend.direction}{
%% ~~Describe \code{legend.direction} here~~
}
\item{legend.nrow}{
%% ~~Describe \code{legend.nrow} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
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
function (data.main, data.meta, reorder = c("none", "row", "col",
"both"), colormap = "spectral", ncolors = 64, zlim = NULL,
barwidth = 0.3, barheight = 4, nbin = 64, legend.title = "Main",
x.rotate = 60, x.text.size = 6, x.text.hjust = 0.5, y.text.size = 6,
legend.text.size = 5, legend.title.size = 6, shape = 19,
size = 2, plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
font.family = "sans", na.color = "grey80", data.label = NULL,
label.size = 1, label.color = "black", meta.colormap = "spectral",
meta.x.rotate = 75, meta.shape.continuous = 15, meta.shape.discrete =
95,
meta.size = 2, meta.location = c("right", "left"), meta.width = 0.5,
gap.width = 0.5, legend.width = NULL, legend.direction = c("vertical",
"horizontal"), legend.nrow = NULL, verbose = TRUE, ...)
{
reorder <- match.arg(reorder)
meta.location <- match.arg(meta.location)
legend.direction <- match.arg(legend.direction)
if (is.null(legend.nrow)) {
if (legend.direction == "vertical") {
legend.nrow <- 3
}
else if (legend.direction == "horizontal") {
legend.nrow <- 6
}
}
gp_main <- xHeatmap(data.main, reorder = reorder, colormap = colormap,
ncolors = ncolors, zlim = zlim, barwidth = barwidth,
barheight = barheight, nbin = nbin, legend.title = legend.title,
x.rotate = x.rotate, x.text.size = x.text.size, x.text.hjust =
x.text.hjust,
y.text.size = y.text.size, legend.text.size = legend.text.size,
legend.title.size = legend.title.size, shape = shape,
size = size, plot.margin = plot.margin, font.family = font.family,
na.color = na.color, data.label = data.label, label.size = label.size,
label.color = label.color, ...)
if (legend.direction == "horizontal") {
gp_main <- gp_main + theme(legend.position = "right",
legend.direction = "horizontal") + guides(color = guide_colorbar(title
= legend.title,
title.position = "top", barwidth = barheight, barheight = barwidth))
}
gp_main_void <- gp_main + theme(legend.position = "none",
axis.title.x = element_blank(), axis.title.y = element_blank())
if (is.null(rownames(data.meta))) {
rownames(data.meta) <- paste("R", 1:nrow(data.meta),
sep = " ")
}
if (is.null(colnames(data.meta))) {
colnames(data.meta) <- paste("C", 1:ncol(data.meta),
sep = " ")
}
if (class(data.meta) == "matrix") {
data.meta <- as.data.frame(data.meta)
}
main_y_texts <- rev(levels(gp_main$data$gene))
main_y_ind <- match(main_y_texts, rownames(data.meta))
meta <- data.meta \%>\% dplyr::arrange(main_y_ind)
rownames(meta) <- rownames(data.meta)[main_y_ind]
ls_gp_meta <- lapply(1:ncol(meta), function(i) {
data <- meta \%>\% dplyr::select(i)
meta.legend.title <- colnames(meta)[i]
colormap <- meta.colormap
meta.shape <- meta.shape.continuous
meta.plot.margin <- unit(c(plot.margin[[1]], 0, plot.margin[[3]],
0), "pt")
if (is.factor(unlist(data))) {
meta.shape <- meta.shape.discrete
if (meta.shape > 25) {
meta.size <- meta.size * 1.5
}
}
gp <- xHeatmap(data, reorder = "none", colormap = meta.colormap,
ncolors = ncolors, zlim = NULL, barwidth = barwidth,
barheight = barheight, nbin = nbin, legend.title = meta.legend.title,
x.rotate = meta.x.rotate, x.text.size = x.text.size,
x.text.hjust = x.text.hjust, legend.text.size = legend.text.size,
legend.title.size = legend.title.size, shape = meta.shape,
size = meta.size, plot.margin = meta.plot.margin,
font.family = font.family, na.color = na.color)
if (legend.direction == "horizontal") {
gp <- gp + theme(legend.position = "right", legend.direction =
"horizontal") +
guides(color = guide_colorbar(title = meta.legend.title,
title.position = "top", barwidth = barheight,
barheight = barwidth))
}
return(gp)
})
ls_gp_meta_void <- lapply(ls_gp_meta, function(gp) {
gp_void <- gp + theme(legend.position = "none", axis.title.x =
element_blank(),
axis.title.y = element_blank(), axis.text.y = element_blank(),
axis.ticks.y = element_blank(), panel.border = element_blank())
})
extract_gp_legend <- function(gp) {
gt <- ggplot_gtable(ggplot_build(gp))
ind <- which(sapply(gt$grobs, function(x) x$name) ==
"guide-box")
legend <- gt$grobs[[ind]]
invisible(legend)
}
ls_gp <- c(list(gp_main), ls_gp_meta)
ls_gp_legend <- lapply(ls_gp, function(gp) {
extract_gp_legend(gp)
})
gt_legend <- gridExtra::arrangeGrob(grobs = ls_gp_legend,
nrow = legend.nrow, as.table = TRUE, top = "", padding = unit(c(1,
0, 0, 0), "line"))
ls_gp_void <- c(list(gp_main_void), ls_gp_meta_void)
ls_gt_void <- lapply(ls_gp_void, function(gp) {
gt <- ggplot_gtable(ggplot_build(gp))
})
ls_heights <- lapply(ls_gt_void, function(gt) {
gt$heights[4:5]
})
maxHeight <- do.call(grid::unit.pmax, ls_heights)
ls_gt_void_aligned <- lapply(ls_gt_void, function(gt) {
gt$heights[4:5] <- maxHeight
invisible(gt)
})
n1 <- ncol(data.main)
n2 <- ncol(data.meta)
gt_gap <- grid::rectGrob(gp = grid::gpar(col = "transparent"))
if (is.null(legend.width)) {
legend.width <- 3 * meta.width * ceiling((n2 + 1)/legend.nrow)
}
if (meta.location == "right") {
grobs <- c(ls_gt_void_aligned[1], list(gt_gap), ls_gt_void_aligned[-1],
list(gt_gap), list(gt_legend))
widths <- c(n1, gap.width, rep(meta.width, n2), gap.width,
legend.width)
}
else if (meta.location == "left") {
grobs <- c(list(gt_gap), ls_gt_void_aligned[-1], ls_gt_void_aligned[1],
list(gt_gap), list(gt_legend))
widths <- c(gap.width, rep(meta.width, n2), n1, gap.width,
legend.width)
}
gt_heatmap <- gridExtra::arrangeGrob(grobs = grobs, nrow = 1,
widths = widths)
if (verbose) {
message(sprintf("data.main: \%d rows X \%d columns", nrow(data.main),
ncol(data.main)), appendLF = TRUE)
n <- sum(sapply(data.meta, is.factor))
m <- ncol(data.meta) - n
message(sprintf("data.meta: \%d rows X \%d columns (\%d continuous +
\%d discrete)",
nrow(data.meta), ncol(data.meta), m, n), appendLF = TRUE)
message(sprintf("widths: \%.2f (meta), \%.2f (blank), \%.2f (legend)",
meta.width, gap.width, legend.width), appendLF = TRUE)
}
invisible(gt_heatmap)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
