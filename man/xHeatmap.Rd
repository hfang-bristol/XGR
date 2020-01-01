\name{xHeatmap}
\alias{xHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xHeatmap(data, reorder = c("none", "row", "col", "both"), colormap =
"spectral", ncolors = 64, zlim = NULL, barwidth = 0.3, barheight =
NULL, nbin = 64, legend.title = "", x.rotate = 90, x.text.size = 6,
x.text.hjust = 0, y.text.size = 6, legend.text.size = 4,
legend.title.size = 6, shape = 19, size = 2, plot.margin = unit(c(5.5,
5.5, 5.5, 5.5), "pt"), font.family = "sans", na.color = "transparent",
data.label = NULL, label.size = 1, label.color = "black", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
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
function (data, reorder = c("none", "row", "col", "both"), colormap =
"spectral",
ncolors = 64, zlim = NULL, barwidth = 0.3, barheight = NULL,
nbin = 64, legend.title = "", x.rotate = 90, x.text.size = 6,
x.text.hjust = 0, y.text.size = 6, legend.text.size = 4,
legend.title.size = 6, shape = 19, size = 2, plot.margin = unit(c(5.5,
5.5, 5.5, 5.5), "pt"), font.family = "sans", na.color = "transparent",
data.label = NULL, label.size = 1, label.color = "black",
...)
{
reorder <- match.arg(reorder)
if (is.null(rownames(data))) {
rownames(data) <- paste("R", 1:nrow(data), sep = "")
}
if (is.null(colnames(data))) {
colnames(data) <- paste("C", 1:ncol(data), sep = "")
}
mat_val <- data
if (!(all(!duplicated(rownames(mat_val))) &
all(!duplicated(colnames(mat_val))))) {
rownames(mat_val) <- make.names(rownames(mat_val), unique = TRUE)
colnames(mat_val) <- make.names(colnames(mat_val), unique = TRUE)
}
flag_factor <- FALSE
if (class(unlist(mat_val)) == "factor") {
lvl <- levels(unlist(mat_val))
indx <- sapply(mat_val, is.factor)
mat_val[indx] <- sapply(mat_val[, indx], function(x) as.numeric(x))
flag_factor <- F
}
if (nrow(mat_val) > 1 & nrow(mat_val) > 1) {
ind_row <- 1:nrow(mat_val)
if (ncol(mat_val) > 1 & (reorder == "row" | reorder ==
"both")) {
mat <- mat_val
colnames(mat) <- 1:ncol(mat)
rownames(mat) <- 1:nrow(mat)
mat[is.na(mat)] <- 0
set.seed(825)
tree_bs <- visTreeBootstrap(mat, visTree = FALSE,
...)
ind_row <- match(tree_bs$tip.label, rownames(mat))
}
ind_col <- 1:ncol(mat_val)
if (nrow(mat_val) > 1 & (reorder == "col" | reorder ==
"both")) {
mat <- mat_val
colnames(mat) <- 1:ncol(mat)
rownames(mat) <- 1:nrow(mat)
mat[is.na(mat)] <- 0
set.seed(825)
tree_bs <- visTreeBootstrap(t(mat), visTree = FALSE,
...)
ind_col <- match(tree_bs$tip.label, colnames(mat))
}
mat_tmp <- as.matrix(mat_val[ind_row, ind_col], ncol = length(ind_col))
rownames(mat_tmp) <- rownames(mat_val)[ind_row]
colnames(mat_tmp) <- colnames(mat_val)[ind_col]
mat_val <- mat_tmp
}
if (class(mat_val) == "matrix") {
mat_val <- as.data.frame(mat_val)
}
if (class(mat_val) == "data.frame") {
if (is.null(zlim)) {
zlim <- c(floor(min(mat_val, na.rm = T) * 10)/10,
ceiling(max(mat_val, na.rm = T) * 10)/10)
if (zlim[1] == zlim[2]) {
zlim[1] <- zlim[2]/2
}
}
mat_val[mat_val <= zlim[1]] <- zlim[1]
mat_val[mat_val >= zlim[2]] <- zlim[2]
gene <- sample <- val <- NULL
df <- mat_val \%>\% dplyr::mutate(gene = rownames(mat_val)) \%>\%
tidyr::gather(sample, val, -gene)
df$gene <- factor(df$gene, levels = rev(rownames(mat_val)))
df$sample <- factor(df$sample, levels = colnames(mat_val))
df <- df \%>\% dplyr::mutate(uid = paste(gene, sample,
sep = ":"))
df$y <- as.numeric(df$gene)
df$x <- as.numeric(df$sample)
if (flag_factor) {
df$val <- factor(df$val, levels = unique(df$val))
}
gp <- ggplot(df, aes(x = sample, y = gene, color = val))
gp <- gp + geom_point(size = size, shape = shape)
if (!flag_factor) {
gp <- gp + scale_colour_gradientn(colors =
xColormap(colormap)(ncolors),
limits = zlim, guide = guide_colorbar(title = legend.title,
title.position = "top", barwidth = barwidth,
barheight = barheight, nbin = nbin, draw.ulim = FALSE,
draw.llim = FALSE), na.value = na.color)
}
else {
gp <- gp + scale_colour_manual(legend.title, values =
xColormap(colormap)(length(lvl)),
labels = lvl)
if (is.null(barheight)) {
barheight <- barwidth
}
gp <- gp + theme(legend.key.width = unit(barwidth,
"pt"), legend.key.height = unit(barheight, "pt")) +
guides(col = guide_legend(ncol = 1))
}
gp <- gp + theme_bw() + theme(legend.position = "right",
axis.title.x = element_blank(), axis.title.y = element_blank(),
axis.text.x = element_text(face = "bold", color = "black",
size = x.text.size, angle = x.rotate, hjust = x.text.hjust),
axis.text.y = element_text(face = "bold", color = "black",
size = y.text.size, angle = 0), panel.background = element_rect(fill =
"transparent")) +
theme(panel.grid.major = element_blank(), panel.grid.minor =
element_blank()) +
theme(plot.margin = plot.margin) + theme(legend.title =
element_text(face = "bold",
color = "black", size = legend.title.size), legend.text =
element_text(face = "bold",
color = "black", size = legend.text.size), legend.title.align = 0.5) +
theme(legend.background = element_rect(fill = "transparent"))
gp <- gp + theme(axis.ticks = element_line(size = 0.25),
axis.ticks.length = unit(0.1, "cm"))
gp <- gp + theme(text = element_text(family = font.family))
gp_main <- gp + scale_x_discrete(position = "top")
if (!is.null(data.label)) {
mat_label <- as.data.frame(data.label[ind_row, ind_col],
stringsAsFactors = F)
gene <- sample <- val <- NULL
df_label <- suppressWarnings(mat_label \%>\% dplyr::mutate(gene =
rownames(mat_label)) \%>\%
tidyr::gather(sample, val, -gene))
df_label$gene <- factor(df_label$gene, levels =
rev(rownames(mat_label)))
df_label$sample <- factor(df_label$sample, levels =
colnames(mat_label))
df_label <- df_label \%>\% dplyr::mutate(uid = paste(gene,
sample, sep = ":"))
df_label$y <- as.numeric(df_label$gene)
df_label$x <- as.numeric(df_label$sample)
gp_main <- gp_main + geom_text(data = df_label, aes(x = sample,
y = gene, label = val), hjust = 0.5, vjust = 0.5,
size = label.size, color = label.color)
}
invisible(gp_main)
}
else {
return(NULL)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
