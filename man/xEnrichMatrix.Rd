\name{xEnrichMatrix}
\alias{xEnrichMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichMatrix(list_eTerm, method = c("ggplot2", "circle", "square",
"color", "pie"), displayBy = c("zscore", "fc", "adjp", "pvalue"),
FDR.cutoff = 0.05, wrap.width = NULL, sharings = NULL, reorder =
c("row", "none", "col", "both"), colormap = "jet", ncolors = 20, zlim =
NULL, slim = NULL, legend.direction = c("horizontal", "vertical"),
title = NULL, flip = FALSE, y.rotate = 45, shape = 19, font.family =
"sans", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{list_eTerm}{
%% ~~Describe \code{list_eTerm} here~~
}
\item{method}{
%% ~~Describe \code{method} here~~
}
\item{displayBy}{
%% ~~Describe \code{displayBy} here~~
}
\item{FDR.cutoff}{
%% ~~Describe \code{FDR.cutoff} here~~
}
\item{wrap.width}{
%% ~~Describe \code{wrap.width} here~~
}
\item{sharings}{
%% ~~Describe \code{sharings} here~~
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
\item{slim}{
%% ~~Describe \code{slim} here~~
}
\item{legend.direction}{
%% ~~Describe \code{legend.direction} here~~
}
\item{title}{
%% ~~Describe \code{title} here~~
}
\item{flip}{
%% ~~Describe \code{flip} here~~
}
\item{y.rotate}{
%% ~~Describe \code{y.rotate} here~~
}
\item{shape}{
%% ~~Describe \code{shape} here~~
}
\item{font.family}{
%% ~~Describe \code{font.family} here~~
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
function (list_eTerm, method = c("ggplot2", "circle", "square",
"color", "pie"), displayBy = c("zscore", "fc", "adjp", "pvalue"),
FDR.cutoff = 0.05, wrap.width = NULL, sharings = NULL, reorder =
c("row",
"none", "col", "both"), colormap = "jet", ncolors = 20,
zlim = NULL, slim = NULL, legend.direction = c("horizontal",
"vertical"), title = NULL, flip = FALSE, y.rotate = 45,
shape = 19, font.family = "sans", ...)
{
method <- match.arg(method)
displayBy <- match.arg(displayBy)
reorder <- match.arg(reorder)
legend.direction <- match.arg(legend.direction)
if (class(list_eTerm) == "list") {
list_eTerm <- base::Filter(base::Negate(is.null), list_eTerm)
list_names <- names(list_eTerm)
if (is.null(list_names)) {
list_names <- paste("Enrichment", 1:length(list_eTerm),
sep = " ")
}
res_ls <- lapply(1:length(list_eTerm), function(i) {
df <- xEnrichViewer(list_eTerm[[i]], top_num = "all",
sortBy = "none")
if (is.null(df)) {
return(NULL)
}
else {
cbind(group = rep(list_names[i], nrow(df)), id = rownames(df),
df, stringsAsFactors = F)
}
})
df_all <- do.call(rbind, res_ls)
rownames(df_all) <- NULL
ind <- which(df_all$adjp < FDR.cutoff)
d <- df_all[ind, c("id", "name", "fc", "adjp", "zscore",
"pvalue", "group")]
}
else if (class(list_eTerm) == "data.frame") {
ind <- colnames(list_eTerm) \%in\% c("group", "name", "adjp",
displayBy)
if (sum(ind) == 4) {
df_all <- list_eTerm
ind <- which(df_all$adjp < FDR.cutoff)
d <- df_all[ind, ]
if (is.factor(d$group)) {
list_names <- levels(d$group)
}
else {
list_names <- unique(d$group)
}
}
else {
return(NULL)
}
}
d$group <- factor(d$group, levels = rev(list_names))
id_ls <- split(x = d$group, f = d$name)
ind <- match(d$name, names(id_ls))
id_full_ls <- id_ls[ind]
nSig <- unlist(lapply(id_full_ls, length))
d$nSig <- nSig
code <- lapply(id_full_ls, function(x) {
res <- rep(0, length(levels(x)))
ind <- match(x, levels(x))
res[ind] <- 1
paste(res, collapse = "-")
})
d$code <- unlist(code)
if (!is.null(wrap.width)) {
width <- as.integer(wrap.width)
res_list <- lapply(d$name, function(x) {
x <- gsub("_", " ", x)
y <- strwrap(x, width = width)
if (length(y) > 1) {
paste0(y[1], "...")
}
else {
y
}
})
d$name <- unlist(res_list)
}
if (!is.null(sharings)) {
sharings <- as.numeric(sharings)
ind <- match(sharings, unique(d$nSig))
found <- sharings[!is.na(ind)]
if (length(found) > 0) {
flag <- match(d$nSig, found)
d <- d[!is.na(flag), ]
nSig <- nSig[!is.na(flag)]
}
}
if (displayBy == "fc") {
d <- d[with(d, order(nSig, group, fc, -adjp)), ]
if (!is.factor(d$name)) {
d$name <- factor(d$name, levels = unique(d$name))
}
d$val <- d$fc
mat_val <- as.matrix(xSparseMatrix(d[, c("name", "group",
"val")], rows = levels(d$name), columns = levels(d$group)))
mat_val[mat_val == 0] <- 1
title_size <- "FC"
}
else if (displayBy == "adjp") {
d$adjp[d$adjp == 0] <- min(d$adjp[d$adjp != 0])
d <- d[with(d, order(nSig, group, -adjp)), ]
if (!is.factor(d$name)) {
d$name <- factor(d$name, levels = unique(d$name))
}
d$val <- -1 * log10(d$adjp)
mat_val <- as.matrix(xSparseMatrix(d[, c("name", "group",
"val")], rows = levels(d$name), columns = levels(d$group)))
title_size <- "-log10(FDR)"
}
else if (displayBy == "zscore") {
d <- d[with(d, order(nSig, group, zscore, -adjp)), ]
if (!is.factor(d$name)) {
d$name <- factor(d$name, levels = unique(d$name))
}
d$val <- d$zscore
mat_val <- as.matrix(xSparseMatrix(d[, c("name", "group",
"val")], rows = levels(d$name), columns = levels(d$group)))
title_size <- "Z-score"
}
else if (displayBy == "pvalue") {
d$pvalue[d$pvalue == 0] <- min(d$pvalue[d$pvalue != 0])
d <- d[with(d, order(nSig, group, -pvalue)), ]
if (!is.factor(d$name)) {
d$name <- factor(d$name, levels = unique(d$name))
}
d$val <- -1 * log10(d$pvalue)
mat_val <- as.matrix(xSparseMatrix(d[, c("name", "group",
"val")], rows = levels(d$name), columns = levels(d$group)))
title_size <- "-log10(p-value)"
}
d$bycol <- -1 * log10(d$adjp)
mat_fdr <- as.matrix(xSparseMatrix(d[, c("name", "group",
"adjp")], rows = levels(d$name), columns = levels(d$group)))
mat_fdr[mat_fdr == 0] <- 1
ind_row <- 1:nrow(mat_val)
if (reorder == "row" | reorder == "both") {
ind_row <- match(levels(d$name), rownames(mat_val))
}
ind_row <- rev(ind_row)
ind_col <- 1:ncol(mat_val)
if (reorder == "col" | reorder == "both") {
mat <- mat_val
colnames(mat) <- 1:ncol(mat)
rownames(mat) <- 1:nrow(mat)
tree_bs <- visTreeBootstrap(t(mat), visTree = FALSE)
ind_col <- match(tree_bs$tip.label, colnames(mat))
}
mat_val <- mat_val[ind_row, ind_col]
mat_fdr <- mat_fdr[ind_row, ind_col]
if (is.null(slim)) {
slim <- c(min(mat_val), max(mat_val))
}
mat_val[mat_val <= slim[1]] <- slim[1]
mat_val[mat_val >= slim[2]] <- slim[2]
if (method != "ggplot2") {
if (is.null(title)) {
corrplot::corrplot(mat_val, method = method, is.cor = FALSE,
col = xColormap(colormap)(ncolors), cl.lim = c(zlim[1],
zlim[2]), p.mat = mat_fdr, sig.level = FDR.cutoff,
insig = "blank", addgrid.col = "transparent",
mar = c(0.1, 0.1, 1, 0.1), ...)
}
else {
corrplot::corrplot(mat_val, method = method, is.cor = FALSE,
col = xColormap(colormap)(ncolors), cl.lim = c(zlim[1],
zlim[2]), p.mat = mat_fdr, sig.level = FDR.cutoff,
insig = "blank", addgrid.col = "transparent",
mar = c(0.1, 0.1, 1, 0.1), title = title, ...)
}
invisible(d)
}
else {
if (flip) {
d$group <- factor(d$group, levels = rev(colnames(mat_val)))
d$name <- factor(d$name, levels = rownames(mat_val))
}
else {
d$group <- factor(d$group, levels = colnames(mat_val))
d$name <- factor(d$name, levels = rev(rownames(mat_val)))
}
d$val[d$val <= slim[1]] <- slim[1]
d$val[d$val >= slim[2]] <- slim[2]
if (is.null(zlim)) {
zlim <- c(0, ceiling(max(d$bycol) * 10)/10)
}
d$bycol[d$bycol >= zlim[2]] <- zlim[2]
group <- name <- val <- zscore <- bycol <- NULL
gp <- ggplot(d, aes(x = group, y = name, color = bycol))
gp <- gp + geom_point(aes(size = val), shape = shape)
if (legend.direction == "vertical") {
gp <- gp + scale_size_continuous(limits = c(floor(min(d$val) *
10)/10, ceiling(max(d$val) * 10)/10), range = c(1,
4), guide = guide_legend(title_size, title.position = "top",
ncol = 1, order = 1))
gp <- gp + scale_colour_gradientn(colors =
xColormap(colormap)(ncolors),
limits = zlim, guide = guide_colorbar(title =
expression(-log[10]("FDR")),
title.position = "top", barwidth = 0.5, order = 2))
}
else if (legend.direction == "horizontal") {
gp <- gp + scale_size_continuous(limits = c(floor(min(d$val) *
10)/10, ceiling(max(d$val) * 10)/10), range = c(1,
4), guide = guide_legend(title_size, title.position = "top",
keywidth = 0.5, keyheight = 0.5, ncol = 3, byrow = T,
order = 1))
gp <- gp + scale_colour_gradientn(colors =
xColormap(colormap)(ncolors),
limits = zlim, guide = guide_colorbar(title =
expression(-log[10]("FDR")),
title.position = "top", , barheight = 0.5,
direction = "horizontal", order = 2))
}
gp <- gp + theme_bw() + theme(legend.position = "right",
axis.title.x = element_blank(), axis.title.y = element_blank(),
axis.text.x = element_text(face = "bold", color = "black",
size = 8, angle = y.rotate), axis.text.y = element_text(face = "bold",
color = "black", size = 8, angle = 0), panel.background =
element_rect(fill = "transparent")) +
theme(panel.grid.major = element_blank(), panel.grid.minor =
element_blank())
gp <- gp + labs(title = title)
gp <- gp + theme(text = element_text(family = font.family))
gp <- gp + scale_x_discrete(position = "top")
if (flip) {
gp <- gp + coord_flip() + scale_x_discrete(position = "bottom") +
scale_y_discrete(position = "right")
gp <- gp + theme(axis.text.x = element_text(face = "bold",
color = "black", size = 8, angle = y.rotate,
hjust = 0))
}
invisible(gp)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
