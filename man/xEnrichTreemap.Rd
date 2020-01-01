\name{xEnrichTreemap}
\alias{xEnrichTreemap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichTreemap(eTerm, top_num = 10, FDR.cutoff = 0.05, CI.one = T,
colormap = "spectral.top", ncolors = 64, zlim = NULL, barwidth = NULL,
barheight = 0.5, wrap.width = NULL, font.family = "sans", drop = F,
details = c("name", "name_FDR", "name_FDR_members"), caption = T,
treemap.grow = F, treemap.reflow = F, treemap.place = "topleft",
treemap.color = "black", treemap.fontface = "bold.italic",
treemap.min.size = 4, area = c("adjp", "or"), area.fill = c("or",
"adjp", "zscore"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{eTerm}{
%% ~~Describe \code{eTerm} here~~
}
\item{top_num}{
%% ~~Describe \code{top_num} here~~
}
\item{FDR.cutoff}{
%% ~~Describe \code{FDR.cutoff} here~~
}
\item{CI.one}{
%% ~~Describe \code{CI.one} here~~
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
\item{wrap.width}{
%% ~~Describe \code{wrap.width} here~~
}
\item{font.family}{
%% ~~Describe \code{font.family} here~~
}
\item{drop}{
%% ~~Describe \code{drop} here~~
}
\item{details}{
%% ~~Describe \code{details} here~~
}
\item{caption}{
%% ~~Describe \code{caption} here~~
}
\item{treemap.grow}{
%% ~~Describe \code{treemap.grow} here~~
}
\item{treemap.reflow}{
%% ~~Describe \code{treemap.reflow} here~~
}
\item{treemap.place}{
%% ~~Describe \code{treemap.place} here~~
}
\item{treemap.color}{
%% ~~Describe \code{treemap.color} here~~
}
\item{treemap.fontface}{
%% ~~Describe \code{treemap.fontface} here~~
}
\item{treemap.min.size}{
%% ~~Describe \code{treemap.min.size} here~~
}
\item{area}{
%% ~~Describe \code{area} here~~
}
\item{area.fill}{
%% ~~Describe \code{area.fill} here~~
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
function (eTerm, top_num = 10, FDR.cutoff = 0.05, CI.one = T,
colormap = "spectral.top", ncolors = 64, zlim = NULL, barwidth = NULL,
barheight = 0.5, wrap.width = NULL, font.family = "sans",
drop = F, details = c("name", "name_FDR", "name_FDR_members"),
caption = T, treemap.grow = F, treemap.reflow = F, treemap.place =
"topleft",
treemap.color = "black", treemap.fontface = "bold.italic",
treemap.min.size = 4, area = c("adjp", "or"), area.fill = c("or",
"adjp", "zscore"))
{
details <- match.arg(details)
area <- match.arg(area)
area.fill <- match.arg(area.fill)
if (is.null(eTerm)) {
warnings("There is no enrichment in the 'eTerm' object.\n")
return(NULL)
}
if (class(eTerm) == "eTerm") {
df <- xEnrichViewer(eTerm, top_num = "all")
if (!CI.one) {
ind <- which(df$CIl > 1 | df$CIu < 1)
df <- df[ind, ]
}
if (top_num == "auto") {
top_num <- sum(df$adjp < FDR.cutoff)
if (top_num <= 1) {
top_num <- 10
}
}
df <- xEnrichViewer(eTerm, top_num = top_num, sortBy = "or",
details = T)
df$group <- "group"
df$ontology <- "ontology"
}
else if (class(eTerm) == "ls_eTerm" | class(eTerm) == "data.frame") {
if (class(eTerm) == "ls_eTerm") {
df <- eTerm$df
}
else if (class(eTerm) == "data.frame") {
if (all(c("group", "ontology", "name", "zscore",
"adjp", "or", "CIl", "CIu", "nOverlap", "members_Overlap") \%in\%
colnames(eTerm))) {
df <- eTerm[, c("group", "ontology", "name",
"zscore", "adjp", "or", "CIl", "CIu", "nOverlap",
"members_Overlap")]
}
else {
details <- "name"
if (all(c("group", "ontology", "name", "zscore",
"adjp", "or", "CIl", "CIu") \%in\% colnames(eTerm))) {
df <- eTerm[, c("group", "ontology", "name",
"zscore", "adjp", "or", "CIl", "CIu")]
}
else if (all(c("group", "name", "zscore", "adjp",
"or", "CIl", "CIu") \%in\% colnames(eTerm))) {
df <- eTerm[, c("group", "name", "zscore",
"adjp", "or", "CIl", "CIu")]
df$ontology <- "ontology"
}
else if (all(c("ontology", "name", "zscore",
"adjp", "or", "CIl", "CIu") \%in\% colnames(eTerm))) {
df <- eTerm[, c("ontology", "name", "zscore",
"adjp", "or", "CIl", "CIu")]
df$group <- "group"
}
else if (all(c("name", "zscore", "adjp", "or",
"CIl", "CIu") \%in\% colnames(eTerm))) {
df <- eTerm[, c("name", "zscore", "adjp", "or",
"CIl", "CIu")]
df$group <- "group"
df$ontology <- "ontology"
}
else {
warnings("The input data.frame does not contain required columns:
c('group','ontology','name','zscore','adjp','or','CIl','CIu').\n")
return(NULL)
}
}
}
df$group <- factor(df$group, levels = unique(df$group))
if (!CI.one) {
ind <- which(df$CIl > 1 | df$CIu < 1)
df <- df[ind, ]
}
or <- group <- ontology <- rank <- NULL
df <- df \%>\% dplyr::arrange(-or)
if (top_num == "auto") {
df <- subset(df, df$adjp < FDR.cutoff)
}
else {
top_num <- as.integer(top_num)
df_tmp <- as.data.frame(df \%>\% dplyr::group_by(group,
ontology) \%>\% dplyr::group_by(rank = rank(-or,
decreasing = T), add = TRUE) \%>\% dplyr::filter(rank <=
top_num))
df <- subset(df, df$name \%in\% df_tmp$name)
df <- subset(df, df$adjp < FDR.cutoff)
}
}
if (nrow(df) == 0) {
return(NULL)
}
if (!is.null(wrap.width)) {
width <- as.integer(wrap.width)
res_list <- lapply(df$name, function(x) {
x <- gsub("_", " ", x)
y <- strwrap(x, width = width)
if (length(y) > 1) {
paste0(y[1], "...")
}
else {
y
}
})
df$name <- unlist(res_list)
}
name <- fdr <- or <- CIl <- CIu <- NULL
group <- ontology <- NULL
label <- NULL
lor <- zscore <- NULL
df$fdr <- -log10(df$adjp)
df$lor <- log2(df$or)
if (area.fill == "adjp") {
title <- expression(-log[10]("FDR"))
if (is.null(zlim)) {
tmp <- df$fdr
zlim <- c(floor(min(tmp)), ceiling(max(tmp[!is.infinite(tmp)])))
}
df$fdr[df$fdr <= zlim[1]] <- zlim[1]
df$fdr[df$fdr >= zlim[2]] <- zlim[2]
}
else if (area.fill == "or") {
title <- expression(log[2]("odds ratio"))
if (is.null(zlim)) {
tmp <- df$lor
zlim <- c(floor(min(tmp)), ceiling(max(tmp[!is.infinite(tmp)])))
}
df$lor[df$lor <= zlim[1]] <- zlim[1]
df$lor[df$lor >= zlim[2]] <- zlim[2]
}
else if (area.fill == "zscore") {
title <- "Z-score"
if (is.null(zlim)) {
tmp <- df$zscore
zlim <- c(floor(min(tmp)), ceiling(max(tmp[!is.infinite(tmp)])))
}
df$zscore[df$zscore <= zlim[1]] <- zlim[1]
df$zscore[df$zscore >= zlim[2]] <- zlim[2]
}
df <- df[with(df, order(group, ontology, or, fdr)), ]
df$name <- factor(df$name, levels = unique(df$name))
if (details == "name") {
df$label <- df$name
}
else if (details == "name_FDR") {
df$label <- paste0(df$name, "\n[OR=", df$or, "; FDR=",
df$adjp, "; n=", df$nOverlap, "]")
}
else if (details == "name_FDR_members") {
treemap.grow <- T
treemap.reflow <- T
df$label <- paste0(df$name, "\n[OR=", df$or, "; FDR=",
df$adjp, "; n=", df$nOverlap, "]\n(", df$members_Overlap,
")")
}
if (area.fill == "adjp") {
bp <- ggplot(df, aes(area = lor, fill = fdr, label = label))
}
else if (area.fill == "or") {
bp <- ggplot(df, aes(area = fdr, fill = lor, label = label))
}
else if (area.fill == "zscore") {
if (area == "adjp") {
bp <- ggplot(df, aes(area = fdr, fill = zscore, label = label))
}
else if (area == "or") {
bp <- ggplot(df, aes(area = lor, fill = zscore, label = label))
}
}
bp <- bp + treemapify::geom_treemap() +
treemapify::geom_treemap_text(fontface = treemap.fontface,
color = treemap.color, place = treemap.place, grow = treemap.grow,
reflow = treemap.reflow, min.size = treemap.min.size)
bp <- bp + theme_bw() + theme(legend.position = "bottom")
bp <- bp + scale_fill_gradientn(colors = xColormap(colormap)(ncolors),
limits = zlim, guide = guide_colorbar(title = title,
title.position = "left", barwidth = barwidth, barheight = barheight,
draw.ulim = FALSE, draw.llim = FALSE))
if (caption) {
if (area == "adjp") {
bp <- bp + labs(caption = "The area is proportional to FDR")
}
else if (area == "or") {
bp <- bp + labs(caption = "The area is proportional to odds ratio")
}
bp <- bp + theme(plot.caption = element_text(hjust = 1,
face = "bold.italic", size = 8, colour = "#002147"))
}
bp <- bp + theme(text = element_text(family = font.family))
ngroup <- length(unique(df$group))
nonto <- length(unique(df$ontology))
if (ngroup != 1 | nonto != 1) {
scales <- "free_y"
space <- "free_y"
if (ngroup == 1) {
bp <- bp + facet_grid(ontology ~ ., scales = scales,
space = space, drop = drop)
}
else if (nonto == 1) {
bp <- bp + facet_grid(. ~ group, scales = scales,
space = space, drop = drop)
}
else {
bp <- bp + facet_grid(ontology ~ group, scales = scales,
space = space, drop = drop)
}
bp <- bp + theme(strip.background = element_rect(fill = "transparent",
color = "transparent"), strip.text = element_text(size = 8,
face = "bold.italic"))
}
invisible(bp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
