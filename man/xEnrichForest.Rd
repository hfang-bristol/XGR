\name{xEnrichForest}
\alias{xEnrichForest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichForest(eTerm, top_num = 10, FDR.cutoff = 0.05, CI.one = T,
colormap = "ggplot2.top", ncolors = 64, zlim = NULL, barwidth = 0.5,
barheight = NULL, wrap.width = NULL, font.family = "sans", signature =
FALSE, drop = F, sortBy = c("or", "adjp", "fdr", "pvalue", "zscore",
"fc", "nAnno", "nOverlap", "none"))
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
\item{signature}{
%% ~~Describe \code{signature} here~~
}
\item{drop}{
%% ~~Describe \code{drop} here~~
}
\item{sortBy}{
%% ~~Describe \code{sortBy} here~~
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
colormap = "ggplot2.top", ncolors = 64, zlim = NULL, barwidth = 0.5,
barheight = NULL, wrap.width = NULL, font.family = "sans",
signature = FALSE, drop = F, sortBy = c("or", "adjp", "fdr",
"pvalue", "zscore", "fc", "nAnno", "nOverlap", "none"))
{
sortBy <- match.arg(sortBy)
if (is.null(eTerm)) {
warnings("There is no enrichment in the 'eTerm' object.\n")
return(NULL)
}
if (any(class(eTerm) \%in\% "eTerm")) {
df <- xEnrichViewer(eTerm, top_num = "all")
if (!CI.one) {
ind <- which(df$CIl > 1 | df$CIu < 1)
df <- df[ind, ]
}
if (top_num == "auto") {
top_num <- sum(df$adjp < FDR.cutoff)
if (top_num <= 1) {
warnings("There is no enrichment in the 'eTerm' object.\n")
top_num <- 10
message(sprintf("0 or 1 term found; instead the top \%d terms sorted by
'\%s' are shown!",
top_num, sortBy), appendLF = TRUE)
}
}
df <- xEnrichViewer(eTerm, top_num = top_num, sortBy = sortBy)
df$group <- "group"
df$ontology <- "ontology"
}
else if (any(class(eTerm) \%in\% "ls_eTerm") | any(class(eTerm) \%in\%
c("data.frame", "tbl_df", "tbl"))) {
if (any(class(eTerm) \%in\% "ls_eTerm")) {
df <- eTerm$df
}
else if (any(class(eTerm) \%in\% c("data.frame", "tbl_df",
"tbl"))) {
eTerm <- as.data.frame(eTerm)
if (all(c("group", "ontology", "name", "adjp", "or",
"CIl", "CIu") \%in\% colnames(eTerm))) {
df <- eTerm[, c("group", "ontology", "name",
"adjp", "or", "CIl", "CIu")]
}
else if (all(c("group", "name", "adjp", "or", "CIl",
"CIu") \%in\% colnames(eTerm))) {
df <- eTerm[, c("group", "name", "adjp", "or",
"CIl", "CIu")]
df$ontology <- "ontology"
}
else if (all(c("ontology", "name", "adjp", "or",
"CIl", "CIu") \%in\% colnames(eTerm))) {
df <- eTerm[, c("ontology", "name", "adjp", "or",
"CIl", "CIu")]
df$group <- "group"
}
else if (all(c("name", "adjp", "or", "CIl", "CIu") \%in\%
colnames(eTerm))) {
df <- eTerm[, c("name", "adjp", "or", "CIl",
"CIu")]
df$group <- "group"
df$ontology <- "ontology"
}
else {
warnings("The input data.frame does not contain required columns:
c('group','ontology','name','adjp','or','CIl','CIu').\n")
return(NULL)
}
}
df$group <- factor(df$group, levels = unique(df$group))
if (!CI.one) {
ind <- which(df$CIl > 1 | df$CIu < 1)
df <- df[ind, ]
}
or <- group <- ontology <- rank <- adjp <- NULL
df <- df \%>\% dplyr::arrange(-or)
if (top_num == "auto") {
df <- subset(df, df$adjp < FDR.cutoff)
}
else {
top_num <- as.integer(top_num)
df <- as.data.frame(df \%>\% dplyr::group_by(group,
ontology) \%>\% dplyr::group_by(rank = rank(-or),
add = TRUE) \%>\% dplyr::filter(rank <= top_num &
adjp < FDR.cutoff))
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
df$fdr <- -log10(df$adjp)
if (is.null(zlim)) {
tmp <- df$fdr
zlim <- c(floor(min(tmp)), ceiling(max(tmp[!is.infinite(tmp)])))
}
df$fdr[df$fdr <= zlim[1]] <- zlim[1]
df$fdr[df$fdr >= zlim[2]] <- zlim[2]
if (class(eTerm) == "eTerm" & sortBy != "or") {
df <- df[rev(1:nrow(df)), ]
}
else {
df <- df[with(df, order(group, ontology, or, -fdr)),
]
}
df$name <- factor(df$name, levels = unique(df$name))
bp <- ggplot(df, aes(x = name, y = log2(or), ymin = log2(CIl),
ymax = log2(CIu), color = fdr))
bp <- bp + geom_pointrange(size = 0.3) + ylab(expression(log[2]("odds
ratio")))
bp <- bp + geom_hline(yintercept = 0, color = "black", linetype =
"dashed") +
coord_flip()
bp <- bp + theme_bw() + theme(legend.position = "right",
axis.title.y = element_blank(), axis.text.y = element_text(size = 8,
color = "black"), axis.title.x = element_text(size = 10,
color = "black"))
bp <- bp + theme(panel.grid.major = element_blank(), panel.grid.minor =
element_blank())
bp <- bp + scale_colour_gradientn(colors =
xColormap(colormap)(ncolors),
limits = zlim, guide = guide_colorbar(title =
expression(-log[10]("FDR")),
title.position = "top", barwidth = barwidth, barheight = barheight,
draw.ulim = FALSE, draw.llim = FALSE))
if (signature) {
caption <- paste("Created by xEnrichForest from XGR version",
utils::packageVersion("XGR"))
bp <- bp + labs(caption = caption) + theme(plot.caption =
element_text(hjust = 1,
face = "bold.italic", size = 8, colour = "#002147"))
}
bp <- bp + theme(text = element_text(family = font.family))
bp <- bp + theme(axis.line.x = element_line(arrow = arrow(angle = 30,
length = unit(0.25, "cm"), type = "open")))
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
