\name{xEnrichChord}
\alias{xEnrichChord}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichChord(eTerm, top_num = 5, FDR.cutoff = 0.05, colormap.group =
"ggplot2", colormap.ontology = NULL, wrap.width = NULL, text.size =
0.6, legend = NULL, vline = F, ...)
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
\item{colormap.group}{
%% ~~Describe \code{colormap.group} here~~
}
\item{colormap.ontology}{
%% ~~Describe \code{colormap.ontology} here~~
}
\item{wrap.width}{
%% ~~Describe \code{wrap.width} here~~
}
\item{text.size}{
%% ~~Describe \code{text.size} here~~
}
\item{legend}{
%% ~~Describe \code{legend} here~~
}
\item{vline}{
%% ~~Describe \code{vline} here~~
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
function (eTerm, top_num = 5, FDR.cutoff = 0.05, colormap.group =
"ggplot2",
colormap.ontology = NULL, wrap.width = NULL, text.size = 0.6,
legend = NULL, vline = F, ...)
{
if (is.null(eTerm)) {
warnings("There is no enrichment in the 'eTerm' object.\n")
return(NULL)
}
if (class(eTerm) == "eTerm") {
df <- xEnrichViewer(eTerm, top_num = "all")
if (top_num == "auto") {
top_num <- sum(df$adjp < FDR.cutoff)
if (top_num <= 1) {
top_num <- 10
}
}
df <- xEnrichViewer(eTerm, top_num = top_num, sortBy = "adjp")
df$group <- "group"
df$ontology <- "ontology"
}
else if (class(eTerm) == "ls_eTerm" | class(eTerm) == "data.frame") {
if (class(eTerm) == "ls_eTerm") {
df <- eTerm$df[, c("group", "ontology", "name", "adjp",
"zscore")]
}
else if (class(eTerm) == "data.frame") {
if (all(c("group", "ontology", "name", "adjp", "zscore") \%in\%
colnames(eTerm))) {
df <- eTerm[, c("group", "ontology", "name",
"adjp", "zscore")]
}
else if (all(c("group", "name", "adjp", "zscore") \%in\%
colnames(eTerm))) {
df <- eTerm[, c("group", "name", "adjp", "zscore")]
df$ontology <- "ontology"
}
else if (all(c("ontology", "name", "adjp", "zscore") \%in\%
colnames(eTerm))) {
df <- eTerm[, c("ontology", "name", "adjp", "zscore")]
df$group <- "group"
}
else if (all(c("name", "adjp", "zscore") \%in\% colnames(eTerm))) {
df <- eTerm[, c("name", "adjp", "zscore")]
df$group <- "group"
df$ontology <- "ontology"
}
else {
warnings("The input data.frame does not contain required columns:
c('group','ontology','name','adjp','zscore').\n")
return(NULL)
}
}
group <- ontology <- rank <- adjp <- NULL
df <- df \%>\% dplyr::arrange(adjp)
if (top_num == "auto") {
df <- subset(df, df$adjp < FDR.cutoff)
}
else {
top_num <- as.integer(top_num)
df <- as.data.frame(df \%>\% dplyr::group_by(group,
ontology) \%>\% dplyr::group_by(rank = rank(adjp),
add = TRUE) \%>\% dplyr::filter(rank <= top_num &
adjp < FDR.cutoff))
}
}
if (!is.null(wrap.width)) {
width <- as.integer(wrap.width)
res_list <- lapply(df$name, function(x) {
x <- gsub("_", " ", x)
y <- strwrap(x, width = width)
if (length(y) > 1) {
paste(y, collapse = "\n")
}
else {
y
}
})
df$name <- unlist(res_list)
}
name <- n <- ontology <- zscore <- group <- NULL
df <- as.data.frame(df \%>\% dplyr::group_by(name) \%>\%
dplyr::group_by(n = n(),
add = T) \%>\% dplyr::arrange(ontology, n, zscore))
order <- c(sort(unique(df$group)), unique(df$name))
if (is.null(colormap.ontology)) {
if (length(table(df$ontology)) == 1) {
colormap.ontology <- "grey-grey"
}
else {
colormap.ontology <- "spectral"
}
}
circlize::circos.clear()
circlize::circos.par(start.degree = 90, clock.wise = FALSE)
tmp <- table(df$group)
color_group <- xColormap(colormap.group)(length(tmp))
names(color_group) <- names(tmp)
tmp <- table(df$ontology)
color_ontology <- xColormap(colormap.ontology)(length(tmp))
names(color_ontology) <- names(tmp)
color_name <- color_ontology[df$ontology]
names(color_name) <- df$name
grid.col <- c(color_group, color_name)
df_data <- df[, c("group", "name", "zscore")]
circlize::chordDiagram(df_data, annotationTrack = "grid",
grid.col = grid.col, preAllocateTracks = 1, ...)
circlize::circos.trackPlotRegion(track.index = 1, panel.fun =
function(x,
y) {
xlim = circlize::get.cell.meta.data("xlim")
ylim = circlize::get.cell.meta.data("ylim")
sector.name = circlize::get.cell.meta.data("sector.index")
circlize::circos.text(mean(xlim), ylim[1] + 0.1, sector.name,
facing = "clockwise", niceFacing = TRUE, adj = c(0,
0.5), cex = text.size)
}, bg.border = NA)
if (vline) {
abline(v = 0, lty = 2, col = "#00000080")
}
if (any((is.null(legend) & length(color_group) > 1), legend)) {
legend("topleft", title = "Groups", legend = names(color_group),
border = color_group, fill = color_group, horiz = F,
box.col = "transparent", cex = text.size)
}
if (any((is.null(legend) & length(color_ontology) > 1), legend)) {
legend("topright", title = "Ontologies", legend =
names(color_ontology),
border = color_ontology, fill = color_ontology, horiz = F,
box.col = "transparent", cex = text.size)
}
invisible(df)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
