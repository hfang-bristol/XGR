\name{xEnrichD3}
\alias{xEnrichD3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichD3(eTerm, top_num = 10, FDR.cutoff = 0.05, type = c("sankey",
"force", "radial", "diagonal"), colormap = "ggplot2", filename =
"xEnrichD3", ...)
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
\item{type}{
%% ~~Describe \code{type} here~~
}
\item{colormap}{
%% ~~Describe \code{colormap} here~~
}
\item{filename}{
%% ~~Describe \code{filename} here~~
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
function (eTerm, top_num = 10, FDR.cutoff = 0.05, type = c("sankey",
"force", "radial", "diagonal"), colormap = "ggplot2", filename =
"xEnrichD3",
...)
{
type <- match.arg(type)
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
df <- eTerm$df
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
df_nodes_group <- data.frame(name = unique(df$group), category =
"group",
stringsAsFactors = FALSE)
df_nodes_ontology <- unique(data.frame(name = df$name, category =
df$ontology,
stringsAsFactors = FALSE))
df_nodes <- rbind(df_nodes_group, df_nodes_ontology)
tmp <- table(df_nodes$category)
color_category <- xColormap(colormap)(length(tmp))
names(color_category) <- names(tmp)
df_nodes$color <- color_category[df_nodes$category]
df_edges <- df[, c("group", "name", "zscore")]
g <- igraph::graph.data.frame(d = df_edges, directed = TRUE,
vertices = df_nodes)
if (type == "sankey") {
d3 <- networkD3::igraph_to_networkD3(g, group = df_nodes$category)
colors <- paste(color_category, collapse = "\", \"")
colourScale <- paste("d3.scaleOrdinal([\"", colors, "\"])")
res <- networkD3::sankeyNetwork(Links = d3$links, Nodes = d3$nodes,
Source = "source", Target = "target", Value = "value",
NodeID = "name", NodeGroup = "group", colourScale = colourScale,
fontSize = 12, sinksRight = F, ...)
}
else if (type == "force") {
d3 <- networkD3::igraph_to_networkD3(g, group = df_nodes$category)
d3 <- networkD3::igraph_to_networkD3(g, group = df_nodes$category)
d3$links$value <- 1 + 4 * (d3$links$value -
min(d3$links$value))/(max(d3$links$value) -
min(d3$links$value))
linkWidth <- networkD3::JS("function(d) { return d.value; }")
colors <- paste(color_category, collapse = "\", \"")
colourScale <- paste("d3.scaleOrdinal([\"", colors, "\"])")
res <- networkD3::forceNetwork(Links = d3$links, Nodes = d3$nodes,
Source = "source", Target = "target", Value = "value",
NodeID = "name", Group = "group", colourScale = colourScale,
zoom = T, linkDistance = 200, linkWidth = linkWidth,
...)
}
else if (type \%in\% c("radial", "diagonal")) {
if (0) {
}
data <- xConverter(g, from = "igraph", to = "lol", verbose = F)
if (type == "radial") {
res <- networkD3::radialNetwork(data, ...)
}
else {
res <- networkD3::diagonalNetwork(data, ...)
}
}
if (all(!is.null(filename), !is.na(filename), filename !=
"")) {
filename <- gsub(".html$", "", filename)
filename <- paste0(filename, ".html")
res \%>\% networkD3::saveNetwork(file = filename)
}
res$ig <- g
invisible(res)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
