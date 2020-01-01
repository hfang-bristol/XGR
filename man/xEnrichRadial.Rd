\name{xEnrichRadial}
\alias{xEnrichRadial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichRadial(eTerm, ig = NULL, fixed = T, node.color = c("or", "adjp",
"zscore"), colormap = "grey-orange-darkred", zlim = NULL, node.size =
c("adjp", "or", "zscore"), slim = NULL, node.size.range = c(0.5, 3.5),
edge.color = "skyblue", edge.color.alpha = 0.5, edge.curve = 0.1,
edge.arrow.gap = 0.02, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{eTerm}{
%% ~~Describe \code{eTerm} here~~
}
\item{ig}{
%% ~~Describe \code{ig} here~~
}
\item{fixed}{
%% ~~Describe \code{fixed} here~~
}
\item{node.color}{
%% ~~Describe \code{node.color} here~~
}
\item{colormap}{
%% ~~Describe \code{colormap} here~~
}
\item{zlim}{
%% ~~Describe \code{zlim} here~~
}
\item{node.size}{
%% ~~Describe \code{node.size} here~~
}
\item{slim}{
%% ~~Describe \code{slim} here~~
}
\item{node.size.range}{
%% ~~Describe \code{node.size.range} here~~
}
\item{edge.color}{
%% ~~Describe \code{edge.color} here~~
}
\item{edge.color.alpha}{
%% ~~Describe \code{edge.color.alpha} here~~
}
\item{edge.curve}{
%% ~~Describe \code{edge.curve} here~~
}
\item{edge.arrow.gap}{
%% ~~Describe \code{edge.arrow.gap} here~~
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
function (eTerm, ig = NULL, fixed = T, node.color = c("or", "adjp",
"zscore"), colormap = "grey-orange-darkred", zlim = NULL,
node.size = c("adjp", "or", "zscore"), slim = NULL, node.size.range =
c(0.5,
3.5), edge.color = "skyblue", edge.color.alpha = 0.5,
edge.curve = 0.1, edge.arrow.gap = 0.02, ...)
{
node.color <- match.arg(node.color)
node.size <- match.arg(node.size)
gp_data <- NULL
gp_code <- NULL
gp_table <- NULL
if (is.null(eTerm)) {
warnings("There is no enrichment in the 'eTerm' object.\n")
return(NULL)
}
if (class(eTerm) == "eTerm") {
df_enrichment_group <- xEnrichViewer(eTerm, top_num = "all",
sortBy = "or")
df_enrichment_group$group <- "group"
}
else if (class(eTerm) == "ls_eTerm" | class(eTerm) == "data.frame") {
if (class(eTerm) == "ls_eTerm") {
df_enrichment_group <- eTerm$df
}
else if (class(eTerm) == "data.frame") {
if (all(c("group", "name", "adjp", "or", "zscore") \%in\%
colnames(eTerm))) {
df_enrichment_group <- eTerm[, c("group", "name",
"adjp", "or", "zscore")]
}
else if (all(c("name", "adjp", "or", "zscore") \%in\%
colnames(eTerm))) {
df_enrichment_group <- eTerm[, c("name", "adjp",
"or", "zscore")]
df_enrichment_group$group <- "group"
}
}
}
if (class(df_enrichment_group$group) == "factor") {
if (length(unique(df_enrichment_group$group)) !=
length(levels(df_enrichment_group$group))) {
df_enrichment_group$group <- factor(df_enrichment_group$group,
levels = sort(unique(df_enrichment_group$group)))
}
}
if (class(ig) != "igraph") {
if (class(eTerm) == "eTerm") {
ig <- eTerm$g
}
else {
return(NULL)
}
}
if (!fixed) {
ind <- match(V(ig)$term_name, df_enrichment_group$name)
nodes_query <- V(ig)$name[!is.na(ind)]
if (class(suppressWarnings(try(subg <- dnet::dDAGinduce(ig,
nodes_query, path.mode = "all_paths"), T))) == "try-error") {
subg <- NULL
}
}
else {
subg <- ig
}
if (class(subg) == "igraph") {
gp_code_table <- xOBOcode(g = subg, node.level = "term_distance",
node.level.value = 2, node.label.color = "black",
node.shape = 21, node.size.range = 4, edge.color.alpha = 0.2,
table.base.size = 7, table.row.space = 2, table.nrow =
min(vcount(subg),
55), ...)
gp_code <- gp_code_table$code
gp_table <- gp_code_table$table
df_enrichment_group$or[is.infinite(df_enrichment_group$or)] <-
max(df_enrichment_group$or[!is.infinite(df_enrichment_group$or)])
ls_df <- split(x = df_enrichment_group[, c("name", "zscore",
"adjp", "or")], f = df_enrichment_group$group)
ls_ig <- lapply(ls_df, function(df_enrichment) {
V(subg)$zscore <- 0
ind <- match(V(subg)$term_name, df_enrichment$name)
V(subg)$zscore[!is.na(ind)] <- df_enrichment$zscore[ind[!is.na(ind)]]
V(subg)$adjp <- 0
ind <- match(V(subg)$term_name, df_enrichment$name)
V(subg)$adjp[!is.na(ind)] <- -1 *
log10(df_enrichment$adjp)[ind[!is.na(ind)]]
V(subg)$or <- 0
ind <- match(V(subg)$term_name, df_enrichment$name)
V(subg)$or[!is.na(ind)] <- log2(df_enrichment$or)[ind[!is.na(ind)]]
subg
})
if (node.color == "or") {
node.color.title <- expression(log[2](OR))
if (is.null(zlim)) {
zlim <- c(0, ceiling(max(log2(df_enrichment_group$or))))
}
}
else if (node.color == "adjp") {
node.color.title <- expression(-log[10](FDR))
if (is.null(zlim)) {
zlim <- c(0, ceiling(max(-1 * log10(df_enrichment_group$adjp))))
}
}
else if (node.color == "zscore") {
node.color.title <- "Z-score"
if (is.null(zlim)) {
zlim <- c(0, ceiling(max(df_enrichment_group$zscore)))
}
}
if (node.size == "or") {
node.size.title <- expression(log[2](OR))
if (is.null(slim)) {
slim <- c(0, ceiling(max(log2(df_enrichment_group$or))))
}
}
else if (node.size == "adjp") {
node.size.title <- expression(-log[10](FDR))
if (is.null(slim)) {
slim <- c(0, ceiling(max(-1 * log10(df_enrichment_group$adjp))))
}
}
else if (node.size == "zscore") {
node.size.title <- "Z-score"
if (is.null(slim)) {
slim <- c(0, ceiling(max(df_enrichment_group$zscore)))
}
}
gp_data <- xGGnetwork(g = ls_ig, node.color = node.color,
node.color.title = node.color.title, colormap = colormap,
zlim = zlim, node.size = node.size, node.size.title = node.size.title,
slim = slim, node.size.range = node.size.range, edge.color =
edge.color,
edge.color.alpha = edge.color.alpha, edge.curve = edge.curve,
edge.arrow.gap = edge.arrow.gap, ...)
}
ls_gp <- list(code = gp_code, table = gp_table, data = gp_data,
ig = ls_ig)
return(ls_gp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
