\name{xEnrichGGraph}
\alias{xEnrichGGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichGGraph(eTerm, ig = NULL, fixed = T, node.color = c("zscore",
"adjp", "or"), colormap = "grey-orange-darkred", zlim = NULL, node.size
= c("adjp", "zscore", "or"), slim = NULL, node.size.range = c(0.5, 4),
node.label.size = 2, leave = T, ncolumns = NULL, ...)
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
\item{node.label.size}{
%% ~~Describe \code{node.label.size} here~~
}
\item{leave}{
%% ~~Describe \code{leave} here~~
}
\item{ncolumns}{
%% ~~Describe \code{ncolumns} here~~
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
function (eTerm, ig = NULL, fixed = T, node.color = c("zscore",
"adjp", "or"), colormap = "grey-orange-darkred", zlim = NULL,
node.size = c("adjp", "zscore", "or"), slim = NULL, node.size.range =
c(0.5,
4), node.label.size = 2, leave = T, ncolumns = NULL,
...)
{
node.color <- match.arg(node.color)
node.size <- match.arg(node.size)
if (is.null(eTerm)) {
warnings("There is no enrichment in the 'eTerm' object.\n")
return(NULL)
}
if (0) {
if (class(ig) != "igraph") {
warnings("The 'ig' object must be provided.\n")
return(NULL)
}
else {
if (!("name" \%in\% igraph::vertex_attr_names(ig))) {
warnings("The 'ig' object must contain a node attribute 'name'.\n")
return(NULL)
}
}
}
if (class(eTerm) == "eTerm") {
df_enrichment_group <- xEnrichViewer(eTerm, top_num = "all")
df_enrichment_group$group <- "group"
}
else if (class(eTerm) == "ls_eTerm" | class(eTerm) == "data.frame") {
if (class(eTerm) == "ls_eTerm") {
df_enrichment_group <- eTerm$df[, c("group", "name",
"adjp", "or", "zscore")]
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
V(ig)$name <- V(ig)$term_name
}
else {
return(NULL)
}
}
if (!fixed) {
ind <- match(V(ig)$name, df_enrichment_group$name)
nodes_query <- V(ig)$name[!is.na(ind)]
if (class(suppressWarnings(try(ig <- dnet::dDAGinduce(ig,
nodes_query, path.mode = "all_paths"), T))) == "try-error") {
ig <- NULL
}
}
gp <- NULL
if (class(ig) == "igraph") {
df_enrichment_group$or[is.infinite(df_enrichment_group$or)] <-
max(df_enrichment_group$or[!is.infinite(df_enrichment_group$or)])
ls_df <- split(x = df_enrichment_group[, c("name", "zscore",
"adjp", "or")], f = df_enrichment_group$group)
ls_igg <- lapply(ls_df, function(df_enrichment) {
igg <- ig
V(igg)$zscore <- 0
ind <- match(V(igg)$name, df_enrichment$name)
V(igg)$zscore[!is.na(ind)] <- df_enrichment$zscore[ind[!is.na(ind)]]
V(igg)$adjp <- 0
ind <- match(V(igg)$name, df_enrichment$name)
V(igg)$adjp[!is.na(ind)] <- -1 *
log10(df_enrichment$adjp)[ind[!is.na(ind)]]
V(igg)$or <- 0
ind <- match(V(igg)$name, df_enrichment$name)
V(igg)$or[!is.na(ind)] <- log2(df_enrichment$or)[ind[!is.na(ind)]]
igg
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
x <- y <- leaf <- x_group <- y_group <- group <- group_id <- color <-
size <- NULL
if (length(ls_igg) > 1) {
gp_template <- xGGraph(ig, leave = leave, node.label.size =
node.label.size,
...)
gp <- xGGraph(ig, leave = leave, node.label.size = 0,
node.size = 0, ...)
}
else {
gp <- xGGraph(ig, leave = leave, node.label.size = node.label.size,
node.size = 0, ...)
}
if (leave & !is.null(gp$data$leaf)) {
df_data <- subset(gp$data, leaf == T)
}
else {
df_data <- gp$data
}
ls_df <- lapply(1:length(ls_igg), function(i) {
g <- ls_igg[[i]]
group <- names(ls_igg)[i]
ind <- match(V(g)$name, df_data$name)
if (node.color == "or") {
color <- V(g)$or[!is.na(ind)]
}
else if (node.color == "adjp") {
color <- V(g)$adjp[!is.na(ind)]
}
else if (node.color == "zscore") {
color <- V(g)$zscore[!is.na(ind)]
}
color[color <= zlim[1]] <- zlim[1]
color[color >= zlim[2]] <- zlim[2]
if (node.size == "or") {
size <- V(g)$or[!is.na(ind)]
}
else if (node.size == "adjp") {
size <- V(g)$adjp[!is.na(ind)]
}
else if (node.size == "zscore") {
size <- V(g)$zscore[!is.na(ind)]
}
size[size <= slim[1]] <- slim[1]
size[size >= slim[2]] <- slim[2]
df <- data.frame(df_data[ind[!is.na(ind)], ], color = color,
size = size, group = group, group_id = i, stringsAsFactors = F)
})
df <- do.call(rbind, ls_df)
gp <- gp + geom_point(data = df, aes(x = x * 1, y = y *
1, color = color, size = size))
gp <- gp + scale_size_continuous(limits = slim, range =
node.size.range,
guide = guide_legend(node.size.title, title.position = "top",
ncol = 1)) + scale_colour_gradientn(colors = xColormap(colormap)(64),
limits = zlim, guide = guide_colorbar(title = node.color.title,
title.position = "top", barwidth = 0.5))
if (length(ls_igg) > 1) {
if (is.null(ncolumns)) {
ncolumns <- ceiling(sqrt(length(ls_igg)))
}
group <- NULL
gp <- gp + facet_wrap(~group, ncol = ncolumns)
gp <- gp + theme(strip.background = element_rect(fill = "transparent",
color = "transparent"), strip.text = element_text(size = 12,
face = "bold"), strip.placement = "inside", panel.spacing = unit(0,
"lines"))
gp$gp_template <- gp_template
}
gp$data_enrichment <- df
}
return(gp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
