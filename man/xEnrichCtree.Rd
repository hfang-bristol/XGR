\name{xEnrichCtree}
\alias{xEnrichCtree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichCtree(eTerm, ig, FDR.cutoff = NULL, node.color = c("zscore",
"adjp", "or", "nOverlap"), colormap = "brewer.Reds", zlim = NULL,
node.size = c("adjp", "zscore", "or", "nOverlap"), slim = NULL,
node.size.range = c(0.5, 4.5), group.gap = 0.08, group.color =
"lightblue", group.size = 0.2, group.label.size = 2, group.label.color
= "black", legend.direction = c("auto", "horizontal", "vertical"),
leave.label.orientation = c("inwards", "outwards"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{eTerm}{
%% ~~Describe \code{eTerm} here~~
}
\item{ig}{
%% ~~Describe \code{ig} here~~
}
\item{FDR.cutoff}{
%% ~~Describe \code{FDR.cutoff} here~~
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
\item{group.gap}{
%% ~~Describe \code{group.gap} here~~
}
\item{group.color}{
%% ~~Describe \code{group.color} here~~
}
\item{group.size}{
%% ~~Describe \code{group.size} here~~
}
\item{group.label.size}{
%% ~~Describe \code{group.label.size} here~~
}
\item{group.label.color}{
%% ~~Describe \code{group.label.color} here~~
}
\item{legend.direction}{
%% ~~Describe \code{legend.direction} here~~
}
\item{leave.label.orientation}{
%% ~~Describe \code{leave.label.orientation} here~~
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
function (eTerm, ig, FDR.cutoff = NULL, node.color = c("zscore",
"adjp", "or", "nOverlap"), colormap = "brewer.Reds", zlim = NULL,
node.size = c("adjp", "zscore", "or", "nOverlap"), slim = NULL,
node.size.range = c(0.5, 4.5), group.gap = 0.08, group.color =
"lightblue",
group.size = 0.2, group.label.size = 2, group.label.color = "black",
legend.direction = c("auto", "horizontal", "vertical"),
leave.label.orientation = c("inwards",
"outwards"), ...)
{
node.color <- match.arg(node.color)
node.size <- match.arg(node.size)
legend.direction <- match.arg(legend.direction)
leave.label.orientation <- match.arg(leave.label.orientation)
if (is.null(eTerm)) {
warnings("There is no enrichment in the 'eTerm' object.\n")
return(NULL)
}
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
if (class(eTerm) == "eTerm") {
df_enrichment_group <- xEnrichViewer(eTerm, top_num = "all")
df_enrichment_group$group <- "group"
}
else if (class(eTerm) == "ls_eTerm" | class(eTerm) == "data.frame") {
if (class(eTerm) == "ls_eTerm") {
df_enrichment_group <- eTerm$df
}
else if (class(eTerm) == "data.frame") {
if (all(c("group", "name", "adjp", "or", "zscore",
"nOverlap") \%in\% colnames(eTerm))) {
df_enrichment_group <- eTerm[, c("group", "name",
"adjp", "or", "zscore", "nOverlap")]
}
else if (all(c("name", "adjp", "or", "zscore", "nOverlap") \%in\%
colnames(eTerm))) {
df_enrichment_group <- eTerm[, c("name", "adjp",
"or", "zscore", "nOverlap")]
df_enrichment_group$group <- "group"
}
}
}
if (!is.null(FDR.cutoff)) {
df_enrichment_group$adjp[df_enrichment_group$adjp >=
FDR.cutoff] <- 1
}
if (class(df_enrichment_group$group) == "factor") {
if (length(unique(df_enrichment_group$group)) !=
length(levels(df_enrichment_group$group))) {
df_enrichment_group$group <- factor(df_enrichment_group$group,
levels = sort(unique(df_enrichment_group$group)))
}
}
gp <- NULL
if (class(ig) == "igraph") {
df_enrichment_group$or[is.infinite(df_enrichment_group$or)] <-
max(df_enrichment_group$or[!is.infinite(df_enrichment_group$or)])
ls_df <- split(x = df_enrichment_group[, c("name", "zscore",
"adjp", "or", "nOverlap")], f = df_enrichment_group$group)
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
V(igg)$nOverlap <- 0
ind <- match(V(igg)$name, df_enrichment$name)
V(igg)$nOverlap[!is.na(ind)] <-
df_enrichment$nOverlap[ind[!is.na(ind)]]
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
else if (node.color == "nOverlap") {
node.color.title <- "# genes"
if (is.null(zlim)) {
zlim <- c(0, ceiling(max(df_enrichment_group$nOverlap)))
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
else if (node.size == "nOverlap") {
node.size.title <- "Num of genes"
if (is.null(slim)) {
slim <- c(0, ceiling(max(df_enrichment_group$nOverlap)))
}
}
x <- y <- leaf <- x_group <- y_group <- group <- group_id <- color <-
size <- NULL
if (length(ls_igg) > 1) {
leave.label.orientation <- "inwards"
}
gp <- xCtree(ig, leave.label.orientation = leave.label.orientation,
...)
df_tips <- subset(gp$data, leaf == T)
ls_df <- lapply(1:length(ls_igg), function(i) {
g <- ls_igg[[i]]
if (length(ls_igg) < 10) {
group <- paste0(i, ": ", names(ls_igg)[i])
group_id <- i
}
else {
group <- sprintf("\%02d: \%s", i, names(ls_igg)[i])
group_id <- sprintf("\%02d", i)
}
x_group <- 1.05 + (i - 1) * group.gap
y_group <- 1.05 + (i - 1) * group.gap
ind <- match(V(g)$name, df_tips$name)
if (node.color == "or") {
color <- V(g)$or[!is.na(ind)]
}
else if (node.color == "adjp") {
color <- V(g)$adjp[!is.na(ind)]
}
else if (node.color == "zscore") {
color <- V(g)$zscore[!is.na(ind)]
}
else if (node.color == "nOverlap") {
color <- V(g)$nOverlap[!is.na(ind)]
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
else if (node.size == "nOverlap") {
size <- V(g)$nOverlap[!is.na(ind)]
}
size[size <= slim[1]] <- slim[1]
size[size >= slim[2]] <- slim[2]
df <- data.frame(df_tips[ind[!is.na(ind)], ], color = color,
size = size, group = group, group_id = group_id,
x_group = x_group, y_group = y_group, stringsAsFactors = F)
})
df <- do.call(rbind, ls_df)
if (length(ls_igg) > 1) {
gp <- gp + geom_polygon(data = df, aes(x = x * x_group,
y = y * y_group, group = group, color = group),
fill = "transparent", size = group.size) + scale_colour_manual(values =
xColormap(paste0(group.color,
"-", group.color))(length(unique(df$group))),
guide = guide_legend(title = "Circles", keywidth = 0.5,
keyheight = 0.6, order = 1, nrow = min(20,
length(unique(df$group)))))
tipid <- NULL
gp <- gp + geom_text(data = subset(df, tipid == 1),
aes(x = x * x_group + 0.1, y = y * y_group, label = group_id),
color = group.label.color, size = group.label.size)
}
gp <- gp + geom_point(data = df, aes(x = x * x_group,
y = y * y_group, fill = color, size = size), alpha = 1,
color = group.color, shape = 21)
if (legend.direction == "auto") {
if (length(ls_igg) > 1) {
legend.direction <- "horizontal"
}
else {
legend.direction <- "vertical"
}
}
if (legend.direction == "vertical") {
gp <- gp + scale_size_continuous(limits = slim, range =
node.size.range,
guide = guide_legend(node.size.title, title.position = "top",
ncol = 1, order = 2))
gp <- gp + scale_fill_gradientn(colors = xColormap(colormap)(64),
limits = zlim, guide = guide_colorbar(title = node.color.title,
title.position = "top", barwidth = 0.5, order = 3))
}
else if (legend.direction == "horizontal") {
gp <- gp + scale_size_continuous(limits = slim, range =
node.size.range,
guide = guide_legend(node.size.title, title.position = "top",
keywidth = 0.5, keyheight = 0.5, ncol = 3,
byrow = T, order = 2))
gp <- gp + scale_fill_gradientn(colors = xColormap(colormap)(64),
limits = zlim, guide = guide_colorbar(title = node.color.title,
title.position = "top", barheight = 0.5, direction = "horizontal",
order = 3))
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
