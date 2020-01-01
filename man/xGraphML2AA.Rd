\name{xGraphML2AA}
\alias{xGraphML2AA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGraphML2AA(data = NULL, org = c("human", "mouse"), query =
"AA:hsa04672", curation = c("manual", "automatic", "any"), node.label =
"label", node.color = "lfc", colormap =
"deepskyblue-lightyellow-darkorange", ncolors = 64, nlegend = 9, zlim =
NULL, legend.title = "", title.thispath = NULL, node.tooltip =
"tooltip", node.highlight = "fdr", node.highlight.cutoff = 0.05,
edge.color = "#00000033", edge.width = 1, color.gene = "#dddddd",
color.thispath = "#dddddd", color.otherpath = "#eeeeee", size.gene =
10, size.gene.found = 11, size.gene.highlight = 12, filename =
"xGraphML2AA", verbose = TRUE, RData.location =
"http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{org}{
%% ~~Describe \code{org} here~~
}
\item{query}{
%% ~~Describe \code{query} here~~
}
\item{curation}{
%% ~~Describe \code{curation} here~~
}
\item{node.label}{
%% ~~Describe \code{node.label} here~~
}
\item{node.color}{
%% ~~Describe \code{node.color} here~~
}
\item{colormap}{
%% ~~Describe \code{colormap} here~~
}
\item{ncolors}{
%% ~~Describe \code{ncolors} here~~
}
\item{nlegend}{
%% ~~Describe \code{nlegend} here~~
}
\item{zlim}{
%% ~~Describe \code{zlim} here~~
}
\item{legend.title}{
%% ~~Describe \code{legend.title} here~~
}
\item{title.thispath}{
%% ~~Describe \code{title.thispath} here~~
}
\item{node.tooltip}{
%% ~~Describe \code{node.tooltip} here~~
}
\item{node.highlight}{
%% ~~Describe \code{node.highlight} here~~
}
\item{node.highlight.cutoff}{
%% ~~Describe \code{node.highlight.cutoff} here~~
}
\item{edge.color}{
%% ~~Describe \code{edge.color} here~~
}
\item{edge.width}{
%% ~~Describe \code{edge.width} here~~
}
\item{color.gene}{
%% ~~Describe \code{color.gene} here~~
}
\item{color.thispath}{
%% ~~Describe \code{color.thispath} here~~
}
\item{color.otherpath}{
%% ~~Describe \code{color.otherpath} here~~
}
\item{size.gene}{
%% ~~Describe \code{size.gene} here~~
}
\item{size.gene.found}{
%% ~~Describe \code{size.gene.found} here~~
}
\item{size.gene.highlight}{
%% ~~Describe \code{size.gene.highlight} here~~
}
\item{filename}{
%% ~~Describe \code{filename} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
}
\item{RData.location}{
%% ~~Describe \code{RData.location} here~~
}
\item{guid}{
%% ~~Describe \code{guid} here~~
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
function (data = NULL, org = c("human", "mouse"), query =
"AA:hsa04672",
curation = c("manual", "automatic", "any"), node.label = "label",
node.color = "lfc", colormap = "deepskyblue-lightyellow-darkorange",
ncolors = 64, nlegend = 9, zlim = NULL, legend.title = "",
title.thispath = NULL, node.tooltip = "tooltip", node.highlight =
"fdr",
node.highlight.cutoff = 0.05, edge.color = "#00000033", edge.width = 1,
color.gene = "#dddddd", color.thispath = "#dddddd", color.otherpath =
"#eeeeee",
size.gene = 10, size.gene.found = 11, size.gene.highlight = 12,
filename = "xGraphML2AA", verbose = TRUE, RData.location =
"http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
startT <- Sys.time()
if (verbose) {
message(paste(c("Start at ", as.character(startT)), collapse = ""),
appendLF = TRUE)
message("", appendLF = TRUE)
}
curation <- match.arg(curation)
if (any(class(data) \%in\% c("tbl_df", "tbl"))) {
data <- as.data.frame(data)
}
if (all(c(node.label, node.color) \%in\% colnames(data))) {
df <- data.frame(Symbol = data[, node.label], LFC = signif(data[,
node.color], digits = 2), stringsAsFactors = FALSE)
}
else {
df <- data.frame(Symbol = NA, LFC = 0, stringsAsFactors = FALSE)
data <- df
}
df$FDR <- rep(1, nrow(df))
if (!is.null(node.highlight)) {
if (node.highlight \%in\% colnames(data)) {
df$FDR <- signif(data[, node.highlight], digits = 2)
}
}
df <- df[!is.na(as.numeric(df$LFC)), ]
df <- df[!is.na(as.numeric(df$FDR)), ]
df_res <- xSymbol2GeneID(df$Symbol, org = org, details = TRUE,
verbose = verbose, RData.location = RData.location, guid = guid)
df$GeneID <- df_res$GeneID
df$description <- df_res$description
if (node.highlight \%in\% colnames(data)) {
df$tooltip <- paste0("Symbol: ", df$Symbol, "\nName: ",
df$description, "\nColor: ", df$LFC, "\nHightlight: ",
df$FDR)
}
else {
df$tooltip <- paste0("Symbol: ", df$Symbol, "\nName: ",
df$description, "\nColor: ", df$LFC)
}
if (!is.null(node.tooltip)) {
if (node.tooltip \%in\% colnames(data)) {
df$tooltip <- data[, node.tooltip]
}
}
df_legends <- NULL
if (!is.null(node.color)) {
ind <- match(node.color, colnames(data))
if (!is.na(ind)) {
pattern <- as.numeric(df$LFC)
pattern_nona <- pattern[!is.na(pattern)]
pattern_nona <- as.numeric(pattern_nona)
if (is.null(zlim)) {
vmin <- floor(stats::quantile(pattern_nona, 0.05))
vmax <- ceiling(stats::quantile(pattern_nona,
0.95))
if (vmin < 0 & vmax > 0) {
vsym <- abs(min(vmin, vmax))
vmin <- -1 * vsym
vmax <- vsym
}
zlim <- c(vmin, vmax)
}
vec2color <- function(vec, colormap = colormap, ncolors = ncolors,
zlim = zlim) {
palette.name <- xColormap(colormap = colormap)
colors <- palette.name(ncolors)
if (max(zlim) != min(zlim)) {
scale <- length(colors)/(max(zlim) - min(zlim))
sapply(1:length(vec), function(x) {
if (is.na(vec[x])) {
"transparent"
}
else {
ind <- floor(1 + (vec[x] - min(zlim)) *
scale)
colors[max(1, min(ncolors, ind))]
}
})
}
else {
colors[rep(ncolors, length(vec))]
}
}
node.color <- vec2color(pattern, colormap = colormap,
ncolors = ncolors, zlim = zlim)
if (length(unique(node.color)) > 1) {
colors <- xColormap(colormap = colormap)(ncolors)
legend_colors <- colors[round(seq(1, ncolors,
length.out = nlegend))]
df_legends <- data.frame(name = paste0("Legend",
(nlegend - 1):0), colors = legend_colors, labels =
signif(seq(min(zlim),
max(zlim), length.out = nlegend), digits = 2),
stringsAsFactors = FALSE)
}
}
else {
node.color <- rep("#BFFFBF", nrow(data))
}
}
else {
node.color <- rep("#BFFFBF", nrow(data))
}
df$node.color <- node.color
manual_ind <- NULL
manual_ind_at <- NULL
if (curation \%in\% c("any", "manual")) {
AA.template <- xRDataLoader("AA.template", verbose = verbose,
RData.location = RData.location, guid = guid)
info <- AA.template$info
path <- gsub("^AA:", "", info$id)
query <- gsub("^AA:", "", query)
query <- gsub("^path:", "", query)
manual_ind <- match(query, path)
if (is.na(manual_ind)) {
manual_ind <- grep(query, info$name)
if (length(manual_ind) > 1) {
warning(sprintf("Manual curation: \%d found for queried \%s: only 1st
kept",
length(manual_ind), query), appendLF = TRUE)
manual_ind <- manual_ind[1]
}
else if (length(manual_ind) == 0) {
warning(sprintf("Manual curation: no found for queried '\%s'",
query), appendLF = TRUE)
AT.path <- xRDataLoader("AT.path", verbose = verbose,
RData.location = RData.location, guid = guid)
info <- AT.path$info
path <- gsub("^AT:", "", info$id)
query <- gsub("^AT:", "", query)
query <- gsub("^path:", "", query)
manual_ind_at <- match(query, path)
if (is.na(manual_ind_at)) {
manual_ind_at <- grep(query, info$name)
if (length(manual_ind_at) > 1) {
warning(sprintf("Automatic curation: \%d found for queried \%s: only
1st kept",
length(manual_ind_at), query), appendLF = TRUE)
manual_ind_at <- manual_ind_at[1]
}
else if (length(manual_ind_at) == 0) {
warning(sprintf("Automatic curation: no found for queried '\%s'",
query), appendLF = TRUE)
if (curation == "manual") {
return(NULL)
}
}
}
}
}
}
if (length(manual_ind) != 0 | length(manual_ind_at) != 0) {
if (length(manual_ind) != 0) {
if (verbose) {
now <- Sys.time()
message(sprintf("Manual curation: visualising '\%s: \%s' (\%s) ...",
AA.template$info$id[manual_ind], AA.template$info$name[manual_ind],
as.character(now)), appendLF = TRUE)
}
detail <- AA.template$detail[[manual_ind]]
df_nodes <- detail$nodes
df_edges <- detail$edges
}
else {
if (verbose) {
now <- Sys.time()
message(sprintf("Automatic curation: visualising '\%s: \%s' (\%s) ...",
AT.path$info$id[manual_ind_at], AT.path$info$name[manual_ind_at],
as.character(now)), appendLF = TRUE)
}
detail <- AT.path$detail[[manual_ind_at]]
df_nodes <- detail$nodes
df_nodes$node_id <- gsub("::.*", "", df_nodes$node_id)
df_edges <- detail$edges
df_edges$source <- gsub("::.*", "", df_edges$source)
df_edges$target <- gsub("::.*", "", df_edges$target)
}
if (org[1] == "mouse") {
Human2Mouse <- xRDataLoader("Human2Mouse", RData.location =
RData.location,
guid = guid)
ind <- match(df_nodes$name, Human2Mouse$Human)
df_nodes$name[!is.na(ind)] <- Human2Mouse$Mouse[ind[!is.na(ind)]]
}
output.head <- "<?xml version=\"1.0\" encoding=\"UTF-8\"
standalone=\"no\"?>\n\t<graphml
xmlns=\"http://graphml.graphdrawing.org/xmlns\"
xmlns:java=\"http://www.yworks.com/xml/yfiles-common/1.0/java\"
xmlns:sys=\"http://www.yworks.com/xml/yfiles-common/markup/primitives/2.0\"
xmlns:x=\"http://www.yworks.com/xml/yfiles-common/markup/2.0\"
xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"
xmlns:y=\"http://www.yworks.com/xml/graphml\"
xmlns:yed=\"http://www.yworks.com/xml/yed/3\"
xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns
http://www.yworks.com/xml/schema/graphml/1.1/ygraphml.xsd\">\n\t
<!--Created by yEd 3.17-->\n\t <key for=\"node\" id=\"d1\"
attr.name=\"url\" attr.type=\"string\"/>\n\t <key for=\"node\"
id=\"d2\" attr.name=\"description\" attr.type=\"string\"/>\n\t <key
for=\"node\" id=\"d3\" yfiles.type=\"nodegraphics\"/>\n\t <key
for=\"edge\" id=\"d4\" attr.name=\"description\"
attr.type=\"string\"/>\n\t <key for=\"edge\" id=\"d5\"
yfiles.type=\"edgegraphics\"/>\n\t <key for=\"graphml\" id=\"d6\"
yfiles.type=\"resources\"/>\n\t <graph edgedefault=\"directed\"
id=\"G\">"
n_total <- 0
n_matched <- 0
for (i in 1:nrow(df_nodes)) {
if (grepl("Legend", df_nodes$name[i])) {
if (is.null(df_legends)) {
next
}
}
ind_found <- match(df_nodes$name[i], df$Symbol)
ind_legend <- match(df_nodes$name[i], df_legends$name)
flag_gene <- grepl("#BFFFBF", df_nodes$fill[i])
flag_thispath <- grepl("#00FF00", df_nodes$fill[i])
flag_otherpath <- grepl("#C0C0C0", df_nodes$fill[i])
if (flag_gene & is.na(ind_legend)) {
n_total <- n_total + 1
if (!is.na(ind_found)) {
n_matched <- n_matched + 1
}
}
}
if (verbose) {
now <- Sys.time()
message(sprintf("\t \%d genes in total and \%d genes matched (\%s)",
n_total, n_matched, as.character(now)), appendLF = TRUE)
}
ls_nodes <- lapply(1:nrow(df_nodes), function(i) {
if (grepl("Legend", df_nodes$name[i])) {
if (is.null(df_legends)) {
return(NULL)
}
}
ind_found <- match(df_nodes$name[i], df$Symbol)
ind_legend <- match(df_nodes$name[i], df_legends$name)
flag_gene <- grepl("#BFFFBF", df_nodes$fill[i])
flag_thispath <- grepl("#00FF00", df_nodes$fill[i])
flag_otherpath <- grepl("#C0C0C0", df_nodes$fill[i])
k <- 0
vec <- vector()
k <- k + 1
vec[k] <- paste0("<node id=\"", df_nodes$node_id[i],
"\">")
if (!is.na(ind_found)) {
k <- k + 1
vec[k] <- paste0("<data key=\"d1\"><![CDATA[",
"https://www.ncbi.nlm.nih.gov/gene/", df$GeneID[ind_found],
"]]></data>")
k <- k + 1
vec[k] <- paste0("<data key=\"d2\"><![CDATA[",
df$tooltip[ind_found], "]]></data>")
}
else {
if (!is.na(ind_legend)) {
k <- k + 1
vec[k] <- paste0("<data key=\"d2\"><![CDATA[",
df_legends$labels[ind_legend], "]]></data>")
}
else {
if (grepl("Legend", df_nodes$name[i])) {
return(NULL)
}
else if (flag_gene) {
k <- k + 1
vec[k] <- paste0("<data key=\"d1\"><![CDATA[",
"http://www.genecards.org/cgi-bin/carddisp.pl?gene=",
df_nodes$name[i], "]]></data>")
k <- k + 1
vec[k] <- paste0("<data key=\"d2\"><![CDATA[",
df_nodes$name[i], "]]></data>")
}
}
}
k <- k + 1
vec[k] <- paste0("<data key=\"d3\">")
k <- k + 1
vec[k] <- paste0("<y:ShapeNode>")
k <- k + 1
vec[k] <- paste0("<y:Geometry ", df_nodes$geometry[i],
"/>")
fill <- df_nodes$fill[i]
if (!is.na(ind_found)) {
fill <- gsub("#BFFFBF", df$node.color[ind_found],
fill)
}
else {
if (!is.na(ind_legend)) {
fill <- gsub("#BFFFBF", df_legends$colors[ind_legend],
fill)
}
else {
fill <- gsub("#BFFFBF", color.gene, fill)
fill <- gsub("#00FF00", color.thispath, fill)
fill <- gsub("#C0C0C0", color.otherpath, fill)
}
}
k <- k + 1
vec[k] <- paste0("<y:Fill ", fill, "/>")
borderstyle <- df_nodes$borderstyle[i]
if (!is.na(ind_legend)) {
borderstyle <- gsub("color=\"#\\w+\"", "color=\"#cccccc\"",
borderstyle)
borderstyle <- gsub("width=\".*\"", "width=\"0.8\"",
borderstyle)
}
if (!is.na(ind_found)) {
borderstyle <- gsub("color=\"#\\w+\"", "color=\"#ffffff00\"",
borderstyle)
borderstyle <- gsub("width=\".*\"", "width=\"0\"",
borderstyle)
}
else if (flag_gene) {
borderstyle <- gsub("color=\"#\\w+\"", "color=\"#ffffff00\"",
borderstyle)
borderstyle <- gsub("width=\".*\"", "width=\"0\"",
borderstyle)
}
k <- k + 1
vec[k] <- paste0("<y:BorderStyle ", borderstyle,
"/>")
nodelabel <- df_nodes$nodelabel[i]
nodelabel <- gsub("fontFamily=\"\\w+\"", "fontFamily=\"Arial\"",
nodelabel)
if (!is.na(ind_found)) {
if (df$FDR[ind_found] < node.highlight.cutoff) {
nodelabel <- gsub("fontStyle=\"\\w+\"", "fontStyle=\"bolditalic\"",
nodelabel)
tmp <- paste0("fontSize=\"", size.gene.highlight,
"\"")
nodelabel <- gsub("fontSize=\"\\w+\"", tmp,
nodelabel)
}
else {
nodelabel <- gsub("fontStyle=\"\\w+\"", "fontStyle=\"italic\"",
nodelabel)
tmp <- paste0("fontSize=\"", size.gene.found,
"\"")
nodelabel <- gsub("fontSize=\"\\w+\"", tmp,
nodelabel)
}
}
else if (flag_gene & is.na(ind_legend)) {
nodelabel <- gsub("fontStyle=\"\\w+\"", "fontStyle=\"italic\"",
nodelabel)
tmp <- paste0("fontSize=\"", size.gene, "\"")
nodelabel <- gsub("fontSize=\"\\w+\"", tmp, nodelabel)
}
name_display <- df_nodes$name[i]
name_display <- gsub("@@", "\n", name_display)
if (flag_thispath) {
if (!is.null(title.thispath)) {
name_display <- paste0(name_display, " ", title.thispath)
}
else {
name_display <- paste0(name_display, " (AA)")
}
}
if (!is.na(ind_legend)) {
name_display <- df_legends$labels[ind_legend]
}
else {
if (name_display \%in\% c("Pi rating", "5-star\nrating")) {
name_display <- legend.title
nodelabel <- gsub("modelPosition=\"\\w+\"",
"modelPosition=\"t\"", nodelabel)
}
}
k <- k + 1
vec[k] <- paste0("<y:NodeLabel ", nodelabel, ">",
name_display, "</y:NodeLabel>")
k <- k + 1
vec[k] <- paste0("<y:Shape ", df_nodes$shape[i],
"/>")
k <- k + 1
vec[k] <- paste0("</y:ShapeNode>")
k <- k + 1
vec[k] <- paste0("</data>")
k <- k + 1
vec[k] <- paste0("</node>")
paste(vec, collapse = "\n")
})
vec_nodes <- unlist(ls_nodes)
output.nodes <- paste(vec_nodes, collapse = "\n")
ls_edges <- lapply(1:nrow(df_edges), function(i) {
k <- 0
vec <- vector()
k <- k + 1
vec[k] <- paste0("<edge id=\"", df_edges$edge_id[i],
"\" source=\"", df_edges$source[i], "\" target=\"",
df_edges$target[i], "\">")
k <- k + 1
vec[k] <- paste0("<data key=\"d5\">")
k <- k + 1
vec[k] <- paste0("<y:GenericEdge configuration=\"DEFAULT\">")
k <- k + 1
vec[k] <- paste0("<y:Path ", df_edges$path[i], "/>")
linestyle <- df_edges$linestyle[i]
linestyle <- gsub("#000000", "#888888", linestyle)
k <- k + 1
vec[k] <- paste0("<y:LineStyle ", linestyle, "/>")
k <- k + 1
vec[k] <- paste0("<y:Arrows ", df_edges$arrows[i],
"/>")
k <- k + 1
vec[k] <- paste0("</y:GenericEdge>")
k <- k + 1
vec[k] <- paste0("</data>")
k <- k + 1
vec[k] <- paste0("</edge>")
paste(vec, collapse = "\n")
})
vec_edges <- unlist(ls_edges)
output.edges <- paste(vec_edges, collapse = "\n")
output.tail <- "</graph>\n\t <data key=\"d6\">\n\t\t<y:Resources/>\n\t
</data>\n\t</graphml>"
output <- paste0(output.head, "\n", output.nodes, "\n",
output.edges, "\n", output.tail, "\n")
if (!is.null(filename)) {
filename <- gsub(".graphml$", "", filename)
outputfile <- paste0(filename, ".graphml")
fileConn <- base::file(outputfile)
base::writeLines(output, fileConn)
base::close(fileConn)
if (verbose) {
message(sprintf("Congratulations! A file '\%s' (in the directory \%s)
has been created!",
outputfile, getwd()), appendLF = T)
}
}
}
else if (curation \%in\% c("automatic", "any")) {
ls_ig <- xRDataLoader("ig.KEGG.list", verbose = verbose,
RData.location = RData.location, guid = guid)
kegg <- gsub("^path:", "", sapply(ls_ig, function(x) x$path))
query <- gsub("^AA:", "", query)
query <- gsub("^AT:", "", query)
query <- gsub("^path:", "", query)
automatic_ind <- match(query, kegg)
if (is.na(automatic_ind)) {
automatic_ind <- grep(query, names(kegg))
if (length(automatic_ind) == 0) {
warning(sprintf("Automated: no found for queried '\%s'",
query), appendLF = TRUE)
return(NULL)
}
else if (length(automatic_ind) > 1) {
warning(sprintf("Automated: \%d found for queried \%s: only 1st kept",
length(automatic_ind), query), appendLF = TRUE)
automatic_ind <- automatic_ind[1]
}
}
if (verbose) {
now <- Sys.time()
message(sprintf("Automated: visualising '\%s: \%s' (\%s) ...",
kegg[automatic_ind], names(kegg[automatic_ind]),
as.character(now)), appendLF = TRUE)
}
ig <- ls_ig[[automatic_ind]]
if (0) {
ind <- match(V(ig)$name, df$Symbol)
nodes_query <- V(ig)$name[!is.na(ind)]
order <- vcount(ig)
order <- 2
neighs.out <- igraph::neighborhood(ig, order = 2,
nodes = nodes_query, mode = "in")
neighbors <- unique(names(unlist(neighs.out)))
subg <- dnet::dNetInduce(ig, nodes_query = neighbors,
knn = 0, remove.loops = TRUE, largest.comp = FALSE)
}
else {
subg <- ig
}
if (verbose) {
now <- Sys.time()
message(sprintf("\t \%d genes and \%d edges (\%s)",
vcount(subg), ecount(subg), as.character(now)),
appendLF = TRUE)
}
glayout <- igraph::layout_with_kk(subg)
V(subg)$xcoord <- glayout[, 1]
V(subg)$ycoord <- glayout[, 2]
ind <- match(V(subg)$name, df$Symbol)
V(subg)$LFC <- df$LFC[ind]
V(subg)$tooltip <- df$tooltip[ind]
V(subg)$tooltip[is.na(V(subg)$tooltip)] <-
V(subg)$description[is.na(V(subg)$tooltip)]
V(subg)$FDR <- df$FDR[ind]
V(subg)$node.size <- ifelse(V(subg)$FDR < node.highlight.cutoff &
!is.na(V(subg)$FDR), 25, 15)
output <- xGraphML(g = subg, node.label = "name", node.label.size = 10,
node.tooltip = "tooltip", node.xcoord = "xcoord",
node.ycoord = "ycoord", node.color.na = color.gene,
node.color = "LFC", node.link =
"http://www.genecards.org/cgi-bin/carddisp.pl?gene=",
nlegend = nlegend, node.size = "node.size", node.coord.scale = 300,
zlim = zlim, colormap = colormap, filename = filename)
}
endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, "\%Y-\%m-\%d
\%H:\%M:\%S"),
strptime(startT, "\%Y-\%m-\%d \%H:\%M:\%S"), units = "secs"))
if (verbose) {
message(paste(c("\nEnd at ", as.character(endT)), collapse = ""),
appendLF = TRUE)
message(paste(c("Runtime in total is: ", runTime, " secs\n"),
collapse = ""), appendLF = TRUE)
}
invisible(output)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
