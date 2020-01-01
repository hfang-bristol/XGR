\name{xCombineNet}
\alias{xCombineNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xCombineNet(list_ig, combineBy = c("union", "intersect"), attrBy =
c("intersect", "union"), keep.all.vertices = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{list_ig}{
%% ~~Describe \code{list_ig} here~~
}
\item{combineBy}{
%% ~~Describe \code{combineBy} here~~
}
\item{attrBy}{
%% ~~Describe \code{attrBy} here~~
}
\item{keep.all.vertices}{
%% ~~Describe \code{keep.all.vertices} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
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
function (list_ig, combineBy = c("union", "intersect"), attrBy =
c("intersect",
"union"), keep.all.vertices = FALSE, verbose = TRUE)
{
combineBy <- match.arg(combineBy)
attrBy <- match.arg(attrBy)
if (any(class(list_ig) \%in\% c("igraph"))) {
ls_ig <- list(list_ig)
}
else if (class(list_ig) == "list") {
ls_ig <- base::Filter(base::Negate(is.null), list_ig)
if (length(ls_ig) == 0) {
return(NULL)
}
}
else {
stop("The function must apply to 'list' of 'igraph' objects or a
'igraph' object.\n")
}
flag_direct <- "undirect"
vec <- sapply(ls_ig, igraph::is_directed)
if (sum(vec) == length(vec)) {
flag_direct <- "direct"
}
ls_node_attr <- lapply(1:length(ls_ig), function(i) {
ig <- ls_ig[[i]]
igraph::vertex_attr_names(ig)
})
if (attrBy == "intersect") {
node_attr <- base::Reduce(intersect, ls_node_attr)
}
else if (attrBy == "union") {
node_attr <- base::Reduce(union, ls_node_attr)
}
ls_node <- lapply(1:length(ls_ig), function(i) {
ig <- ls_ig[[i]]
nodes <- igraph::get.data.frame(ig, what = "vertices")[,
node_attr]
})
df_node <- unique(do.call(rbind, ls_node))
if (any(duplicated(df_node$name))) {
name <- NULL
df_node <- df_node[!duplicated(df_node$name), ] \%>\%
dplyr::select(name)
node_attr <- "name"
}
ls_edge <- lapply(1:length(ls_ig), function(i) {
ig <- ls_ig[[i]]
relations <- igraph::get.data.frame(ig, what = "edges")[,
c(1, 2)]
if (flag_direct == "undirect") {
ind <- which(relations[, 1] > relations[, 2])
relations[ind, c(1:2)] <- relations[ind, c(2, 1)]
relations <- unique(relations)
}
return(relations)
})
df_edge <- do.call(rbind, ls_edge)
if (combineBy == "intersect") {
from <- to <- count <- NULL
ftc <- df_edge \%>\% dplyr::group_by(from, to) \%>\%
dplyr::summarize(count = n()) \%>\%
dplyr::filter(count == length(ls_ig)) \%>\% dplyr::select(from,
to)
df_edge <- as.data.frame(ftc)
if (!keep.all.vertices) {
tmp_nodes <- unique(c(df_edge[, 1], df_edge[, 2]))
ind <- match(df_node$name, tmp_nodes)
df_node <- df_node[!is.na(ind), ]
}
}
else if (combineBy == "union") {
df_edge <- unique(df_edge)
}
if (flag_direct == "direct") {
combined_ig <- igraph::graph.data.frame(d = df_edge,
directed = TRUE, vertices = df_node)
}
else {
combined_ig <- igraph::graph.data.frame(d = df_edge,
directed = FALSE, vertices = df_node)
}
if (verbose) {
message(sprintf("\%d network(s) are combined (via '\%s') into a '\%s'
network (\%d nodes and \%d edges) with \%d node attributes (via
'\%s')",
length(ls_ig), combineBy, flag_direct, vcount(combined_ig),
ecount(combined_ig), length(node_attr), attrBy, appendLF = TRUE))
}
invisible(combined_ig)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
