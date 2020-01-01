\name{xDefineNet}
\alias{xDefineNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xDefineNet(network = c("STRING_highest", "STRING_high",
"STRING_medium", "STRING_low", "PCommonsUN_high", "PCommonsUN_medium",
"PCommonsDN_high", "PCommonsDN_medium", "PCommonsDN_Reactome",
"PCommonsDN_KEGG", "PCommonsDN_HumanCyc", "PCommonsDN_PID",
"PCommonsDN_PANTHER", "PCommonsDN_ReconX", "PCommonsDN_TRANSFAC",
"PCommonsDN_PhosphoSite", "PCommonsDN_CTD", "KEGG", "KEGG_metabolism",
"KEGG_genetic", "KEGG_environmental", "KEGG_cellular",
"KEGG_organismal", "KEGG_disease", "REACTOME", "TRRUST"), STRING.only =
c(NA, "neighborhood_score", "fusion_score", "cooccurence_score",
"coexpression_score", "experimental_score", "database_score",
"textmining_score")[1], weighted = FALSE, verbose = TRUE,
RData.location = "http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{network}{
%% ~~Describe \code{network} here~~
}
\item{STRING.only}{
%% ~~Describe \code{STRING.only} here~~
}
\item{weighted}{
%% ~~Describe \code{weighted} here~~
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
function (network = c("STRING_highest", "STRING_high", "STRING_medium",
"STRING_low", "PCommonsUN_high", "PCommonsUN_medium",
"PCommonsDN_high",
"PCommonsDN_medium", "PCommonsDN_Reactome", "PCommonsDN_KEGG",
"PCommonsDN_HumanCyc", "PCommonsDN_PID", "PCommonsDN_PANTHER",
"PCommonsDN_ReconX", "PCommonsDN_TRANSFAC", "PCommonsDN_PhosphoSite",
"PCommonsDN_CTD", "KEGG", "KEGG_metabolism", "KEGG_genetic",
"KEGG_environmental", "KEGG_cellular", "KEGG_organismal",
"KEGG_disease", "REACTOME", "TRRUST"), STRING.only = c(NA,
"neighborhood_score", "fusion_score", "cooccurence_score",
"coexpression_score", "experimental_score", "database_score",
"textmining_score")[1], weighted = FALSE, verbose = TRUE,
RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
network <- match.arg(network)
if (verbose) {
now <- Sys.time()
message(sprintf("Load the network \%s (\%s) ...", network,
as.character(now)), appendLF = TRUE)
}
if (length(grep("STRING", network, perl = TRUE)) > 0) {
g <- xRDataLoader(RData.customised = "org.Hs.string",
RData.location = RData.location, guid = guid, verbose = verbose)
flag <- unlist(strsplit(network, "_"))[2]
if (flag == "highest") {
eval(parse(text = "g <- igraph::subgraph.edges(g,
eids=E(g)[combined_score>=900])"))
}
else if (flag == "high") {
eval(parse(text = "g <- igraph::subgraph.edges(g,
eids=E(g)[combined_score>=700])"))
}
else if (flag == "medium") {
eval(parse(text = "g <- igraph::subgraph.edges(g,
eids=E(g)[combined_score>=400])"))
}
else if (flag == "low") {
eval(parse(text = "g <- igraph::subgraph.edges(g,
eids=E(g)[combined_score>=150])"))
}
default.STRING.only <- c("neighborhood_score", "fusion_score",
"cooccurence_score", "coexpression_score", "experimental_score",
"database_score", "textmining_score")
ind <- match(default.STRING.only, STRING.only)
STRING.only <- default.STRING.only[!is.na(ind)]
if (length(STRING.only) > 0) {
x <- sapply(STRING.only, function(x) paste0(x, ">0"))
x <- paste0(x, collapse = " | ")
x <- paste0("g <- igraph::subgraph.edges(g, eids=E(g)[",
x, "])")
eval(parse(text = x))
}
V(g)$name <- V(g)$symbol
if (weighted) {
relations <- igraph::get.data.frame(g, what = "edges")[,
c(1, 2, 10)]
colnames(relations) <- c("from", "to", "weight")
}
else {
relations <- igraph::get.data.frame(g, what = "edges")[,
c(1, 2)]
colnames(relations) <- c("from", "to")
relations$weight <- rep(1, nrow(relations))
}
if (0) {
ind <- which(!(grepl("HLA-", relations$from) & grepl("HLA-",
relations$to)))
relations <- relations[ind, ]
}
g <- igraph::delete_vertex_attr(g, "name")
g <- igraph::delete_vertex_attr(g, "seqid")
g <- igraph::delete_vertex_attr(g, "geneid")
nodes <- igraph::get.data.frame(g, what = "vertices")
nodes <- nodes[!duplicated(nodes), ]
g <- igraph::graph.data.frame(d = relations, directed = FALSE,
vertices = nodes)
}
else if (length(grep("PCommonsUN", network, perl = TRUE)) >
0) {
g <- xRDataLoader(RData.customised = "org.Hs.PCommons_UN",
RData.location = RData.location, guid = guid, verbose = verbose)
flag <- unlist(strsplit(network, "_"))[2]
if (flag == "high") {
eval(parse(text = "g <- igraph::subgraph.edges(g,
eids=E(g)[in_complex_with>=102 | interacts_with>=102])"))
}
else if (flag == "medium") {
eval(parse(text = "g <- igraph::subgraph.edges(g,
eids=E(g)[in_complex_with>=101 | interacts_with>=101])"))
}
relations <- igraph::get.data.frame(g, what = "edges")[,
c(1, 2)]
colnames(relations) <- c("from", "to")
nodes <- igraph::get.data.frame(g, what = "vertices")[,
c(3, 4)]
g <- igraph::graph.data.frame(d = relations, directed = FALSE,
vertices = nodes)
}
else if (length(grep("PCommonsDN", network, perl = TRUE)) >
0) {
flag <- unlist(strsplit(network, "_"))[2]
if (flag == "high") {
g <- xRDataLoader(RData.customised = "org.Hs.PCommons_DN",
RData.location = RData.location, guid = guid,
verbose = verbose)
eval(parse(text = "g <- igraph::subgraph.edges(g,
eids=E(g)[catalysis_precedes>=102 | controls_expression_of>=102 |
controls_phosphorylation_of>=102 | controls_state_change_of>=102 |
controls_transport_of>=102])"))
}
else if (flag == "medium") {
g <- xRDataLoader(RData.customised = "org.Hs.PCommons_DN",
RData.location = RData.location, guid = guid,
verbose = verbose)
eval(parse(text = "g <- igraph::subgraph.edges(g,
eids=E(g)[catalysis_precedes>=101 | controls_expression_of>=101 |
controls_phosphorylation_of>=101 | controls_state_change_of>=101 |
controls_transport_of>=101])"))
}
else {
g <- xRDataLoader(RData.customised = "org.Hs.PCommons_DN.source",
RData.location = RData.location, guid = guid,
verbose = verbose)
g <- g[[flag]]
eval(parse(text = "g <- igraph::subgraph.edges(g,
eids=E(g)[catalysis_precedes>=101 | controls_expression_of>=101 |
controls_phosphorylation_of>=101 | controls_state_change_of>=101 |
controls_transport_of>=101])"))
}
relations <- igraph::get.data.frame(g, what = "edges")[,
c(1, 2)]
colnames(relations) <- c("from", "to")
relations$weight <- rep(1, nrow(relations))
nodes <- igraph::get.data.frame(g, what = "vertices")[,
c(3, 4)]
g <- igraph::graph.data.frame(d = relations, directed = TRUE,
vertices = nodes)
}
else if (network == "KEGG") {
g <- xRDataLoader(RData.customised = "ig.KEGG.merged",
RData.location = RData.location, guid = guid, verbose = verbose)
g <- igraph::delete_vertex_attr(g, "hsa")
g <- igraph::delete_vertex_attr(g, "GeneID")
g <- igraph::delete_vertex_attr(g, "Symbol")
E(g)$weight <- 1
if (1) {
v <- V(g)[grepl("HLA-", V(g)$name)]
g <- delete_vertices(g, v)
}
}
else if (length(grep("KEGG_", network, perl = TRUE)) > 0) {
ls_ig <- xRDataLoader(RData.customised = "ig.KEGG.mergedCategory",
RData.location = RData.location, guid = guid, verbose = verbose)
if (network == "KEGG_metabolism") {
g <- ls_ig[["Metabolism"]]
}
else if (network == "KEGG_genetic") {
g <- ls_ig[["Genetic Information Processing"]]
}
else if (network == "KEGG_environmental") {
g <- ls_ig[["Environmental Information Processing"]]
}
else if (network == "KEGG_cellular") {
g <- ls_ig[["Cellular Processes"]]
}
else if (network == "KEGG_organismal") {
g <- ls_ig[["Organismal Systems"]]
}
else if (network == "KEGG_disease") {
g <- ls_ig[["Human Diseases"]]
}
g <- igraph::delete_vertex_attr(g, "hsa")
g <- igraph::delete_vertex_attr(g, "GeneID")
g <- igraph::delete_vertex_attr(g, "Symbol")
E(g)$weight <- 1
if (1) {
v <- V(g)[grepl("HLA-", V(g)$name)]
g <- delete_vertices(g, v)
}
}
else if (network == "REACTOME") {
g <- xRDataLoader(RData.customised = "ig.REACTOME.merged",
RData.location = RData.location, guid = guid, verbose = verbose)
g <- igraph::delete_vertex_attr(g, "geneid")
g <- igraph::delete_vertex_attr(g, "symbol")
E(g)$weight <- 1
}
else if (network == "TRRUST") {
g <- xRDataLoader(RData.customised = "ig.TRRUST", RData.location =
RData.location,
guid = guid, verbose = verbose)
E(g)$weight <- 1
}
invisible(g)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
