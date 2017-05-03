\name{xSubneterGenes}
\alias{xSubneterGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSubneterGenes(data, network = c("STRING_highest", "STRING_high",
"STRING_medium", "STRING_low", "PCommonsUN_high", "PCommonsUN_medium",
"PCommonsDN_high", "PCommonsDN_medium", "PCommonsDN_Reactome",
"PCommonsDN_KEGG", "PCommonsDN_HumanCyc", "PCommonsDN_PID",
"PCommonsDN_PANTHER", "PCommonsDN_ReconX", "PCommonsDN_TRANSFAC",
"PCommonsDN_PhosphoSite", "PCommonsDN_CTD"), STRING.only = c(NA,
"neighborhood_score", "fusion_score", "cooccurence_score",
"coexpression_score", "experimental_score", "database_score",
"textmining_score")[1], network.customised = NULL, seed.genes = T,
subnet.significance = 0.01, subnet.size = NULL, verbose = T,
RData.location = "http://galahad.well.ox.ac.uk/bigdata")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{network}{
%% ~~Describe \code{network} here~~
}
\item{STRING.only}{
%% ~~Describe \code{STRING.only} here~~
}
\item{network.customised}{
%% ~~Describe \code{network.customised} here~~
}
\item{seed.genes}{
%% ~~Describe \code{seed.genes} here~~
}
\item{subnet.significance}{
%% ~~Describe \code{subnet.significance} here~~
}
\item{subnet.size}{
%% ~~Describe \code{subnet.size} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
}
\item{RData.location}{
%% ~~Describe \code{RData.location} here~~
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
function (data, network = c("STRING_highest", "STRING_high",
"STRING_medium", "STRING_low", "PCommonsUN_high", "PCommonsUN_medium",
"PCommonsDN_high", "PCommonsDN_medium", "PCommonsDN_Reactome",
"PCommonsDN_KEGG", "PCommonsDN_HumanCyc", "PCommonsDN_PID",
"PCommonsDN_PANTHER", "PCommonsDN_ReconX", "PCommonsDN_TRANSFAC",
"PCommonsDN_PhosphoSite", "PCommonsDN_CTD"), STRING.only = c(NA,
"neighborhood_score", "fusion_score", "cooccurence_score",
"coexpression_score", "experimental_score", "database_score",
"textmining_score")[1], network.customised = NULL, seed.genes = T,
subnet.significance = 0.01, subnet.size = NULL, verbose = T,
RData.location = "http://galahad.well.ox.ac.uk/bigdata")
{
startT <- Sys.time()
if (verbose) {
message(paste(c("Start at ", as.character(startT)), collapse = ""),
appendLF = T)
message("", appendLF = T)
}
network <- match.arg(network)
if (is.null(data)) {
stop("The input data must be not NULL.\n")
}
if (is.vector(data)) {
if (length(data) > 1) {
if (is.null(names(data))) {
stop("The input data must have names with attached gene symbols.\n")
}
}
else {
data <- utils::read.delim(file = data, header = F,
row.names = NULL, stringsAsFactors = F)
}
}
if (is.vector(data)) {
pval <- data
}
else if (is.matrix(data) | is.data.frame(data)) {
data <- as.matrix(data)
data_list <- split(x = data[, 2], f = as.character(data[,
1]))
res_list <- lapply(data_list, function(x) {
x <- as.numeric(x)
x <- x[!is.na(x)]
if (length(x) > 0) {
min(x)
}
else {
NULL
}
})
pval <- unlist(res_list)
}
if (!is.null(network.customised) && class(network.customised) ==
"igraph") {
if (verbose) {
now <- Sys.time()
message(sprintf("Load the customised network (\%s) ...",
as.character(now)), appendLF = T)
}
g <- network.customised
}
else {
if (verbose) {
now <- Sys.time()
message(sprintf("Load the network \%s (\%s) ...", network,
as.character(now)), appendLF = T)
}
g <- xDefineNet(network = network, STRING.only = STRING.only,
weighted = FALSE, verbose = FALSE, RData.location = RData.location)
}
if (verbose) {
message(sprintf("The network you choose has \%d nodes and \%d edges",
vcount(g), ecount(g)), appendLF = T)
}
if (seed.genes) {
ind <- match(V(g)$name, names(pval))
nodes_mapped <- V(g)$name[!is.na(ind)]
g <- dnet::dNetInduce(g = g, nodes_query = nodes_mapped,
knn = 0, remove.loops = F, largest.comp = T)
}
else {
ind <- match(V(g)$name, names(pval))
nodes_not_mapped <- V(g)$name[is.na(ind)]
pval_not_mapped <- rep(1, length(nodes_not_mapped))
names(pval_not_mapped) <- nodes_not_mapped
pval <- c(pval, pval_not_mapped)
}
if (verbose) {
message(sprintf("Restricted to data/nodes of interest, the network
(with the largest interconnected component) has \%d nodes and \%d
edges",
vcount(g), ecount(g)), appendLF = T)
}
if (verbose) {
now <- Sys.time()
message(sprintf("\n#######################################################",
appendLF = T))
message(sprintf("Start to identify a subnetwork (\%s):",
as.character(now)), appendLF = T)
message(sprintf("#######################################################",
appendLF = T))
}
subnet <- dNetPipeline(g = g, pval = pval, method = "customised",
significance.threshold = subnet.significance, nsize = subnet.size,
plot = F, verbose = verbose)
if (ecount(subnet) > 0 && class(subnet) == "igraph") {
relations <- igraph::get.data.frame(subnet, what = "edges")[,
c(1, 2)]
nodes <- igraph::get.data.frame(subnet, what = "vertices")
nodes <- cbind(name = nodes$name, description = nodes$description,
significance = pval[rownames(nodes)], score = nodes$score)
if (is.directed(subnet)) {
subg <- igraph::graph.data.frame(d = relations, directed = TRUE,
vertices = nodes)
}
else {
subg <- igraph::graph.data.frame(d = relations, directed = FALSE,
vertices = nodes)
}
}
else {
subg <- NULL
}
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("The subnetwork has been identified (\%s)!",
as.character(now)), appendLF = T)
message(sprintf("#######################################################\n",
appendLF = T))
}
endT <- Sys.time()
if (verbose) {
message(paste(c("\nFinish at ", as.character(endT)),
collapse = ""), appendLF = T)
}
runTime <- as.numeric(difftime(strptime(endT, "\%Y-\%m-\%d
\%H:\%M:\%S"),
strptime(startT, "\%Y-\%m-\%d \%H:\%M:\%S"), units = "secs"))
message(paste(c("Runtime in total is: ", runTime, " secs\n"),
collapse = ""), appendLF = T)
return(subg)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
