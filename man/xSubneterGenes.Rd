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
"PCommonsDN_PhosphoSite", "PCommonsDN_CTD", "KEGG", "KEGG_metabolism",
"KEGG_genetic", "KEGG_environmental", "KEGG_cellular",
"KEGG_organismal", "KEGG_disease", "REACTOME"), STRING.only = c(NA,
"neighborhood_score", "fusion_score", "cooccurence_score",
"coexpression_score", "experimental_score", "database_score",
"textmining_score")[1], network.customised = NULL, seed.genes = T,
subnet.significance = 0.01, subnet.size = NULL, test.permutation = F,
num.permutation = 100, respect = c("none", "degree"), aggregateBy =
c("Ztransform", "fishers", "logistic", "orderStatistic"), verbose = T,
silent = F, RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
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
\item{test.permutation}{
%% ~~Describe \code{test.permutation} here~~
}
\item{num.permutation}{
%% ~~Describe \code{num.permutation} here~~
}
\item{respect}{
%% ~~Describe \code{respect} here~~
}
\item{aggregateBy}{
%% ~~Describe \code{aggregateBy} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
}
\item{silent}{
%% ~~Describe \code{silent} here~~
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
function (data, network = c("STRING_highest", "STRING_high",
"STRING_medium", "STRING_low", "PCommonsUN_high", "PCommonsUN_medium",
"PCommonsDN_high", "PCommonsDN_medium", "PCommonsDN_Reactome",
"PCommonsDN_KEGG", "PCommonsDN_HumanCyc", "PCommonsDN_PID",
"PCommonsDN_PANTHER", "PCommonsDN_ReconX", "PCommonsDN_TRANSFAC",
"PCommonsDN_PhosphoSite", "PCommonsDN_CTD", "KEGG", "KEGG_metabolism",
"KEGG_genetic", "KEGG_environmental", "KEGG_cellular",
"KEGG_organismal",
"KEGG_disease", "REACTOME"), STRING.only = c(NA, "neighborhood_score",
"fusion_score", "cooccurence_score", "coexpression_score",
"experimental_score", "database_score", "textmining_score")[1],
network.customised = NULL, seed.genes = T, subnet.significance = 0.01,
subnet.size = NULL, test.permutation = F, num.permutation = 100,
respect = c("none", "degree"), aggregateBy = c("Ztransform",
"fishers", "logistic", "orderStatistic"), verbose = T,
silent = F, RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
startT <- Sys.time()
if (!silent) {
message(paste(c("Start at ", as.character(startT)), collapse = ""),
appendLF = TRUE)
message("", appendLF = TRUE)
}
else {
verbose <- FALSE
}
network <- match.arg(network)
respect <- match.arg(respect)
aggregateBy <- match.arg(aggregateBy)
if (is.null(data)) {
return(NULL)
}
if (is.vector(data)) {
if (length(data) > 1) {
if (is.null(names(data))) {
return(NULL)
}
}
else {
data <- utils::read.delim(file = data, header = F,
row.names = NULL, stringsAsFactors = F)
}
}
if (is.vector(data)) {
pval <- data[!is.na(data)]
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
weighted = FALSE, verbose = FALSE, RData.location = RData.location,
guid = guid)
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
if (class(suppressWarnings(try(subnet <- dnet::dNetPipeline(g = g,
pval = pval, method = "customised", significance.threshold =
subnet.significance,
nsize = subnet.size, plot = F, verbose = verbose), T))) ==
"try-error") {
subnet <- NULL
}
else {
if (test.permutation & !is.null(subnet.size)) {
if (verbose) {
message(sprintf("Estimate the significance of the identified network
(\%d nodes) based on \%d permutation test respecting '\%s' (\%s) ...",
vcount(subnet), num.permutation, respect, as.character(Sys.time())),
appendLF = T)
}
if (respect == "degree") {
nbin <- 10
vec_degree <- degree(g)
breaks <- unique(stats::quantile(vec_degree,
seq(0, 1, 1/nbin)))
cut_index <- as.numeric(cut(vec_degree, breaks = breaks))
cut_index[is.na(cut_index)] <- 1
names(cut_index) <- V(g)$name
ind <- match(names(pval), names(cut_index))
pval <- pval[!is.na(ind)]
pval_degree <- cut_index[ind[!is.na(ind)]]
B <- num.permutation
set.seed(825)
ls_df <- lapply(1:length(pval), function(i) {
x <- pval[i]
ind <- match(names(x), names(pval_degree))
all_to_sample <- which(pval_degree == pval_degree[ind])
ind_sampled <- base::sample(all_to_sample,
B, replace = T)
res <- data.frame(name = names(pval[i]), pval = pval[ind_sampled],
B = 1:B, stringsAsFactors = F)
})
df_ind_B <- do.call(rbind, ls_df)
ls_index <- split(x = df_ind_B[, c("name", "pval")],
f = df_ind_B$B)
ls_subnet_permutated <- lapply(1:length(ls_index),
function(j) {
if (verbose & j\%\%10 == 0) {
message(sprintf("\t\%d (out of \%d) (\%s) ...",
j, B, as.character(Sys.time())), appendLF = T)
}
pval_permutated <- ls_index[[j]]$pval
names(pval_permutated) <- ls_index[[j]]$name
if (class(suppressWarnings(try(subnet_permutated <-
suppressMessages(dnet::dNetPipeline(g = g,
pval = pval_permutated, method = "customised",
significance.threshold = subnet.significance,
nsize = subnet.size, plot = F, verbose = F)),
T))) == "try-error") {
return(NULL)
}
else {
return(subnet_permutated)
}
})
ls_subnet_permutated <- base::Filter(base::Negate(is.null),
ls_subnet_permutated)
}
else {
B <- num.permutation
ls_subnet_permutated <- list()
set.seed(825)
for (j in 1:B) {
if (verbose & j\%\%10 != 0) {
message(sprintf("\t\%d (out of \%d) (\%s) ...",
j, B, as.character(Sys.time())), appendLF = T)
}
ind <- base::sample(1:length(pval), replace = T)
pval_permutated <- pval[ind]
names(pval_permutated) <- names(pval)
if (class(suppressWarnings(try(ls_subnet_permutated[[j]] <-
suppressMessages(dnet::dNetPipeline(g = g,
pval = pval_permutated, method = "customised",
significance.threshold = subnet.significance,
nsize = subnet.size, plot = F, verbose = F)),
T))) == "try-error") {
ls_subnet_permutated[[j]] <- NULL
}
}
ls_subnet_permutated <- base::Filter(base::Negate(is.null),
ls_subnet_permutated)
}
subnet <- dNetConfidence(target = subnet, sources =
ls_subnet_permutated,
plot = F)
E(subnet)$edgeConfidence <- E(subnet)$edgeConfidence/100
p_combined <- dnet::dPvalAggregate(pmatrix =
matrix(E(subnet)$edgeConfidence,
nrow = 1), method = aggregateBy)
subnet$combinedP <- p_combined
if (verbose) {
message(sprintf("\t'\%s' combined p-value (\%1.2e) of the identified
network (\%d nodes) based on \%d permutation test respecting '\%s'
(\%s)",
aggregateBy, subnet$combinedP, vcount(subnet),
num.permutation, respect, as.character(Sys.time())),
appendLF = T)
}
}
}
if (class(subnet) == "igraph") {
relations <- igraph::get.data.frame(subnet, what = "edges")[,
c(1, 2)]
if (!is.null(subnet$combinedP)) {
relations$edgeConfidence <- igraph::get.data.frame(subnet,
what = "edges")[, "edgeConfidence"]
}
nodes <- igraph::get.data.frame(subnet, what = "vertices")
nodes <- cbind(name = nodes$name, description = nodes$description,
significance = pval[rownames(nodes)], score = nodes$score,
type = nodes$type)
if (igraph::is.directed(subnet)) {
subg <- igraph::graph.data.frame(d = relations, directed = TRUE,
vertices = nodes)
}
else {
subg <- igraph::graph.data.frame(d = relations, directed = FALSE,
vertices = nodes)
}
if (!is.null(subnet$combinedP)) {
subg$combinedP <- subnet$combinedP
}
subg$threshold <- subnet$threshold
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
runTime <- as.numeric(difftime(strptime(endT, "\%Y-\%m-\%d
\%H:\%M:\%S"),
strptime(startT, "\%Y-\%m-\%d \%H:\%M:\%S"), units = "secs"))
if (!silent) {
message(paste(c("\nEnd at ", as.character(endT)), collapse = ""),
appendLF = TRUE)
message(paste(c("Runtime in total (xSubneterGenes): ",
runTime, " secs\n"), collapse = ""), appendLF = TRUE)
}
return(subg)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
