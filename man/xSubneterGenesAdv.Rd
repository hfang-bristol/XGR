\name{xSubneterGenesAdv}
\alias{xSubneterGenesAdv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSubneterGenesAdv(data, network = c("STRING_highest", "STRING_high",
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
subnet.size = 50, test.permutation = F, num.permutation = 100, respect
= c("none", "degree"), aggregateBy = c("Ztransform", "fishers",
"logistic", "orderStatistic"), num.subnets = NULL, verbose = T, silent
= F, RData.location = "http://galahad.well.ox.ac.uk/bigdata", guid =
NULL)
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
\item{num.subnets}{
%% ~~Describe \code{num.subnets} here~~
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
network.customised = NULL, seed.genes = T, subnet.size = 50,
test.permutation = F, num.permutation = 100, respect = c("none",
"degree"), aggregateBy = c("Ztransform", "fishers", "logistic",
"orderStatistic"), num.subnets = NULL, verbose = T, silent = F,
RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
network <- match.arg(network)
respect <- match.arg(respect)
aggregateBy <- match.arg(aggregateBy)
subnet.size <- as.integer(subnet.size)
if (subnet.size <= 0) {
return(NULL)
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
if (is.null(num.subnets)) {
num.subnets <- 100
}
num.subnets <- as.integer(num.subnets)
ig_tmp <- g
ls_subg <- list()
k <- 1
flag <- T
while (flag && k <= num.subnets) {
subg_tmp <- xSubneterGenes(data = data, network.customised = ig_tmp,
seed.genes = seed.genes, subnet.significance = NULL,
subnet.size = subnet.size, test.permutation = test.permutation,
num.permutation = num.permutation, respect = respect,
aggregateBy = aggregateBy, verbose = verbose, silent = T,
RData.location = RData.location, guid = guid)
if (verbose) {
message(sprintf("Iteration \%d: thresholded at \%1.2e (\%s) ...",
k, subg_tmp$threshold, as.character(Sys.time())),
appendLF = T)
}
if (subg_tmp$threshold < 0.05) {
ls_subg[[k]] <- subg_tmp
k <- k + 1
ig_tmp <- difference(ig_tmp, subg_tmp)
}
else {
flag <- F
}
}
isubg <- list(g = g, ls_subg = ls_subg, call = match.call())
class(isubg) <- "iSubg"
invisible(isubg)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
