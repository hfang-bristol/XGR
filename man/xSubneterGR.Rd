\name{xSubneterGR}
\alias{xSubneterGR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSubneterGR(data, significance.threshold = 5e-05, score.cap = 10,
build.conversion = c(NA, "hg38.to.hg19", "hg18.to.hg19"), distance.max
= 50000, decay.kernel = c("slow", "linear", "rapid", "constant"),
decay.exponent = 2, GR.Gene = c("UCSC_knownGene",
"UCSC_knownCanonical"), scoring.scheme = c("max", "sum", "sequential"),
network = c("STRING_highest", "STRING_high", "STRING_medium",
"STRING_low", "PCommonsUN_high", "PCommonsUN_medium",
"PCommonsDN_high", "PCommonsDN_medium", "PCommonsDN_Reactome",
"PCommonsDN_KEGG", "PCommonsDN_HumanCyc", "PCommonsDN_PID",
"PCommonsDN_PANTHER", "PCommonsDN_ReconX", "PCommonsDN_TRANSFAC",
"PCommonsDN_PhosphoSite", "PCommonsDN_CTD"), network.customised = NULL,
seed.genes = T, subnet.significance = 5e-05, subnet.size = NULL,
verbose = T, RData.location = "http://galahad.well.ox.ac.uk/bigdata")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{significance.threshold}{
%% ~~Describe \code{significance.threshold} here~~
}
\item{score.cap}{
%% ~~Describe \code{score.cap} here~~
}
\item{build.conversion}{
%% ~~Describe \code{build.conversion} here~~
}
\item{distance.max}{
%% ~~Describe \code{distance.max} here~~
}
\item{decay.kernel}{
%% ~~Describe \code{decay.kernel} here~~
}
\item{decay.exponent}{
%% ~~Describe \code{decay.exponent} here~~
}
\item{GR.Gene}{
%% ~~Describe \code{GR.Gene} here~~
}
\item{scoring.scheme}{
%% ~~Describe \code{scoring.scheme} here~~
}
\item{network}{
%% ~~Describe \code{network} here~~
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
function (data, significance.threshold = 5e-05, score.cap = 10,
build.conversion = c(NA, "hg38.to.hg19", "hg18.to.hg19"),
distance.max = 50000, decay.kernel = c("slow", "linear",
"rapid", "constant"), decay.exponent = 2, GR.Gene = c("UCSC_knownGene",
"UCSC_knownCanonical"), scoring.scheme = c("max", "sum",
"sequential"), network = c("STRING_highest", "STRING_high",
"STRING_medium", "STRING_low", "PCommonsUN_high", "PCommonsUN_medium",
"PCommonsDN_high", "PCommonsDN_medium", "PCommonsDN_Reactome",
"PCommonsDN_KEGG", "PCommonsDN_HumanCyc", "PCommonsDN_PID",
"PCommonsDN_PANTHER", "PCommonsDN_ReconX", "PCommonsDN_TRANSFAC",
"PCommonsDN_PhosphoSite", "PCommonsDN_CTD"), network.customised = NULL,
seed.genes = T, subnet.significance = 5e-05, subnet.size = NULL,
verbose = T, RData.location = "http://galahad.well.ox.ac.uk/bigdata")
{
startT <- Sys.time()
if (verbose) {
message(paste(c("Start at ", as.character(startT)), collapse = ""),
appendLF = T)
message("", appendLF = T)
}
build.conversion <- match.arg(build.conversion)
decay.kernel <- match.arg(decay.kernel)
scoring.scheme <- match.arg(scoring.scheme)
network <- match.arg(network)
if (verbose) {
now <- Sys.time()
message(sprintf("\n#######################################################",
appendLF = T))
message(sprintf("'xGR2GeneScores' is being called to score seed genes
(\%s):",
as.character(now)), appendLF = T)
message(sprintf("#######################################################",
appendLF = T))
}
mSeed <- xGR2GeneScores(data = data, significance.threshold =
significance.threshold,
build.conversion = build.conversion, distance.max = distance.max,
decay.kernel = decay.kernel, decay.exponent = decay.exponent,
GR.Gene = GR.Gene, scoring.scheme = scoring.scheme, verbose = verbose,
RData.location = RData.location)
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("'xGR2GeneScores' has been finished (\%s)!",
as.character(now)), appendLF = T)
message(sprintf("#######################################################\n",
appendLF = T))
}
df_Gene <- mSeed$Gene
pval <- df_Gene$Pval
names(pval) <- df_Gene$Gene
if (verbose) {
now <- Sys.time()
message(sprintf("\t\t minimum p-value: \%1.2e; maximum p-value:
\%1.2e",
min(pval), max(pval)), appendLF = T)
}
if (verbose) {
now <- Sys.time()
message(sprintf("\n#######################################################",
appendLF = T))
message(sprintf("xSubneterGenes is being called (\%s):",
as.character(now)), appendLF = T)
message(sprintf("#######################################################",
appendLF = T))
}
subg <- xSubneterGenes(data = pval, network = network,
network.customised = network.customised,
seed.genes = seed.genes, subnet.significance = subnet.significance,
subnet.size = subnet.size, verbose = verbose, RData.location =
RData.location)
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("xSubneterGenes has finished (\%s)!",
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
