\name{xEnricherGenesAdv}
\alias{xEnricherGenesAdv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnricherGenesAdv(list_vec, background = NULL, check.symbol.identity =
F, ontologies = NA, size.range = c(10, 2000), min.overlap = 5,
which.distance = NULL, test = c("fisher", "hypergeo", "binomial"),
background.annotatable.only = NULL, p.tail = c("one-tail",
"two-tails"), p.adjust.method = c("BH", "BY", "bonferroni", "holm",
"hochberg", "hommel"), ontology.algorithm = c("none", "pc", "elim",
"lea"), elim.pvalue = 0.01, lea.depth = 2, path.mode = c("all_paths",
"shortest_paths", "all_shortest_paths"), true.path.rule = F, verbose =
F, silent = FALSE, plot = TRUE, fdr.cutoff = 0.05, displayBy =
c("zscore", "fdr", "pvalue", "fc", "or"), RData.location =
"http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{list_vec}{
%% ~~Describe \code{list_vec} here~~
}
\item{background}{
%% ~~Describe \code{background} here~~
}
\item{check.symbol.identity}{
%% ~~Describe \code{check.symbol.identity} here~~
}
\item{ontologies}{
%% ~~Describe \code{ontologies} here~~
}
\item{size.range}{
%% ~~Describe \code{size.range} here~~
}
\item{min.overlap}{
%% ~~Describe \code{min.overlap} here~~
}
\item{which.distance}{
%% ~~Describe \code{which.distance} here~~
}
\item{test}{
%% ~~Describe \code{test} here~~
}
\item{background.annotatable.only}{
%% ~~Describe \code{background.annotatable.only} here~~
}
\item{p.tail}{
%% ~~Describe \code{p.tail} here~~
}
\item{p.adjust.method}{
%% ~~Describe \code{p.adjust.method} here~~
}
\item{ontology.algorithm}{
%% ~~Describe \code{ontology.algorithm} here~~
}
\item{elim.pvalue}{
%% ~~Describe \code{elim.pvalue} here~~
}
\item{lea.depth}{
%% ~~Describe \code{lea.depth} here~~
}
\item{path.mode}{
%% ~~Describe \code{path.mode} here~~
}
\item{true.path.rule}{
%% ~~Describe \code{true.path.rule} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
}
\item{silent}{
%% ~~Describe \code{silent} here~~
}
\item{plot}{
%% ~~Describe \code{plot} here~~
}
\item{fdr.cutoff}{
%% ~~Describe \code{fdr.cutoff} here~~
}
\item{displayBy}{
%% ~~Describe \code{displayBy} here~~
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
function (list_vec, background = NULL, check.symbol.identity = F,
ontologies = NA, size.range = c(10, 2000), min.overlap = 5,
which.distance = NULL, test = c("fisher", "hypergeo", "binomial"),
background.annotatable.only = NULL, p.tail = c("one-tail",
"two-tails"), p.adjust.method = c("BH", "BY", "bonferroni",
"holm", "hochberg", "hommel"), ontology.algorithm = c("none",
"pc", "elim", "lea"), elim.pvalue = 0.01, lea.depth = 2,
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
true.path.rule = F, verbose = F, silent = FALSE, plot = TRUE,
fdr.cutoff = 0.05, displayBy = c("zscore", "fdr", "pvalue",
"fc", "or"), RData.location = "http://galahad.well.ox.ac.uk/bigdata",
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
test <- match.arg(test)
p.tail <- match.arg(p.tail)
p.adjust.method <- match.arg(p.adjust.method)
ontology.algorithm <- match.arg(ontology.algorithm)
path.mode <- match.arg(path.mode)
p.tail <- match.arg(p.tail)
displayBy <- match.arg(displayBy)
if (class(list_vec) == "cModule") {
list_vec <- split(x = list_vec$mem$nodes, f = list_vec$mem$modules)
}
if (length(list_vec) == 0) {
return(NULL)
}
if (is.vector(list_vec) & class(list_vec) != "list") {
list_vec <- list(list_vec)
}
else if (class(list_vec) == "list") {
list_vec <- base::Filter(base::Negate(is.null), list_vec)
if (length(list_vec) == 0) {
return(NULL)
}
}
else {
stop("The input data must be a vector or a list of vectors.\n")
}
list_names <- names(list_vec)
if (is.null(list_names)) {
list_names <- paste0("G", 1:length(list_vec))
names(list_vec) <- list_names
}
ls_df <- lapply(1:length(list_vec), function(i) {
if (verbose) {
message(sprintf("Analysing group \%d ('\%s') (\%s) ...",
i, names(list_vec)[i], as.character(Sys.time())),
appendLF = T)
}
data <- list_vec[[i]]
ls_df <- lapply(1:length(ontologies), function(j) {
if (verbose) {
message(sprintf("\tontology \%d ('\%s') (\%s) ...",
j, ontologies[j], as.character(Sys.time())),
appendLF = T)
}
ontology <- ontologies[j]
eTerm <- xEnricherGenes(data = data, background = background,
check.symbol.identity = check.symbol.identity,
ontology = ontology, size.range = size.range,
min.overlap = min.overlap, which.distance = which.distance,
test = test, background.annotatable.only = background.annotatable.only,
p.tail = p.tail, p.adjust.method = p.adjust.method,
ontology.algorithm = ontology.algorithm, elim.pvalue = elim.pvalue,
lea.depth = lea.depth, path.mode = path.mode,
true.path.rule = true.path.rule, verbose = verbose,
silent = !verbose, RData.location = RData.location,
guid = guid)
df <- xEnrichViewer(eTerm, top_num = "all", sortBy = "or",
details = TRUE)
if (is.null(df)) {
return(NULL)
}
else {
cbind(group = rep(names(list_vec)[i], nrow(df)),
ontology = rep(ontology, nrow(df)), id = rownames(df),
df, stringsAsFactors = F)
}
})
df <- do.call(rbind, ls_df)
})
df_all <- do.call(rbind, ls_df)
df_all$group <- factor(df_all$group, levels = names(list_vec))
if (plot & !is.null(df_all)) {
gp <- xEnrichHeatmap(list_eTerm = df_all, fdr.cutoff = fdr.cutoff,
displayBy = displayBy, colormap = NULL, zlim = NULL,
reorder = "none")
mat <- gp$mat
gp$mat <- NULL
}
else {
mat <- NULL
gp <- NULL
}
ls_eTerm <- list(df = df_all, mat = mat, gp = gp)
class(ls_eTerm) <- "ls_eTerm"
endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, "\%Y-\%m-\%d
\%H:\%M:\%S"),
strptime(startT, "\%Y-\%m-\%d \%H:\%M:\%S"), units = "secs"))
if (!silent) {
message(paste(c("\nEnd at ", as.character(endT)), collapse = ""),
appendLF = TRUE)
message(paste(c("Runtime in total (xEnricherGenesAdv): ",
runTime, " secs\n"), collapse = ""), appendLF = TRUE)
}
invisible(ls_eTerm)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
