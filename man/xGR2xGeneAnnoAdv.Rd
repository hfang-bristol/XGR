\name{xGR2xGeneAnnoAdv}
\alias{xGR2xGeneAnnoAdv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGR2xGeneAnnoAdv(list_vec, background = NULL, build.conversion = c(NA,
"hg38.to.hg19", "hg18.to.hg19"), crosslink = c("genehancer",
"PCHiC_combined", "GTEx_V6p_combined", "nearby"), crosslink.customised
= NULL, crosslink.top = NULL, nearby.distance.max = 50000,
nearby.decay.kernel = c("rapid", "slow", "linear", "constant"),
nearby.decay.exponent = 2, ontologies = NA, size.range = c(10, 2000),
min.overlap = 5, which.distance = NULL, test = c("hypergeo", "fisher",
"binomial"), background.annotatable.only = NULL, p.tail = c("one-tail",
"two-tails"), p.adjust.method = c("BH", "BY", "bonferroni", "holm",
"hochberg", "hommel"), ontology.algorithm = c("none", "pc", "elim",
"lea"), elim.pvalue = 0.01, lea.depth = 2, path.mode = c("all_paths",
"shortest_paths", "all_shortest_paths"), true.path.rule = F, verbose =
F, silent = F, plot = T, fdr.cutoff = 0.05, displayBy = c("zscore",
"fdr", "pvalue", "fc", "or"), RData.location =
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
\item{build.conversion}{
%% ~~Describe \code{build.conversion} here~~
}
\item{crosslink}{
%% ~~Describe \code{crosslink} here~~
}
\item{crosslink.customised}{
%% ~~Describe \code{crosslink.customised} here~~
}
\item{crosslink.top}{
%% ~~Describe \code{crosslink.top} here~~
}
\item{nearby.distance.max}{
%% ~~Describe \code{nearby.distance.max} here~~
}
\item{nearby.decay.kernel}{
%% ~~Describe \code{nearby.decay.kernel} here~~
}
\item{nearby.decay.exponent}{
%% ~~Describe \code{nearby.decay.exponent} here~~
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
function (list_vec, background = NULL, build.conversion = c(NA,
"hg38.to.hg19", "hg18.to.hg19"), crosslink = c("genehancer",
"PCHiC_combined", "GTEx_V6p_combined", "nearby"), crosslink.customised
= NULL,
crosslink.top = NULL, nearby.distance.max = 50000, nearby.decay.kernel
= c("rapid",
"slow", "linear", "constant"), nearby.decay.exponent = 2,
ontologies = NA, size.range = c(10, 2000), min.overlap = 5,
which.distance = NULL, test = c("hypergeo", "fisher", "binomial"),
background.annotatable.only = NULL, p.tail = c("one-tail",
"two-tails"), p.adjust.method = c("BH", "BY", "bonferroni",
"holm", "hochberg", "hommel"), ontology.algorithm = c("none",
"pc", "elim", "lea"), elim.pvalue = 0.01, lea.depth = 2,
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
true.path.rule = F, verbose = F, silent = F, plot = T, fdr.cutoff =
0.05,
displayBy = c("zscore", "fdr", "pvalue", "fc", "or"), RData.location =
"http://galahad.well.ox.ac.uk/bigdata",
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
build.conversion <- match.arg(build.conversion)
nearby.decay.kernel <- match.arg(nearby.decay.kernel)
test <- match.arg(test)
p.tail <- match.arg(p.tail)
p.adjust.method <- match.arg(p.adjust.method)
ontology.algorithm <- match.arg(ontology.algorithm)
path.mode <- match.arg(path.mode)
displayBy <- match.arg(displayBy)
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
fast <- T
if (!fast) {
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
eTerm <- xGR2xGeneAnno(data = data, background = background,
format = "chr:start-end", build.conversion = build.conversion,
crosslink = crosslink, crosslink.customised = crosslink.customised,
crosslink.top = crosslink.top, nearby.distance.max =
nearby.distance.max,
nearby.decay.kernel = nearby.decay.kernel,
nearby.decay.exponent = nearby.decay.exponent,
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
}
else {
Score <- Gene <- NULL
if (verbose) {
now <- Sys.time()
message(sprintf("First, import the background (\%s) ...",
as.character(now)), appendLF = T)
}
bGR <- xGR(data = background, format = "chr:start-end",
build.conversion = build.conversion, verbose = verbose,
RData.location = RData.location, guid = guid)
if (verbose) {
now <- Sys.time()
message(sprintf("Second, define crosslinked genes based on '\%s' (\%s)
...",
crosslink, as.character(now)), appendLF = T)
}
df_xGenes_background <- xGR2xGenes(data = bGR, format = "GRanges",
crosslink = crosslink, crosslink.customised = crosslink.customised,
cdf.function = "original", scoring = TRUE, scoring.scheme = "max",
scoring.rescale = F, nearby.distance.max = nearby.distance.max,
nearby.decay.kernel = nearby.decay.kernel, nearby.decay.exponent =
nearby.decay.exponent,
verbose = verbose, silent = !verbose, RData.location = RData.location,
guid = guid)
if (!is.null(df_xGenes_background)) {
bGR_genes <- (df_xGenes_background \%>\% dplyr::arrange(-Score))$Gene
}
else {
bGR_genes <- NULL
}
ls_df <- lapply(1:length(list_vec), function(i) {
if (verbose) {
message(sprintf("Analysing group \%d ('\%s') (\%s) ...",
i, names(list_vec)[i], as.character(Sys.time())),
appendLF = T)
}
data <- list_vec[[i]]
dGR <- xGR(data = data, format = "chr:start-end",
build.conversion = build.conversion, verbose = verbose,
RData.location = RData.location, guid = guid)
df_xGenes_data <- xGR2xGenes(data = dGR, format = "GRanges",
crosslink = crosslink, crosslink.customised = crosslink.customised,
cdf.function = "original", scoring = TRUE, scoring.scheme = "max",
scoring.rescale = F, nearby.distance.max = nearby.distance.max,
nearby.decay.kernel = nearby.decay.kernel, nearby.decay.exponent =
nearby.decay.exponent,
verbose = verbose, silent = !verbose, RData.location = RData.location,
guid = guid)
df_xGenes_data <- df_xGenes_data \%>\% dplyr::arrange(-Score)
if (is.null(crosslink.top)) {
crosslink.top <- nrow(df_xGenes_data)
}
if (crosslink.top > nrow(df_xGenes_data)) {
crosslink.top <- nrow(df_xGenes_data)
}
crosslink.top <- as.integer(crosslink.top)
crosslink.cutoff <- df_xGenes_data[crosslink.top,
"Score"]
dGR_genes <- df_xGenes_data$Gene[df_xGenes_data$Score >=
crosslink.cutoff]
if (verbose) {
if (is.null(bGR_genes)) {
message(sprintf("\t\%d (out of \%d crosslinked genes) are used.",
length(dGR_genes), nrow(df_xGenes_data),
as.character(Sys.time())), appendLF = T)
}
else {
message(sprintf("\t\%d (out of \%d crosslinked genes) and \%d
background genes are used.",
length(dGR_genes), nrow(df_xGenes_data),
length(bGR_genes), as.character(Sys.time())),
appendLF = T)
}
}
ls_df <- lapply(1:length(ontologies), function(j) {
if (verbose) {
message(sprintf("\tontology \%d ('\%s') (\%s) ...",
j, ontologies[j], as.character(Sys.time())),
appendLF = T)
}
ontology <- ontologies[j]
if (verbose) {
now <- Sys.time()
message(sprintf("\n#######################################################",
appendLF = T))
message(sprintf("'xEnricherGenes' is being called (\%s):",
as.character(now)), appendLF = T)
message(sprintf("#######################################################",
appendLF = T))
}
eTerm <- xEnricherGenes(data = dGR_genes, background = bGR_genes,
ontology = ontology, size.range = size.range,
min.overlap = min.overlap, which.distance = which.distance,
test = test, background.annotatable.only = background.annotatable.only,
p.tail = p.tail, p.adjust.method = p.adjust.method,
ontology.algorithm = ontology.algorithm, elim.pvalue = elim.pvalue,
lea.depth = lea.depth, path.mode = path.mode,
true.path.rule = true.path.rule, verbose = verbose,
silent = !verbose, RData.location = RData.location,
guid = guid)
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("'xEnricherGenes' has been finished (\%s)!",
as.character(now)), appendLF = T)
message(sprintf("#######################################################\n",
appendLF = T))
}
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
}
if (plot & !is.null(df_all)) {
adjp <- NULL
gp <- NULL
mat <- NULL
ls_df <- split(x = df_all[, -12], f = df_all$ontology)
ls_df <- ls_df[unique(df_all$ontology)]
ls_mat <- lapply(1:length(ls_df), function(i) {
df <- ls_df[[i]]
ind <- which(df$adjp < fdr.cutoff)
if (length(ind) >= 1) {
df <- as.data.frame(df \%>\% dplyr::filter(adjp <
fdr.cutoff))
if (displayBy == "fdr") {
mat <- as.matrix(xSparseMatrix(df[, c("name",
"group", "adjp")], rows = unique(df$name),
columns = names(list_vec)))
mat[mat == 0] <- NA
mat <- -log10(mat)
}
else if (displayBy == "pvalue") {
mat <- as.matrix(xSparseMatrix(df[, c("name",
"group", "pvalue")], rows = unique(df$name),
columns = names(list_vec)))
mat[mat == 0] <- NA
mat <- -log10(mat)
}
else if (displayBy == "zscore") {
mat <- as.matrix(xSparseMatrix(df[, c("name",
"group", "zscore")], rows = unique(df$name),
columns = names(list_vec)))
mat[mat == 0] <- NA
}
else if (displayBy == "fc") {
mat <- as.matrix(xSparseMatrix(df[, c("name",
"group", "fc")], rows = unique(df$name),
columns = names(list_vec)))
mat[mat == 0] <- NA
mat <- log2(mat)
}
else if (displayBy == "or") {
mat <- as.matrix(xSparseMatrix(df[, c("name",
"group", "or")], rows = unique(df$name),
columns = names(list_vec)))
mat[mat == 0] <- NA
mat <- log2(mat)
}
if (nrow(mat) == 1) {
df_mat <- mat
}
else {
rname_ordered <- rownames(mat)[order(-nchar(rownames(mat)))]
if (names(ls_df)[i] == "PS2") {
df_tmp <- unique(df[, c("id", "name")])
df_tmp <- df_tmp[with(df_tmp, order(as.numeric(df_tmp$id))),
]
rname_ordered <- df_tmp$name
}
ind <- match(rname_ordered, rownames(mat))
df_mat <- as.matrix(mat[ind, ], ncol = ncol(mat))
colnames(df_mat) <- colnames(mat)
colnames(df_mat) <- colnames(mat)
}
return(df_mat)
}
else {
return(NULL)
}
})
mat <- do.call(rbind, ls_mat)
if (!is.null(mat)) {
if (displayBy == "fdr" | displayBy == "pvalue") {
colormap <- "grey100-darkorange"
zlim <- c(0, ceiling(max(mat[!is.na(mat)])))
if (displayBy == "fdr") {
legend.title <- expression(-log[10]("FDR"))
}
else if (displayBy == "pvalue") {
legend.title <- expression(-log[10]("p-value"))
}
}
else if (displayBy == "fc" | displayBy == "zscore" |
displayBy == "or") {
tmp_max <- ceiling(max(mat[!is.na(mat)]))
tmp_min <- floor(min(mat[!is.na(mat)]))
if (tmp_max > 0 & tmp_min < 0) {
colormap <- "deepskyblue-grey100-darkorange"
tmp <- max(tmp_max, abs(tmp_min))
zlim <- c(-tmp, tmp)
}
else if (tmp_max <= 0) {
colormap <- "deepskyblue-grey100"
zlim <- c(tmp_min, 0)
}
else if (tmp_min >= 0) {
colormap <- "grey100-darkorange"
zlim <- c(0, tmp_max)
}
if (displayBy == "fc") {
legend.title <- expression(log[2]("FC"))
}
else if (displayBy == "zscore") {
legend.title <- ("Z-score")
}
else if (displayBy == "or") {
legend.title <- expression(log[2]("OR"))
}
}
gp <- xHeatmap(mat, reorder = "none", colormap = colormap,
ncolors = 64, zlim = zlim, legend.title = legend.title,
barwidth = 0.4, x.rotate = 60, shape = 19, size = 2,
x.text.size = 6, y.text.size = 6, na.color = "transparent",
barheight = max(3, min(5, nrow(mat))))
gp <- gp + theme(legend.title = element_text(size = 8))
}
}
else {
mat <- NULL
gp <- NULL
}
ls_eTerm <- NULL
if (!is.null(df_all)) {
ls_eTerm <- list(df = df_all, mat = mat, gp = gp)
class(ls_eTerm) <- "ls_eTerm"
}
endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, "\%Y-\%m-\%d
\%H:\%M:\%S"),
strptime(startT, "\%Y-\%m-\%d \%H:\%M:\%S"), units = "secs"))
if (!silent) {
message(paste(c("\nEnd at ", as.character(endT)), collapse = ""),
appendLF = TRUE)
message(paste(c("Runtime in total (xGR2xGeneAnnoAdv): ",
runTime, " secs\n"), collapse = ""), appendLF = TRUE)
}
invisible(ls_eTerm)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
