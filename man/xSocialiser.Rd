\name{xSocialiser}
\alias{xSocialiser}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSocialiser(data, annotation, g, measure = c("BM.average", "BM.max",
"BM.complete", "average", "max"), method.term = c("Resnik", "Lin",
"Schlicker", "Jiang", "Pesquita"), rescale = TRUE, force = TRUE, fast =
TRUE, parallel = TRUE, multicores = NULL, path.mode = c("all_paths",
"shortest_paths", "all_shortest_paths"), true.path.rule = TRUE, verbose
= T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{annotation}{
%% ~~Describe \code{annotation} here~~
}
\item{g}{
%% ~~Describe \code{g} here~~
}
\item{measure}{
%% ~~Describe \code{measure} here~~
}
\item{method.term}{
%% ~~Describe \code{method.term} here~~
}
\item{rescale}{
%% ~~Describe \code{rescale} here~~
}
\item{force}{
%% ~~Describe \code{force} here~~
}
\item{fast}{
%% ~~Describe \code{fast} here~~
}
\item{parallel}{
%% ~~Describe \code{parallel} here~~
}
\item{multicores}{
%% ~~Describe \code{multicores} here~~
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
function (data, annotation, g, measure = c("BM.average", "BM.max",
"BM.complete", "average", "max"), method.term = c("Resnik",
"Lin", "Schlicker", "Jiang", "Pesquita"), rescale = TRUE,
force = TRUE, fast = TRUE, parallel = TRUE, multicores = NULL,
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
true.path.rule = TRUE, verbose = T)
{
measure <- match.arg(measure)
method.term <- match.arg(method.term)
path.mode <- match.arg(path.mode)
if (class(annotation) == "GS") {
originAnnos <- annotation$gs
}
else if (class(annotation) == "list") {
originAnnos <- annotation
}
else if (class(annotation) == "dgCMatrix") {
D <- annotation
originAnnos <- sapply(1:ncol(D), function(j) {
names(which(D[, j] != 0))
})
names(originAnnos) <- colnames(annotation)
}
else {
stop("The input annotation must be either 'GS' or 'list' or 'dgCMatrix'
object.\n")
}
annotation <- originAnnos
ig <- g
if (class(ig) != "igraph") {
stop("The function must apply to the 'igraph' object.\n")
}
else {
if (verbose) {
now <- Sys.time()
message(sprintf("First, generate a subgraph induced (via '\%s' mode) by
the annotation data (\%s) ...",
path.mode, as.character(now)), appendLF = T)
}
ig <- xDAGanno(g = ig, annotation = annotation, path.mode = path.mode,
true.path.rule = true.path.rule, verbose = verbose)
}
anno <- V(ig)$anno
allSNPs <- sort(unique(unlist(anno)))
SNPs <- data[!is.na(data)]
if (is.null(SNPs) || is.na(SNPs)) {
}
else {
flag <- SNPs \%in\% allSNPs
if (sum(flag) != 0) {
SNPs <- SNPs[flag]
}
else {
SNPs <- SNPs[flag]
}
}
if (length(SNPs) <= 1 || is.na(SNPs)) {
res <- NULL
warning("The function returns NULL as no similarity is found at
all.\n")
return(res)
}
allterms <- 1:length(anno)
sGT <- Matrix::Matrix(0, nrow = length(SNPs), ncol = length(allterms),
sparse = T)
for (j in 1:length(allterms)) {
ind <- match(anno[[j]], SNPs)
flag <- ind[!is.na(ind)]
if (length(flag) != 0) {
sGT[flag, j] <- 1
}
}
colnames(sGT) <- V(ig)$name
rownames(sGT) <- SNPs
if (verbose) {
message(sprintf("\tthere are \%d inputs amongst \%d annotatable",
length(SNPs), length(allSNPs)), appendLF = T)
}
SNPs2terms <- lapply(1:length(SNPs), function(x) {
res <- names(which(sGT[x, ] == 1))
if (force) {
subg <- dnet::dDAGinduce(ig, nodes_query = res, path.mode =
"all_paths")
res <- dnet::dDAGtip(subg)
}
return(res)
})
names(SNPs2terms) <- SNPs
terms <- unique(unlist(SNPs2terms))
SNPs2terms_index <- lapply(SNPs2terms, function(x) {
match(x, terms)
})
progress_indicate <- function(i, B, step, flag = F) {
if (i\%\%ceiling(B/step) == 0 | i == B | i == 1) {
if (flag & verbose) {
message(sprintf("\t\%d out of \%d (\%s)", i, B,
as.character(Sys.time())), appendLF = T)
}
}
}
if (verbose) {
now <- Sys.time()
if (force) {
message(sprintf("Next, pre-compute semantic similarity between \%d
terms (forced to be the most specific for each gene/SNP) using '\%s'
method (\%s)...",
length(terms), method.term, as.character(now)),
appendLF = T)
}
else {
message(sprintf("Next, pre-compute semantic similarity between \%d
terms using '\%s' method (\%s)...",
length(terms), method.term, as.character(now)),
appendLF = T)
}
}
sim.term <- suppressMessages(xDAGsim(ig, terms = terms, method.term =
method.term,
fast = fast, parallel = parallel, multicores = multicores,
verbose = T))
if (class(sim.term) == "igraph") {
sim.term <- xConverter(sim.term, from = "igraph", to = "dgCMatrix",
verbose = F)
}
if (verbose) {
now <- Sys.time()
message(sprintf("Last, calculate pair-wise semantic similarity between
\%d genes/SNPs using '\%s' measure (\%s)...",
length(SNPs), measure, as.character(now)), appendLF = T)
}
num_SNPs <- length(SNPs2terms)
flag_parallel <- F
if (parallel == TRUE) {
flag_parallel <- xCheckParallel(multicores = multicores,
verbose = verbose)
if (flag_parallel) {
if (measure == "average") {
i <- 1
sim <- foreach::`\%dopar\%`(foreach::foreach(i = 1:(num_SNPs -
1), .inorder = T, .combine = rbind), {
ind1 <- SNPs2terms_index[[i]]
progress_indicate(i, num_SNPs, 10, flag = T)
fast <- T
if (fast) {
js <- (i + 1):num_SNPs
ind_js <- SNPs2terms_index[js]
sim12 <- matrix(sim.term[ind1, unlist(ind_js)],
nrow = length(ind1))
new_ind_js <- rep(1:length(ind_js), sapply(ind_js,
length))
res <- sapply(1:length(ind_js), function(k) {
mean(sim12[, which(new_ind_js == k)])
})
x <- rep(0, num_SNPs)
x[js] <- res
x
}
})
}
else if (measure == "max") {
i <- 1
sim <- foreach::`\%dopar\%`(foreach::foreach(i = 1:(num_SNPs -
1), .inorder = T, .combine = rbind), {
ind1 <- SNPs2terms_index[[i]]
progress_indicate(i, num_SNPs, 10, flag = T)
fast <- T
if (fast) {
js <- (i + 1):num_SNPs
ind_js <- SNPs2terms_index[js]
sim12 <- matrix(sim.term[ind1, unlist(ind_js)],
nrow = length(ind1))
new_ind_js <- rep(1:length(ind_js), sapply(ind_js,
length))
res <- sapply(1:length(ind_js), function(k) {
max(sim12[, which(new_ind_js == k)])
})
x <- rep(0, num_SNPs)
x[js] <- res
x
}
})
}
else if (measure == "BM.average") {
i <- 1
sim <- foreach::`\%dopar\%`(foreach::foreach(i = 1:(num_SNPs -
1), .inorder = T, .combine = rbind), {
ind1 <- SNPs2terms_index[[i]]
progress_indicate(i, num_SNPs, 10, flag = T)
fast <- T
if (fast) {
js <- (i + 1):num_SNPs
ind_js <- SNPs2terms_index[js]
sim12 <- matrix(sim.term[ind1, unlist(ind_js)],
nrow = length(ind1))
new_ind_js <- rep(1:length(ind_js), sapply(ind_js,
length))
res <- sapply(1:length(ind_js), function(k) {
x <- as.matrix(sim12[, which(new_ind_js ==
k)])
0.5 * (mean(apply(x, 1, max)) + mean(apply(x,
2, max)))
})
x <- rep(0, num_SNPs)
x[js] <- res
x
}
})
}
else if (measure == "BM.max") {
i <- 1
sim <- foreach::`\%dopar\%`(foreach::foreach(i = 1:(num_SNPs -
1), .inorder = T, .combine = rbind), {
ind1 <- SNPs2terms_index[[i]]
progress_indicate(i, num_SNPs, 10, flag = T)
fast <- T
if (fast) {
js <- (i + 1):num_SNPs
ind_js <- SNPs2terms_index[js]
sim12 <- matrix(sim.term[ind1, unlist(ind_js)],
nrow = length(ind1))
new_ind_js <- rep(1:length(ind_js), sapply(ind_js,
length))
res <- sapply(1:length(ind_js), function(k) {
x <- as.matrix(sim12[, which(new_ind_js ==
k)])
max(mean(apply(x, 1, max)), mean(apply(x,
2, max)))
})
x <- rep(0, num_SNPs)
x[js] <- res
x
}
})
}
else if (measure == "BM.complete") {
i <- 1
sim <- foreach::`\%dopar\%`(foreach::foreach(i = 1:(num_SNPs -
1), .inorder = T, .combine = rbind), {
ind1 <- SNPs2terms_index[[i]]
progress_indicate(i, num_SNPs, 10, flag = T)
fast <- T
if (fast) {
js <- (i + 1):num_SNPs
ind_js <- SNPs2terms_index[js]
sim12 <- matrix(sim.term[ind1, unlist(ind_js)],
nrow = length(ind1))
new_ind_js <- rep(1:length(ind_js), sapply(ind_js,
length))
res <- sapply(1:length(ind_js), function(k) {
x <- as.matrix(sim12[, which(new_ind_js ==
k)])
min(c(apply(x, 1, max), apply(x, 2, max)))
})
x <- rep(0, num_SNPs)
x[js] <- res
x
}
})
}
sim <- rbind(sim, rep(0, num_SNPs))
sim <- sim + Matrix::t(sim)
sim <- Matrix::Matrix(sim, sparse = T)
}
}
if (flag_parallel == F) {
sim <- Matrix::Matrix(0, nrow = length(SNPs), ncol = length(SNPs),
sparse = T)
if (measure == "average") {
for (i in 1:(num_SNPs - 1)) {
ind1 <- SNPs2terms_index[[i]]
progress_indicate(i, num_SNPs, 10, flag = T)
if (fast) {
js <- (i + 1):num_SNPs
ind_js <- SNPs2terms_index[js]
sim12 <- matrix(sim.term[ind1, unlist(ind_js)],
nrow = length(ind1))
new_ind_js <- rep(1:length(ind_js), sapply(ind_js,
length))
res <- sapply(1:length(ind_js), function(k) {
mean(sim12[, which(new_ind_js == k)])
})
sim[i, js] <- res
}
else {
for (j in (i + 1):num_SNPs) {
ind2 <- SNPs2terms_index[[j]]
sim12 <- as.matrix(sim.term[ind1, ind2])
sim[i, j] <- mean(sim12)
}
}
}
}
else if (measure == "max") {
for (i in 1:(num_SNPs - 1)) {
ind1 <- SNPs2terms_index[[i]]
progress_indicate(i, num_SNPs, 10, flag = T)
if (fast) {
js <- (i + 1):num_SNPs
ind_js <- SNPs2terms_index[js]
sim12 <- matrix(sim.term[ind1, unlist(ind_js)],
nrow = length(ind1))
new_ind_js <- rep(1:length(ind_js), sapply(ind_js,
length))
res <- sapply(1:length(ind_js), function(k) {
max(sim12[, which(new_ind_js == k)])
})
sim[i, js] <- res
}
else {
for (j in (i + 1):num_SNPs) {
ind2 <- SNPs2terms_index[[j]]
sim12 <- as.matrix(sim.term[ind1, ind2])
sim[i, j] <- max(sim12)
}
}
}
}
else if (measure == "BM.average") {
for (i in 1:(num_SNPs - 1)) {
ind1 <- SNPs2terms_index[[i]]
progress_indicate(i, num_SNPs, 10, flag = T)
if (fast) {
js <- (i + 1):num_SNPs
ind_js <- SNPs2terms_index[js]
sim12 <- matrix(sim.term[ind1, unlist(ind_js)],
nrow = length(ind1))
new_ind_js <- rep(1:length(ind_js), sapply(ind_js,
length))
res <- sapply(1:length(ind_js), function(k) {
x <- as.matrix(sim12[, which(new_ind_js ==
k)])
0.5 * (mean(apply(x, 1, max)) + mean(apply(x,
2, max)))
})
sim[i, js] <- res
}
else {
for (j in (i + 1):num_SNPs) {
ind2 <- SNPs2terms_index[[j]]
sim12 <- as.matrix(sim.term[ind1, ind2])
sim[i, j] <- 0.5 * (mean(apply(sim12, 1,
max)) + mean(apply(sim12, 2, max)))
}
}
}
}
else if (measure == "BM.max") {
for (i in 1:(num_SNPs - 1)) {
ind1 <- SNPs2terms_index[[i]]
progress_indicate(i, num_SNPs, 10, flag = T)
if (fast) {
js <- (i + 1):num_SNPs
ind_js <- SNPs2terms_index[js]
sim12 <- matrix(sim.term[ind1, unlist(ind_js)],
nrow = length(ind1))
new_ind_js <- rep(1:length(ind_js), sapply(ind_js,
length))
res <- sapply(1:length(ind_js), function(k) {
x <- as.matrix(sim12[, which(new_ind_js ==
k)])
max(mean(apply(x, 1, max)), mean(apply(x,
2, max)))
})
sim[i, js] <- res
}
else {
for (j in (i + 1):num_SNPs) {
ind2 <- SNPs2terms_index[[j]]
sim12 <- as.matrix(sim.term[ind1, ind2])
sim[i, j] <- max(mean(apply(sim12, 1, max)),
mean(apply(sim12, 2, max)))
}
}
}
}
else if (measure == "BM.complete") {
for (i in 1:(num_SNPs - 1)) {
ind1 <- SNPs2terms_index[[i]]
progress_indicate(i, num_SNPs, 10, flag = T)
if (fast) {
js <- (i + 1):num_SNPs
ind_js <- SNPs2terms_index[js]
sim12 <- matrix(sim.term[ind1, unlist(ind_js)],
nrow = length(ind1))
new_ind_js <- rep(1:length(ind_js), sapply(ind_js,
length))
res <- sapply(1:length(ind_js), function(k) {
x <- as.matrix(sim12[, which(new_ind_js ==
k)])
min(c(apply(x, 1, max), apply(x, 2, max)))
})
sim[i, js] <- res
}
else {
for (j in (i + 1):num_SNPs) {
ind2 <- SNPs2terms_index[[j]]
sim12 <- as.matrix(sim.term[ind1, ind2])
sim[i, j] <- min(c(apply(sim12, 1, max),
apply(sim12, 2, max)))
}
}
}
}
sim <- sim + Matrix::t(sim)
}
rownames(sim) <- colnames(sim) <- SNPs
sim[as.matrix(is.na(sim))] <- 0
if (rescale) {
if (verbose) {
now <- Sys.time()
message(sprintf("Also rescale similarity into the [0,1] range (\%s)",
as.character(now)), appendLF = T)
}
d_min <- min(sim)
d_max <- max(sim)
sim <- apply(sim, 2, function(x) {
(x - d_min)/(d_max - d_min)
})
sim <- Matrix::Matrix(sim, sparse = T)
}
if (class(sim) == "dgCMatrix" | class(sim) == "dsCMatrix") {
res <- xConverter(sim, from = "dgCMatrix", to = "igraph",
verbose = F)
if (class(res) == "igraph") {
if (!is.null(E(res)$weight)) {
E(res)$weight <- as.numeric(E(res)$weight)
}
}
res$dag <- ig
}
if (ecount(res) == 0) {
res <- NULL
warning("The function returns NULL as no similarity is found at
all.\n")
}
invisible(res)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
