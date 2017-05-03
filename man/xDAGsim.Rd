\name{xDAGsim}
\alias{xDAGsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xDAGsim(g, terms = NULL, method.term = c("Resnik", "Lin", "Schlicker",
"Jiang", "Pesquita"), fast = T, parallel = TRUE, multicores = NULL,
verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{g}{
%% ~~Describe \code{g} here~~
}
\item{terms}{
%% ~~Describe \code{terms} here~~
}
\item{method.term}{
%% ~~Describe \code{method.term} here~~
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
function (g, terms = NULL, method.term = c("Resnik", "Lin",
"Schlicker",
"Jiang", "Pesquita"), fast = T, parallel = TRUE, multicores = NULL,
verbose = T)
{
method.term <- match.arg(method.term)
ig <- g
if (class(ig) != "igraph") {
stop("The function must apply to the 'igraph' object.\n")
}
if (is.null(V(ig)$anno) | is.null(V(ig)$IC)) {
stop("The function requires that input graph has already contained
annotation data.  Please first run 'xDAGanno'.\n")
}
IC <- V(ig)$IC
names(IC) <- V(ig)$name
terms <- terms[!is.na(terms)]
if (is.null(terms) || is.na(terms)) {
terms <- V(ig)$name
}
else {
flag <- terms \%in\% V(ig)$name
if (sum(flag) != 0) {
terms <- terms[flag]
}
else {
terms <- V(ig)$name
}
}
progress_indicate <- function(i, B, step, flag = F) {
if (i\%\%ceiling(B/step) == 0 | i == B | i == 1) {
if (flag & verbose) {
message(sprintf("\t\%d out of \%d (\%s)", i, B,
as.character(Sys.time())), appendLF = T)
}
}
}
num_terms <- length(terms)
if (verbose) {
now <- Sys.time()
message(sprintf("Calculate semantic similarity between \%d terms using
'\%s' method (\%s)...",
num_terms, method.term, as.character(now)), appendLF = T)
}
sCP <- dnet::dDAGancestor(ig, term1 = terms, term2 = NULL,
verbose = T)
allterms <- V(ig)$name
ind <- match(terms, allterms)
flag_parallel <- F
if (parallel == TRUE) {
flag_parallel <- xCheckParallel(multicores = multicores,
verbose = verbose)
if (flag_parallel) {
if (method.term == "Resnik") {
i <- 1
sim <- foreach::`\%dopar\%`(foreach::foreach(i = 1:num_terms,
.inorder = T, .combine = rbind), {
ancestor_i <- which(sCP[i, ] == 1)
progress_indicate(i, num_terms, 10, flag = T)
fast <- T
if (fast) {
mat <- sCP[i:num_terms, ]
ancestor_js <- which(matrix(mat == 1, nrow = num_terms -
i + 1), arr.ind = T)
flag <- is.element(ancestor_js[, 2], ancestor_i)
ca_js <- ancestor_js[flag, ]
ca_js_list <- split(ca_js[, 2], ca_js[, 1])
mica_js <- sapply(ca_js_list, function(x) {
x[which.max(IC[x])]
})
js <- as.numeric(names(ca_js_list)) + i -
1
res <- IC[mica_js]
x <- rep(0, num_terms)
x[js] <- res
x
}
})
}
else if (method.term == "Lin") {
i <- 1
sim <- foreach::`\%dopar\%`(foreach::foreach(i = 1:num_terms,
.inorder = T, .combine = rbind), {
ancestor_i <- which(sCP[i, ] == 1)
progress_indicate(i, num_terms, 10, flag = T)
fast <- T
if (fast) {
mat <- sCP[i:num_terms, ]
ancestor_js <- which(matrix(mat == 1, nrow = num_terms -
i + 1), arr.ind = T)
flag <- is.element(ancestor_js[, 2], ancestor_i)
ca_js <- ancestor_js[flag, ]
ca_js_list <- split(ca_js[, 2], ca_js[, 1])
mica_js <- sapply(ca_js_list, function(x) {
x[which.max(IC[x])]
})
js <- as.numeric(names(ca_js_list)) + i -
1
res <- 2 * IC[mica_js]/(IC[ind[i]] + IC[ind[js]])
x <- rep(0, num_terms)
x[js] <- res
x
}
})
}
else if (method.term == "Schlicker") {
i <- 1
sim <- foreach::`\%dopar\%`(foreach::foreach(i = 1:num_terms,
.inorder = T, .combine = rbind), {
ancestor_i <- which(sCP[i, ] == 1)
progress_indicate(i, num_terms, 10, flag = T)
fast <- T
if (fast) {
mat <- sCP[i:num_terms, ]
ancestor_js <- which(matrix(mat == 1, nrow = num_terms -
i + 1), arr.ind = T)
flag <- is.element(ancestor_js[, 2], ancestor_i)
ca_js <- ancestor_js[flag, ]
ca_js_list <- split(ca_js[, 2], ca_js[, 1])
mica_js <- sapply(ca_js_list, function(x) {
x[which.max(IC[x])]
})
js <- as.numeric(names(ca_js_list)) + i -
1
res <- (2 * IC[mica_js]/(IC[ind[i]] + IC[ind[js]])) *
(1 - 10^(-IC[mica_js]))
x <- rep(0, num_terms)
x[js] <- res
x
}
})
}
else if (method.term == "Jiang") {
i <- 1
sim <- foreach::`\%dopar\%`(foreach::foreach(i = 1:num_terms,
.inorder = T, .combine = rbind), {
ancestor_i <- which(sCP[i, ] == 1)
progress_indicate(i, num_terms, 10, flag = T)
fast <- T
if (fast) {
mat <- sCP[i:num_terms, ]
ancestor_js <- which(matrix(mat == 1, nrow = num_terms -
i + 1), arr.ind = T)
flag <- is.element(ancestor_js[, 2], ancestor_i)
ca_js <- ancestor_js[flag, ]
ca_js_list <- split(ca_js[, 2], ca_js[, 1])
mica_js <- sapply(ca_js_list, function(x) {
x[which.max(IC[x])]
})
js <- as.numeric(names(ca_js_list)) + i -
1
tmp <- IC[ind[i]] + IC[ind[js]] - 2 * IC[mica_js]
tmp[tmp > 1] <- 1
res <- 1 - tmp
x <- rep(0, num_terms)
x[js] <- res
x
}
})
}
else if (method.term == "Pesquita") {
i <- 1
sim <- foreach::`\%dopar\%`(foreach::foreach(i = 1:num_terms,
.inorder = T, .combine = rbind), {
ancestor_i <- which(sCP[i, ] == 1)
progress_indicate(i, num_terms, 10, flag = T)
fast <- T
if (fast) {
mat <- sCP[i:num_terms, ]
ancestor_js <- which(matrix(mat == 1, nrow = num_terms -
i + 1), arr.ind = T)
flag <- is.element(ancestor_js[, 2], ancestor_i)
ca_js <- ancestor_js[flag, ]
ca_js_list <- split(ca_js[, 2], ca_js[, 1])
mica_js <- sapply(ca_js_list, function(x) {
x[which.max(IC[x])]
})
js <- as.numeric(names(ca_js_list)) + i -
1
allan_js_list <- split(ancestor_js[, 2],
ancestor_js[, 1])
allan_union <- sapply(allan_js_list, function(x) {
ux <- base::union(x, ancestor_i)
sum(IC[ux])
})
allca_union <- sapply(ca_js_list, function(x) {
sum(IC[x])
})
res <- allca_union/allan_union
x <- rep(0, num_terms)
x[js] <- res
x
}
})
}
sim <- sim + Matrix::t(sim)
sim <- Matrix::Matrix(sim, sparse = T)
}
}
if (flag_parallel == F) {
sim <- Matrix::Matrix(0, nrow = num_terms, ncol = num_terms,
sparse = T)
if (method.term == "Resnik") {
for (i in 1:num_terms) {
ancestor_i <- which(sCP[i, ] == 1)
progress_indicate(i, num_terms, 10, flag = T)
if (fast) {
mat <- sCP[i:num_terms, ]
ancestor_js <- which(matrix(mat == 1, nrow = num_terms -
i + 1), arr.ind = T)
flag <- is.element(ancestor_js[, 2], ancestor_i)
ca_js <- ancestor_js[flag, ]
ca_js_list <- split(ca_js[, 2], ca_js[, 1])
mica_js <- sapply(ca_js_list, function(x) {
x[which.max(IC[x])]
})
js <- as.numeric(names(ca_js_list)) + i - 1
res <- IC[mica_js]
sim[i, js] <- res
}
else {
for (j in i:num_terms) {
ancestor_j <- which(sCP[j, ] == 1)
ancestors <- intersect(ancestor_i, ancestor_j)
mica <- ancestors[which.max(IC[ancestors])]
res <- IC[mica]
sim[i, j] <- res
}
}
}
}
else if (method.term == "Lin") {
for (i in 1:num_terms) {
ancestor_i <- which(sCP[i, ] == 1)
progress_indicate(i, num_terms, 10, flag = T)
if (fast) {
mat <- sCP[i:num_terms, ]
ancestor_js <- which(matrix(mat == 1, nrow = num_terms -
i + 1), arr.ind = T)
flag <- is.element(ancestor_js[, 2], ancestor_i)
ca_js <- ancestor_js[flag, ]
ca_js_list <- split(ca_js[, 2], ca_js[, 1])
mica_js <- sapply(ca_js_list, function(x) {
x[which.max(IC[x])]
})
js <- as.numeric(names(ca_js_list)) + i - 1
res <- 2 * IC[mica_js]/(IC[ind[i]] + IC[ind[js]])
sim[i, js] <- res
}
else {
for (j in i:num_terms) {
ancestor_j <- which(sCP[j, ] == 1)
ancestors <- intersect(ancestor_i, ancestor_j)
mica <- ancestors[which.max(IC[ancestors])]
res <- 2 * IC[mica]/(IC[ind[i]] + IC[ind[j]])
sim[i, j] <- res
}
}
}
}
else if (method.term == "Schlicker") {
for (i in 1:num_terms) {
ancestor_i <- which(sCP[i, ] == 1)
progress_indicate(i, num_terms, 10, flag = T)
if (fast) {
mat <- sCP[i:num_terms, ]
ancestor_js <- which(matrix(mat == 1, nrow = num_terms -
i + 1), arr.ind = T)
flag <- is.element(ancestor_js[, 2], ancestor_i)
ca_js <- ancestor_js[flag, ]
ca_js_list <- split(ca_js[, 2], ca_js[, 1])
mica_js <- sapply(ca_js_list, function(x) {
x[which.max(IC[x])]
})
js <- as.numeric(names(ca_js_list)) + i - 1
res <- (2 * IC[mica_js]/(IC[ind[i]] + IC[ind[js]])) *
(1 - 10^(-IC[mica_js]))
sim[i, js] <- res
}
else {
for (j in i:num_terms) {
ancestor_j <- which(sCP[j, ] == 1)
ancestors <- intersect(ancestor_i, ancestor_j)
mica <- ancestors[which.max(IC[ancestors])]
res <- (2 * IC[mica]/(IC[ind[i]] + IC[ind[j]])) *
(1 - 10^(-IC[mica]))
sim[i, j] <- res
}
}
}
}
else if (method.term == "Jiang") {
for (i in 1:num_terms) {
ancestor_i <- which(sCP[i, ] == 1)
progress_indicate(i, num_terms, 10, flag = T)
if (fast) {
mat <- sCP[i:num_terms, ]
ancestor_js <- which(matrix(mat == 1, nrow = num_terms -
i + 1), arr.ind = T)
flag <- is.element(ancestor_js[, 2], ancestor_i)
ca_js <- ancestor_js[flag, ]
ca_js_list <- split(ca_js[, 2], ca_js[, 1])
mica_js <- sapply(ca_js_list, function(x) {
x[which.max(IC[x])]
})
js <- as.numeric(names(ca_js_list)) + i - 1
tmp <- IC[ind[i]] + IC[ind[js]] - 2 * IC[mica_js]
tmp[tmp > 1] <- 1
res <- 1 - tmp
sim[i, js] <- res
}
else {
for (j in i:num_terms) {
ancestor_j <- which(sCP[j, ] == 1)
ancestors <- intersect(ancestor_i, ancestor_j)
mica <- ancestors[which.max(IC[ancestors])]
res <- 1 - min(1, IC[ind[i]] + IC[ind[j]] -
2 * IC[mica])
sim[i, j] <- res
}
}
}
}
else if (method.term == "Pesquita") {
for (i in 1:num_terms) {
ancestor_i <- which(sCP[i, ] == 1)
progress_indicate(i, num_terms, 10, flag = T)
if (fast) {
mat <- sCP[i:num_terms, ]
ancestor_js <- which(matrix(mat == 1, nrow = num_terms -
i + 1), arr.ind = T)
flag <- is.element(ancestor_js[, 2], ancestor_i)
ca_js <- ancestor_js[flag, ]
ca_js_list <- split(ca_js[, 2], ca_js[, 1])
mica_js <- sapply(ca_js_list, function(x) {
x[which.max(IC[x])]
})
js <- as.numeric(names(ca_js_list)) + i - 1
allan_js_list <- split(ancestor_js[, 2], ancestor_js[,
1])
allan_union <- sapply(allan_js_list, function(x) {
ux <- base::union(x, ancestor_i)
sum(IC[ux])
})
allca_union <- sapply(ca_js_list, function(x) {
sum(IC[x])
})
res <- allca_union/allan_union
sim[i, js] <- res
}
else {
for (j in i:num_terms) {
ancestor_j <- which(sCP[j, ] == 1)
ancestors <- intersect(ancestor_i, ancestor_j)
mica <- ancestors[which.max(IC[ancestors])]
allancestors <- base::union(ancestor_i, ancestor_j)
res <- sum(IC[ancestors])/sum(IC[allancestors])
sim[i, j] <- res
}
}
}
}
sim <- sim + Matrix::t(sim)
}
if (verbose) {
now <- Sys.time()
message(sprintf("Semantic similarity has been calculated (\%s)!",
as.character(now)), appendLF = T)
}
rownames(sim) <- colnames(sim) <- terms
sim[as.matrix(is.na(sim))] <- 0
if (class(sim) == "dgCMatrix") {
res <- xConverter(sim, from = "dgCMatrix", to = "igraph",
verbose = verbose)
}
invisible(res)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
