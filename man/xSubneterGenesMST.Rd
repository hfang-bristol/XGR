\name{xSubneterGenesMST}
\alias{xSubneterGenesMST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSubneterGenesMST(isubg, metric = c("hybrid", "max", "jaccard"),
verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{isubg}{
%% ~~Describe \code{isubg} here~~
}
\item{metric}{
%% ~~Describe \code{metric} here~~
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
function (isubg, metric = c("hybrid", "max", "jaccard"), verbose = T)
{
metric <- match.arg(metric)
if (class(isubg) == "iSubg") {
ls_subg <- isubg$ls_subg
g <- isubg$g
ls_subg <- base::Filter(base::Negate(is.null), ls_subg)
if (length(ls_subg) <= 1) {
return(NULL)
}
}
else {
stop("The function must apply to a 'iSubg' objects.\n")
}
ls_df <- lapply(1:length(ls_subg), function(i) {
x <- ls_subg[[i]]
y <- range(as.numeric(V(x)$significance))
data.frame(iteration = i, num.node = vcount(x), num.edge = ecount(x),
num.linker = sum(V(x)$type == "linker"), threshold = x$threshold,
sig.min = y[1], sig.max = y[2], stringsAsFactors = F)
})
df_summary <- do.call(rbind, ls_df)
ls_df <- lapply(1:length(ls_subg), function(i) {
x <- ls_subg[[i]]
y <- as_data_frame(x, what = "vertices")
data.frame(iteration = i, y, stringsAsFactors = F)
})
df_detail <- do.call(rbind, ls_df)
name <- iteration <- score <- NULL
mat <- df_detail \%>\% dplyr::select(name, iteration, score) \%>\%
tidyr::spread(key = iteration, value = score)
rownames(mat) <- mat[, 1]
mat <- mat[, -1]
mat[is.na(mat)] <- 0
mat[mat < 0] <- -1
mat[mat > 0] <- 1
mat <- data.matrix(mat)
n <- ncol(mat)
mat_max <- matrix(0, nrow = n, ncol = n)
mat_jaccard <- matrix(0, nrow = n, ncol = n)
for (i in 1:(n - 1)) {
name_i <- names(which(mat[, i] != 0))
v <- match(name_i, V(g)$name)
for (j in (i + 1):n) {
name_j <- names(which(mat[, j] != 0))
to <- match(name_j, V(g)$name)
dist_tmp <- igraph::distances(g, v = v, to = to,
mode = "all")
mat_max[i, j] <- max(dist_tmp)
mat_jaccard[i, j] <- 1 - length(intersect(v, to))/length(union(v,
to))
}
}
mat_max <- mat_max + t(mat_max)
colnames(mat_max) <- rownames(mat_max) <- colnames(mat)
mat_jaccard <- mat_jaccard + t(mat_jaccard)
colnames(mat_jaccard) <- rownames(mat_jaccard) <- colnames(mat)
if (metric == "max") {
mat_dist <- mat_max
}
else if (metric == "jaccard") {
mat_dist <- mat_jaccard
}
else if (metric == "hybrid") {
mat_dist <- mat_max * mat_jaccard
}
g_subg <- xConverter(mat_dist, from = "dgCMatrix", to = "igraph",
verbose = verbose)
mst <- igraph::minimum.spanning.tree(g_subg, weights =
E(g_subg)$weight)
w <- E(mst)$weight
E(mst)$weight_scaled <- 1 + 4 * (w - min(w))/(max(w) - min(w))
set.seed(825)
coords <- igraph::layout_with_fr(mst, weights = 5 -
(E(mst)$weight_scaled -
1))
V(mst)$xcoord <- coords[, 1]
V(mst)$ycoord <- coords[, 2]
ind <- match(V(mst)$name, df_summary$iteration)
V(mst)$num_nodes <- df_summary$num.node[ind]
V(mst)$num_edges <- df_summary$num.edge[ind]
mst$summary <- df_summary
mst$detail <- df_detail
mst$matrix <- mat
return(mst)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
