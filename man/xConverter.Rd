\name{xConverter}
\alias{xConverter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xConverter(obj, from = c("dgCMatrix", "igraph"), to = c("igraph",
"dgCMatrix"), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{obj}{
%% ~~Describe \code{obj} here~~
}
\item{from}{
%% ~~Describe \code{from} here~~
}
\item{to}{
%% ~~Describe \code{to} here~~
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
function (obj, from = c("dgCMatrix", "igraph"), to = c("igraph",
"dgCMatrix"), verbose = TRUE)
{
from <- match.arg(from)
to <- match.arg(to)
if (class(obj) != from) {
}
if (from != "igraph" & to != "igraph") {
stop(sprintf("Conversion between '\%s' and '\%s' is not supported.\n",
from, to))
}
if (from == to) {
warnings(sprintf("Since the class '\%s' converted from is the same as
the class '\%s' converted to, it will return exactly what you
input.\n",
from, to))
return(obj)
}
if (from == "igraph") {
data <- igraph::get.data.frame(obj, what = "vertices")
if ("weight" \%in\% list.edge.attributes(obj)) {
E(obj)$weight <- as.numeric(E(obj)$weight)
objConverted <- igraph::as_adjacency_matrix(obj,
type = "both", attr = "weight", edges = F, names = T,
sparse = getIgraphOpt("sparsematrices"))
}
else {
objConverted <- igraph::as_adjacency_matrix(obj,
type = "both", attr = NULL, edges = F, names = T,
sparse = getIgraphOpt("sparsematrices"))
}
}
else if (from == "dgCMatrix") {
nodes <- data.frame(name = rownames(obj))
nodenames <- rownames(obj)
adjM <- obj
tmp <- which(as.matrix(adjM != 0), arr.ind = T)
if (from == "dgCMatrix") {
ind <- which(tmp[, 1] < tmp[, 2])
ttmp <- matrix(0, nrow = length(ind), ncol = 2)
ttmp[1:length(ind), ] <- tmp[ind, ]
tmp <- ttmp
}
weight_flag <- T
if (all(adjM[tmp] == 1)) {
weight_flag <- F
}
if (weight_flag) {
relations <- data.frame(from = nodenames[tmp[, 1]],
to = nodenames[tmp[, 2]], weight = adjM[tmp])
}
else {
relations <- data.frame(from = nodenames[tmp[, 1]],
to = nodenames[tmp[, 2]])
}
if (from == "dgCMatrix") {
objConverted <- igraph::graph_from_data_frame(d = relations,
directed = F, vertices = nodes)
}
}
if (verbose) {
message(sprintf("Your input object of class '\%s' has been converted
into an object of class '\%s'.",
from, to), appendLF = T)
}
return(objConverted)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
