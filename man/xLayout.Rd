\name{xLayout}
\alias{xLayout}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xLayout(g, layout = c("layout_nicely", "layout_randomly",
"layout_in_circle", "layout_on_sphere", "layout_with_fr",
"layout_with_kk", "layout_as_tree", "layout_with_lgl",
"layout_with_graphopt", "layout_with_sugiyama", "layout_with_dh",
"layout_with_drl", "layout_with_gem", "layout_with_mds",
"layout_as_bipartite", "gplot.layout.adj", "gplot.layout.circle",
"gplot.layout.circrand", "gplot.layout.eigen",
"gplot.layout.fruchtermanreingold", "gplot.layout.geodist",
"gplot.layout.hall", "gplot.layout.kamadakawai", "gplot.layout.mds",
"gplot.layout.princoord", "gplot.layout.random", "gplot.layout.rmds",
"gplot.layout.segeo", "gplot.layout.seham", "gplot.layout.spring",
"gplot.layout.springrepulse", "gplot.layout.target",
"graphlayouts.layout_with_stress", "graphlayouts.layout_as_backbone",
"gephi.forceatlas2"), seed = 825, flip = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{g}{
%% ~~Describe \code{g} here~~
}
\item{layout}{
%% ~~Describe \code{layout} here~~
}
\item{seed}{
%% ~~Describe \code{seed} here~~
}
\item{flip}{
%% ~~Describe \code{flip} here~~
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
function (g, layout = c("layout_nicely", "layout_randomly",
"layout_in_circle",
"layout_on_sphere", "layout_with_fr", "layout_with_kk",
"layout_as_tree",
"layout_with_lgl", "layout_with_graphopt", "layout_with_sugiyama",
"layout_with_dh", "layout_with_drl", "layout_with_gem",
"layout_with_mds",
"layout_as_bipartite", "gplot.layout.adj", "gplot.layout.circle",
"gplot.layout.circrand", "gplot.layout.eigen",
"gplot.layout.fruchtermanreingold",
"gplot.layout.geodist", "gplot.layout.hall",
"gplot.layout.kamadakawai",
"gplot.layout.mds", "gplot.layout.princoord", "gplot.layout.random",
"gplot.layout.rmds", "gplot.layout.segeo", "gplot.layout.seham",
"gplot.layout.spring", "gplot.layout.springrepulse",
"gplot.layout.target",
"graphlayouts.layout_with_stress", "graphlayouts.layout_as_backbone",
"gephi.forceatlas2"), seed = 825, flip = F)
{
layout <- layout[1]
if (class(g) == "graphNEL") {
ig <- igraph.from.graphNEL(g)
}
else {
ig <- g
}
if (class(ig) != "igraph") {
stop("The function must apply to either 'igraph' or 'graphNEL'
object.\n")
}
glayout <- NULL
if (grepl("gplot", layout)) {
m <- as.matrix(xConverter(ig, from = "igraph", to = "dgCMatrix",
verbose = F))
set.seed(seed)
eval(parse(text = paste0("glayout <- sna::", layout,
"(m, NULL)")))
}
else if (grepl("graphlayouts", layout)) {
layout <- gsub("graphlayouts.", "", layout)
set.seed(seed)
eval(parse(text = paste0("glayout <- graphlayouts::",
layout, "(ig)")))
if (layout == "layout_as_backbone") {
glayout <- glayout$xy
}
}
else if (layout == "gephi.forceatlas2") {
iterations = 100
linlog = FALSE
pos = NULL
nohubs = FALSE
k = 400
gravity = 1
ks = 0.1
ksmax = 10
delta = 1
center = NULL
tolerance = 0.1
dim = 2
if (length(names(igraph::edge.attributes(ig))) == 0) {
attr_g <- NULL
}
else {
attr_g <- names(igraph::edge.attributes(ig))
}
A <- igraph::get.adjacency(ig, type = "both", attr = attr_g,
edges = F, names = T, sparse = F)
if (is.null(center))
center <- rep(0, dim)
nnodes <- nrow(A)
Binary <- A
Binary[Binary != 0] <- 1
Deg <- rowSums(Binary)
Forces1 <- matrix(0, nrow = dim, ncol = nnodes)
if (is.null(pos)) {
difference <- 2000/(nnodes * dim)
set.seed(seed)
position <- matrix(sample(seq(-1000, 1000, difference),
nnodes * dim), nnodes, dim)
}
else {
position <- pos
}
temp <- which(position[, 1] == center[1])
for (index in 2:ncol(position)) {
temp <- intersect(temp, which(position[, index] ==
center[index]))
}
position[temp, ] <- center + 0.01
displacement <- matrix(0, nrow = dim, ncol = nnodes)
m <- nrow(position)
for (iteration in 1:iterations) {
displacement <- displacement * 0
Forces2 <- Forces1
Forces1 <- matrix(, nrow = dim, ncol = 0)
distances <- as.matrix(stats::dist(position))
distances[which(distances < 0.01)] <- 0.01
mylist <- vector("list", dim)
for (j in 1:dim) {
mylist[[j]] <- (tcrossprod(position[, j], rep(1,
m)) - tcrossprod(rep(1, m), position[, j]))/distances
}
Fr <- k * ((tcrossprod(rep(1, m), Deg) + 1) * (tcrossprod(Deg,
rep(1, m)) + 1))/distances
Fa <- distances
if (linlog) {
Fa <- log(1 + Fa)
}
Fa <- (A^delta) * Fa
if (nohubs) {
Fa <- Fa/(tcrossprod(Deg, rep(1, m)) + 1)
}
Farfunction <- function(x) rowSums(x * (Fr - Fa),
na.rm = T)
Far <- do.call(rbind, lapply(mylist, Farfunction))
uv2 <- apply(matrix(rep(center, m), nrow = m, byrow = T) -
position, 1, function(x) x/sqrt(sum(x^2)))
Fg <- uv2 * matrix(rep(gravity * (Deg + 1), dim),
nrow = dim, byrow = T)
Forces1 <- Far + Fg
Forces1 <- round(Forces1, 2)
swing <- abs(colSums((Forces1 - Forces2)^2)^(1/2))
Global_swing <- sum((Deg + 1) * swing)
if (all(swing == 0)) {
message(sprintf("Convergence reached at step \%d (\%s)",
iteration, as.character(Sys.time())), appendLF = T)
break
}
tra <- abs(colSums((Forces1 + Forces2)^2)^(1/2))/2
Global_tra <- sum((Deg + 1) * tra)
Global_speed <- tolerance * Global_tra/Global_swing
speed <- ks * Global_speed/(1 + Global_speed * (swing)^(1/2))
speed_constrain <- ksmax/abs(colSums((Forces1^2))^(1/2))
speed <- ifelse(speed >= speed_constrain, speed_constrain,
speed)
displacement <- Forces1 * t(matrix(rep(speed, dim),
nnodes, dim))
position <- position + t(displacement)
}
glayout <- position
}
else {
set.seed(seed)
eval(parse(text = paste0("glayout <- ", layout, "(ig)")))
}
if (!is.null(glayout)) {
if (flip) {
node.xcoord <- glayout[, 2]
node.ycoord <- glayout[, 1]
}
else {
node.xcoord <- glayout[, 1]
node.ycoord <- glayout[, 2]
}
if (max(node.xcoord) != min(node.xcoord)) {
node.xcoord <- (node.xcoord - min(node.xcoord))/(max(node.xcoord) -
min(node.xcoord)) * 2 - 1
}
if (max(node.ycoord) != min(node.ycoord)) {
node.ycoord <- (node.ycoord - min(node.ycoord))/(max(node.ycoord) -
min(node.ycoord)) * 2 - 1
}
glayout <- cbind(node.xcoord, node.ycoord)
V(ig)$xcoord <- glayout[, 1]
V(ig)$ycoord <- glayout[, 2]
}
invisible(ig)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
