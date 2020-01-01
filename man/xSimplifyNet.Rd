\name{xSimplifyNet}
\alias{xSimplifyNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSimplifyNet(g, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{g}{
%% ~~Describe \code{g} here~~
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
function (g, verbose = TRUE)
{
if (class(g) == "igraph") {
ig <- g
if (verbose) {
now <- Sys.time()
message(sprintf("The input graph has \%d nodes and \%d edges (\%s)
...",
vcount(ig), ecount(ig), as.character(now)), appendLF = TRUE)
}
}
else {
stop("The function must apply to the 'igraph' object.\n")
}
if (!igraph::is_directed(ig)) {
stop("The function must apply to the direct graph.\n")
}
root <- dnet::dDAGroot(ig)
if (verbose) {
message(sprintf("\tthe number of roots: \%d", length(root)),
appendLF = TRUE)
}
tips <- dnet::dDAGtip(ig)
if (verbose) {
message(sprintf("\tthe number of tips: \%d", length(tips)),
appendLF = TRUE)
}
epaths <- igraph::get.shortest.paths(ig, from = root, to = tips,
output = "epath")
subg <- igraph::subgraph.edges(ig, eids = unlist(epaths$epath),
delete.vertices = TRUE)
if (verbose) {
now <- Sys.time()
message(sprintf("The simplified graph has \%d nodes and \%d edges (\%s)
...",
vcount(subg), ecount(subg), as.character(now)), appendLF = TRUE)
}
invisible(subg)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
