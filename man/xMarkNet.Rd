\name{xMarkNet}
\alias{xMarkNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xMarkNet(ig1, ig2, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ig1}{
%% ~~Describe \code{ig1} here~~
}
\item{ig2}{
%% ~~Describe \code{ig2} here~~
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
function (ig1, ig2, verbose = TRUE)
{
if (!(any(class(ig1) \%in\% "igraph")) & !(any(class(ig2) \%in\%
"igraph"))) {
return(NULL)
}
flag_direct <- "undirect"
vec <- sapply(list(ig1, ig2), igraph::is_directed)
if (sum(vec) == length(vec)) {
flag_direct <- "direct"
}
mark <- name <- NULL
node1 <- ig1 \%>\% igraph::as_data_frame("vertices")
edge1 <- ig1 \%>\% igraph::as_data_frame("edges")
node2 <- ig2 \%>\% igraph::as_data_frame("vertices") \%>\%
dplyr::mutate(mark = 1) \%>\%
dplyr::select(name, mark)
edge2 <- ig2 \%>\% igraph::as_data_frame("edges") \%>\%
dplyr::mutate(mark = 1)
node1 <- node1 \%>\% dplyr::left_join(node2, by = c("name")) \%>\%
dplyr::mutate(mark = ifelse(is.na(mark), 0, mark))
edge1 <- edge1 \%>\% dplyr::left_join(edge2, by = c("from",
"to")) \%>\% dplyr::mutate(mark = ifelse(is.na(mark), 0,
mark))
if (flag_direct == "direct") {
ig <- igraph::graph_from_data_frame(d = edge1, directed = TRUE,
vertices = node1)
}
else {
ig <- igraph::graph_from_data_frame(d = edge1, directed = FALSE,
vertices = node1)
}
if (verbose) {
message(sprintf("The '\%s' network (\%d nodes and \%d edges) has \%d
nodes and \%d edges marked",
flag_direct, vcount(ig), ecount(ig), sum(V(ig)$mark),
sum(E(ig)$mark), appendLF = TRUE))
}
invisible(ig)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
