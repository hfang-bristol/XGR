\name{xGeneID2Symbol}
\alias{xGeneID2Symbol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGeneID2Symbol(data, org = c("human", "mouse"), details = F, verbose =
T, RData.location = "http://galahad.well.ox.ac.uk/bigdata", guid =
NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{org}{
%% ~~Describe \code{org} here~~
}
\item{details}{
%% ~~Describe \code{details} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
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
function (data, org = c("human", "mouse"), details = F, verbose = T,
RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
if (!is.vector(data)) {
stop("The input data must be a vector.\n")
}
GeneID <- as.numeric(data)
if (verbose) {
message(sprintf("\%d GeneID of input data (\%s)", length(data),
as.character(Sys.time())), appendLF = T)
}
if (class(org) == "EG") {
df_eg <- org$gene_info
if (verbose) {
message(sprintf("Customised organism (\%s)", as.character(Sys.time())),
appendLF = T)
}
}
else {
org <- org[1]
if (org == "human") {
df_eg <- xRDataLoader(RData.customised = "org.Hs.eg",
RData.location = RData.location, guid = guid,
verbose = verbose)$gene_info
}
else if (org == "mouse") {
df_eg <- xRDataLoader(RData.customised = "org.Mm.eg",
RData.location = RData.location, guid = guid,
verbose = verbose)$gene_info
}
if (verbose) {
message(sprintf("\%s organism (\%s)", org, as.character(Sys.time())),
appendLF = T)
}
}
ind <- match(data, df_eg$GeneID)
df_res <- df_eg[ind, ]
if (verbose) {
message(sprintf("\%d mappable but \%d left unmappable (\%s)",
sum(!is.na(ind)), sum(is.na(ind)), as.character(Sys.time())),
appendLF = T)
}
if (details) {
df_res <- data.frame(Input = GeneID, df_res, stringsAsFactors = F)
return(df_res)
}
else {
return(df_res$Symbol)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
