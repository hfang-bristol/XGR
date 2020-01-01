\name{xEnrichViewer}
\alias{xEnrichViewer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichViewer(eTerm, top_num = 10, sortBy = c("adjp", "fdr", "pvalue",
"zscore", "fc", "nAnno", "nOverlap", "or", "none"), decreasing = NULL,
details = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{eTerm}{
%% ~~Describe \code{eTerm} here~~
}
\item{top_num}{
%% ~~Describe \code{top_num} here~~
}
\item{sortBy}{
%% ~~Describe \code{sortBy} here~~
}
\item{decreasing}{
%% ~~Describe \code{decreasing} here~~
}
\item{details}{
%% ~~Describe \code{details} here~~
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
function (eTerm, top_num = 10, sortBy = c("adjp", "fdr", "pvalue",
"zscore", "fc", "nAnno", "nOverlap", "or", "none"), decreasing = NULL,
details = F)
{
sortBy <- match.arg(sortBy)
if (sortBy == "fdr") {
sortBy <- "adjp"
}
if (is.null(eTerm)) {
warnings("There is no enrichment in the 'eTerm' object.\n")
return(NULL)
}
if (class(eTerm) == "data.frame") {
warnings("The function apply to a 'data.frame' object.\n")
if (is.null(top_num)) {
top_num <- nrow(eTerm)
}
if (top_num > nrow(eTerm)) {
top_num <- nrow(eTerm)
}
top_num <- as.integer(top_num)
tab <- eTerm
}
else if (class(eTerm) == "eTerm") {
if (is.null(top_num)) {
top_num <- length(eTerm$term_info$id)
}
if (top_num > length(eTerm$term_info$id)) {
top_num <- length(eTerm$term_info$id)
}
top_num <- as.integer(top_num)
if (dim(eTerm$term_info)[1] == 1) {
tab <- data.frame(name = as.character(eTerm$term_info$name),
nAnno = as.numeric(sapply(eTerm$annotation, length)),
nOverlap = as.numeric(sapply(eTerm$overlap, length)),
fc = as.numeric(eTerm$fc), zscore = as.numeric(eTerm$zscore),
pvalue = as.numeric(eTerm$pvalue), adjp = as.numeric(eTerm$adjp),
or = as.numeric(eTerm$or), CIl = as.numeric(eTerm$CIl),
CIu = as.numeric(eTerm$CIu), distance = eTerm$term_info$distance,
namespace = eTerm$term_info$namespace, members_Overlap =
sapply(eTerm$overlap,
function(x) paste(sort(x), collapse = ", ")),
members_Anno = sapply(eTerm$annotation, function(x) paste(sort(x),
collapse = ", ")), stringsAsFactors = F)
}
else {
tab <- data.frame(name = as.character(eTerm$term_info$name),
nAnno = as.numeric(sapply(eTerm$annotation, length)),
nOverlap = as.numeric(sapply(eTerm$overlap, length)),
fc = as.numeric(eTerm$fc), zscore = as.numeric(eTerm$zscore),
pvalue = as.numeric(eTerm$pvalue), adjp = as.numeric(eTerm$adjp),
or = as.numeric(eTerm$or), CIl = as.numeric(eTerm$CIl),
CIu = as.numeric(eTerm$CIu), distance = eTerm$term_info$distance,
namespace = eTerm$term_info$namespace, members_Overlap =
sapply(eTerm$overlap,
function(x) paste(sort(x), collapse = ", ")),
members_Anno = sapply(eTerm$annotation, function(x) paste(sort(x),
collapse = ", ")), stringsAsFactors = F)
}
rownames(tab) <- eTerm$term_info$id
}
if (details == T) {
res <- tab[, c(1:14)]
}
else {
res <- tab[, c(1:10)]
}
if (is.null(decreasing)) {
if (sortBy == "zscore" | sortBy == "nAnno" | sortBy ==
"nOverlap" | sortBy == "or") {
decreasing <- T
}
else {
decreasing <- F
}
}
switch(sortBy, fdr = {
res <- res[with(res, order(adjp, -zscore))[1:top_num],
]
}, adjp = {
res <- res[with(res, order(adjp, -zscore))[1:top_num],
]
}, pvalue = {
res <- res[with(res, order(pvalue, -zscore))[1:top_num],
]
}, zscore = {
res <- res[with(res, order(-zscore, adjp))[1:top_num],
]
}, fc = {
res <- res[with(res, order(-fc, adjp))[1:top_num], ]
}, or = {
res <- res[with(res, order(-or, adjp))[1:top_num], ]
}, nAnno = {
res <- res[with(res, order(-nAnno, adjp, -zscore))[1:top_num],
]
}, nOverlap = {
res <- res[with(res, order(-nOverlap, adjp, -zscore))[1:top_num],
]
}, none = {
res <- res[order(rownames(res), decreasing = decreasing)[1:top_num],
]
})
if (sortBy == "none") {
suppressWarnings(flag <- all(!is.na(as.numeric(rownames(res)))))
if (flag) {
res <- res[order(as.numeric(rownames(res)), decreasing =
decreasing)[1:top_num],
]
}
}
res
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
