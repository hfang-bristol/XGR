\name{xEnrichConciser}
\alias{xEnrichConciser}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichConciser(eTerm, cutoff = c(0.9, 0.5), verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{eTerm}{
%% ~~Describe \code{eTerm} here~~
}
\item{cutoff}{
%% ~~Describe \code{cutoff} here~~
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
function (eTerm, cutoff = c(0.9, 0.5), verbose = T)
{
if (is.null(eTerm)) {
warnings("There is no enrichment in the 'eTerm' object.\n")
return(NULL)
}
if (class(eTerm) == "eTerm") {
cross <- eTerm$cross
df <- xEnrichViewer(eTerm, top_num = "all", sortBy = "pvalue")
ind <- match(rownames(df), rownames(cross))
cross <- cross[ind, ind]
nRedundant_1 <- matrix(0, nrow = ncol(cross), ncol = ncol(cross))
nRedundant_2 <- matrix(0, nrow = ncol(cross), ncol = ncol(cross))
for (j in seq(2, ncol(cross))) {
for (i in seq(1, j - 1)) {
nRedundant_1[i, j] <- cross[i, j] >= cross[j,
j] * cutoff[1]
nRedundant_2[i, j] <- cross[i, j] >= cross[i,
i] * cutoff[2]
}
}
nRedundant <- apply(nRedundant_1 & nRedundant_2, 2, sum)
names(nRedundant) <- colnames(cross)
ind <- match(names(eTerm$adjp), names(nRedundant))
nRedundant <- nRedundant[ind]
flag <- nRedundant == 0
if (verbose) {
now <- Sys.time()
message(sprintf("\tAmong \%d terms, there are \%d non-redundant terms",
length(flag), sum(flag)), appendLF = T)
}
eTerm$term_info <- eTerm$term_info[flag, ]
eTerm$annotation <- eTerm$annotation[flag]
eTerm$overlap <- eTerm$overlap[flag]
eTerm$fc <- eTerm$fc[flag]
eTerm$zscore <- eTerm$zscore[flag]
eTerm$pvalue <- eTerm$pvalue[flag]
eTerm$adjp <- eTerm$adjp[flag]
eTerm$or <- eTerm$or[flag]
eTerm$CIl <- eTerm$CIl[flag]
eTerm$CIu <- eTerm$CIu[flag]
eTerm$cross <- eTerm$cross[flag, flag]
res <- eTerm
}
else {
res <- NULL
message("The provided object cannot be used for removing redundant
terms.\n")
}
res
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
