\name{xGRsort}
\alias{xGRsort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGRsort(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
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
function (data)
{
data <- gsub(",.*", "", data)
input <- do.call(rbind, strsplit(data, ":|-|,"))
if (ncol(input) == 2) {
data <- paste(input[, 1], ":", input[, 2], "-", input[,
2], sep = "")
}
suppressWarnings(gr <- xGR(data, format = "chr:start-end"))
if (is.null(gr)) {
return(NULL)
}
gr <- GenomeInfoDb::sortSeqlevels(gr)
gr <- sort(gr)
df <- data.frame(sid = 1:length(gr), gr = names(gr), stringsAsFactors =
F)
ind <- match(data, df$gr)
df <- df[ind, ]
df$oid <- 1:length(data)
sid <- NULL
df <- df \%>\% dplyr::arrange(sid)
invisible(df$oid)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
