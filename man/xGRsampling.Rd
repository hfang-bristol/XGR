\name{xGRsampling}
\alias{xGRsampling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xGRsampling(GR.data, GR.background, num.samples = 100, gap.max = 50000,
max.distance = NULL, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{GR.data}{
%% ~~Describe \code{GR.data} here~~
}
\item{GR.background}{
%% ~~Describe \code{GR.background} here~~
}
\item{num.samples}{
%% ~~Describe \code{num.samples} here~~
}
\item{gap.max}{
%% ~~Describe \code{gap.max} here~~
}
\item{max.distance}{
%% ~~Describe \code{max.distance} here~~
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
function (GR.data, GR.background, num.samples = 100, gap.max = 50000,
max.distance = NULL, verbose = T)
{
if (is.null(max.distance)) {
max.distance <- gap.max
}
else if (max.distance > gap.max) {
max.distance <- gap.max
}
if (class(GR.data) != "GRanges") {
stop("The function must apply to a 'GRanges' object for input data.\n")
}
if (class(GR.background) == "GRangesList") {
GR.background <- BiocGenerics::unlist(GR.background)
}
else if (class(GR.background) == "list") {
GR.background <-
BiocGenerics::unlist(GenomicRanges::GRangesList(GR.background))
}
if (class(GR.background) != "GRanges") {
stop("The function must apply to a 'GRanges' object for input
background.\n")
}
if (verbose) {
now <- Sys.time()
message(sprintf("First, get non-overlapping regions for both input data
and background (\%s) ...",
as.character(now)), appendLF = T)
}
dGR_reduced <- IRanges::reduce(GR.data)
bGR_reduced <- IRanges::reduce(GR.background)
if (verbose) {
now <- Sys.time()
message(sprintf("\tnon-overlapping regions: \%d for data, \%d for
background",
length(dGR_reduced), length(bGR_reduced)), appendLF = T)
}
if (verbose) {
now <- Sys.time()
message(sprintf("Second, keep only data regions that are within
background regions (\%s) ...",
as.character(now)), appendLF = T)
}
mergeOverlaps <- function(qGR, sGR, maxgap = -1L, minoverlap = 0L) {
hits <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query =
qGR,
subject = sGR, maxgap = maxgap, minoverlap = minoverlap,
type = "any", select = "all", ignore.strand = T)))
qhits <- qGR[hits[, 1]]
shits <- sGR[hits[, 2]]
oGR <- IRanges::pintersect(qhits, shits)
IRanges::reduce(oGR)
}
dGR_reduced <- mergeOverlaps(qGR = dGR_reduced, sGR = bGR_reduced,
maxgap = -1L, minoverlap = 0L)
if (verbose) {
now <- Sys.time()
message(sprintf("\t\%d within background", length(dGR_reduced)),
appendLF = T)
}
if (verbose) {
now <- Sys.time()
message(sprintf("Third, find background islands that contain data
regions (\%s) ...",
as.character(now)), appendLF = T)
}
hits <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query =
dGR_reduced,
subject = GR.background, maxgap = gap.max - 1, minoverlap = 0L,
type = "any", select = "all", ignore.strand = T)))
ind_data <- hits[, 1]
ind_background <- hits[, 2]
dt_ls <- split(x = ind_background, f = ind_data)
if (verbose) {
now <- Sys.time()
message(sprintf("\t\%d background islands",
length(unique(ind_background))),
appendLF = T)
}
if (is.null(max.distance)) {
if (verbose) {
now <- Sys.time()
message(sprintf("Fourth, define the sampling range (\%s) ...",
as.character(now)), appendLF = T)
}
}
else {
max.distance <- as.integer(max.distance)
if (verbose) {
now <- Sys.time()
message(sprintf("Fourth, define the sampling range within \%d bp
distance away from data regions (\%s) ...",
max.distance, as.character(now)), appendLF = T)
}
}
df_data <- GenomicRanges::as.data.frame(dGR_reduced, row.names = NULL)
df_background <- GenomicRanges::as.data.frame(GR.background,
row.names = NULL)
range_ls <- lapply(1:length(dt_ls), function(i) {
df_dt <- df_data[as.numeric(names(dt_ls)[i]), ]
dw <- df_dt$width
df_bg <- df_background[dt_ls[[i]], ]
res <- lapply(1:nrow(df_bg), function(j) {
bs <- df_bg$start[j]
be <- df_bg$end[j]
rs <- bs
re <- be - dw
if (!is.null(max.distance)) {
if (rs < df_dt$start - max.distance) {
rs <- df_dt$start - max.distance
}
if (re > df_dt$end + max.distance) {
re <- df_dt$end + max.distance
}
}
seq(rs, re)
})
unlist(res)
})
if (verbose) {
now <- Sys.time()
message(sprintf("Fifth, do '\%d' sampling for the starting points (\%s)
...",
num.samples, as.character(now)), appendLF = T)
}
res_ls <- lapply(range_ls, function(x) {
base::sample(x, num.samples, replace = T)
})
df_samples <- do.call(rbind, res_ls)
if (verbose) {
now <- Sys.time()
message(sprintf("Finally, construct GR for each sampling (\%s) ...",
as.character(now)), appendLF = T)
}
df_dt_all <- df_data[as.numeric(names(dt_ls)), ]
sGR_list <- lapply(1:ncol(df_samples), function(j) {
sGR <- GenomicRanges::GRanges(seqnames =
S4Vectors::Rle(df_dt_all$seqnames),
ranges = IRanges::IRanges(start = df_samples[, j],
end = df_samples[, j] + df_dt_all$width - 1),
strand = S4Vectors::Rle(df_dt_all$strand))
})
invisible(sGR_list)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
