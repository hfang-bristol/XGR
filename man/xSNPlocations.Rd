\name{xSNPlocations}
\alias{xSNPlocations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSNPlocations(data, GR.SNP = c("dbSNP_GWAS", "dbSNP_Common",
"dbSNP_Single"), verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{GR.SNP}{
%% ~~Describe \code{GR.SNP} here~~
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
function (data, GR.SNP = c("dbSNP_GWAS", "dbSNP_Common",
"dbSNP_Single"),
verbose = T, RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
data <- gsub("_", ":", data, perl = T)
data <- gsub("imm:", "chr", data, perl = T)
data <- unique(data)
if (verbose) {
now <- Sys.time()
message(sprintf("Load positional information for SNPs (\%s) ...",
as.character(now)), appendLF = T)
}
if (class(GR.SNP) == "GRanges") {
pos_SNP <- GR.SNP
}
else {
pos_SNP <- xRDataLoader(GR.SNP[1], verbose = verbose,
RData.location = RData.location, guid = guid)
if (is.null(pos_SNP)) {
GR.SNP <- "dbSNP_GWAS"
if (verbose) {
message(sprintf("Instead, \%s will be used", GR.SNP),
appendLF = T)
}
pos_SNP <- xRDataLoader(GR.SNP, verbose = verbose,
RData.location = RData.location, guid = guid)
}
}
ind <- match(data, names(pos_SNP))
data_rest <- data[is.na(ind)]
ind <- ind[!is.na(ind)]
if (length(ind)) {
gr_SNP <- pos_SNP[ind, ]
GenomicRanges::mcols(gr_SNP) <- NULL
}
else {
gr_SNP <- NULL
}
ind <- grep("^chr\\w+:\\d+", data_rest, perl = T)
if (length(ind) > 0) {
data_rest <- data_rest[ind]
res_ls <- strsplit(data_rest, ":")
res_df <- do.call(rbind, res_ls)
res_gr <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(res_df[,
1]), ranges = IRanges::IRanges(start = as.numeric(res_df[,
2]), end = as.numeric(res_df[, 2]), names = data_rest),
strand = S4Vectors::Rle(rep("*", nrow(res_df))))
}
else {
res_gr <- NULL
}
if (!is.null(gr_SNP)) {
if (!is.null(res_gr)) {
gr_SNP <- c(gr_SNP, res_gr)
}
}
else {
if (!is.null(res_gr)) {
gr_SNP <- res_gr
}
}
if (verbose) {
now <- Sys.time()
message(sprintf("\tOut of \%d input SNPs, \%d SNPs have positional
info",
length(data), length(gr_SNP)), appendLF = T)
}
if (!is.null(gr_SNP)) {
tmp_df <- GenomicRanges::as.data.frame(gr_SNP, row.names = NULL)
mcols_df <- data.frame(variant_id = paste(tmp_df[, 1],
":", tmp_df[, 3], sep = ""), stringsAsFactors = F)
GenomicRanges::mcols(gr_SNP) <- mcols_df
}
invisible(gr_SNP)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
