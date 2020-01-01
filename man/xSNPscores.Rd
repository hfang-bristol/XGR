\name{xSNPscores}
\alias{xSNPscores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSNPscores(data, include.LD = NA, LD.customised = NULL, LD.r2 = 0.8,
significance.threshold = 5e-05, score.cap = 10, verbose = T,
RData.location = "http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{include.LD}{
%% ~~Describe \code{include.LD} here~~
}
\item{LD.customised}{
%% ~~Describe \code{LD.customised} here~~
}
\item{LD.r2}{
%% ~~Describe \code{LD.r2} here~~
}
\item{significance.threshold}{
%% ~~Describe \code{significance.threshold} here~~
}
\item{score.cap}{
%% ~~Describe \code{score.cap} here~~
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
function (data, include.LD = NA, LD.customised = NULL, LD.r2 = 0.8,
significance.threshold = 5e-05, score.cap = 10, verbose = T,
RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
if (is.null(data)) {
stop("The input data must be not NULL.\n")
}
else {
if (class(data) == "DataFrame") {
data <- S4Vectors::as.matrix(data)
}
if (is.vector(data)) {
if (length(data) > 1) {
if (is.null(names(data))) {
stop("The input data must have names with attached dbSNP ID.\n")
}
}
else {
data <- utils::read.delim(file = data, header = F,
row.names = NULL, stringsAsFactors = F)
}
}
if (is.vector(data)) {
pval <- data
}
else if (is.matrix(data) | is.data.frame(data)) {
data <- as.matrix(data)
data_list <- split(x = data[, 2], f = as.character(data[,
1]))
res_list <- lapply(data_list, function(x) {
x <- as.numeric(x)
x <- x[!is.na(x)]
if (length(x) > 0) {
min(x)
}
else {
NULL
}
})
pval <- unlist(res_list)
}
tmp <- min(pval[pval != 0])
pval[pval < tmp] <- tmp
}
tmp <- names(pval)
tmp <- gsub("_", ":", tmp, perl = T)
names(pval) <- gsub("imm:", "chr", tmp, perl = T)
Lead_Sig <- data.frame(SNP = names(pval), Sig = pval, row.names = NULL,
stringsAsFactors = F)
leads <- Lead_Sig[, 1]
sigs <- Lead_Sig[, 2]
if (verbose) {
now <- Sys.time()
message(sprintf("A total of \%d Lead SNPs are input",
length(leads)), appendLF = T)
}
if (LD.r2 >= 0.8 & LD.r2 <= 1) {
default.include.LD <- c("ACB", "AFR", "AMR", "ASW", "BEB",
"CDX", "CEU", "CHB", "CHS", "CLM", "EAS", "ESN",
"EUR", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU",
"JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL",
"PUR", "SAS", "STU", "TSI", "YRI")
ind <- match(default.include.LD, include.LD)
include.LD <- default.include.LD[!is.na(ind)]
}
else {
include.LD <- NULL
}
LLR <- NULL
if (length(include.LD) > 0 & is.null(LD.customised)) {
if (verbose) {
now <- Sys.time()
message(sprintf("Inclusion of LD SNPs is based on population (\%s) with
R2 >= \%f",
paste(include.LD, collapse = ","), LD.r2), appendLF = T)
}
GWAS_LD <- xRDataLoader("GWAS_LD", RData.location = RData.location,
guid = guid, verbose = verbose)
res_list <- lapply(include.LD, function(x) {
data_ld <- ""
eval(parse(text = paste("data_ld <- GWAS_LD$", x,
sep = "")))
ind <- match(rownames(data_ld), leads)
ind_lead <- which(!is.na(ind))
if (length(ind_lead) >= 2) {
ind_ld <- which(Matrix::colSums(data_ld[ind_lead,
] >= LD.r2) > 0)
sLL <- data_ld[ind_lead, ind_ld]
summ <- summary(sLL)
res <- data.frame(Lead = rownames(sLL)[summ$i],
LD = colnames(sLL)[summ$j], R2 = summ$x, stringsAsFactors = F)
}
else if (length(ind_lead) == 1) {
ind_ld <- which(data_ld[ind_lead, ] >= LD.r2)
sLL <- data_ld[ind_lead, ind_ld]
res <- data.frame(Lead = rep(rownames(data_ld)[ind_lead],
length(sLL)), LD = colnames(data_ld)[ind_ld],
R2 = sLL, stringsAsFactors = F)
}
else {
NULL
}
})
LLR <- do.call(rbind, res_list)
if (1) {
ImmunoBase_LD <- xRDataLoader("ImmunoBase_LD", RData.location =
RData.location,
guid = guid, verbose = verbose)
res_list <- lapply(include.LD, function(x) {
data_ld <- ""
eval(parse(text = paste("data_ld <- ImmunoBase_LD$",
x, sep = "")))
ind <- match(rownames(data_ld), leads)
ind_lead <- which(!is.na(ind))
if (length(ind_lead) >= 2) {
ind_ld <- which(Matrix::colSums(data_ld[ind_lead,
] >= LD.r2) > 0)
sLL <- data_ld[ind_lead, ind_ld]
summ <- summary(sLL)
res <- data.frame(Lead = rownames(sLL)[summ$i],
LD = colnames(sLL)[summ$j], R2 = summ$x,
stringsAsFactors = F)
}
else if (length(ind_lead) == 1) {
ind_ld <- which(data_ld[ind_lead, ] >= LD.r2)
sLL <- data_ld[ind_lead, ind_ld]
res <- data.frame(Lead = rep(rownames(data_ld)[ind_lead],
length(sLL)), LD = colnames(data_ld)[ind_ld],
R2 = sLL, stringsAsFactors = F)
}
else {
NULL
}
})
LLR_tmp <- do.call(rbind, res_list)
LLR <- rbind(LLR, LLR_tmp)
}
}
else if (!is.null(LD.customised)) {
if (is.vector(LD.customised)) {
LLR <- utils::read.delim(file = LD.customised, header = F,
row.names = NULL, stringsAsFactors = F)
}
else if (is.matrix(LD.customised) | is.data.frame(LD.customised)) {
LLR <- LD.customised
}
if (!is.null(LLR)) {
flag <- LLR[, 3] >= LD.r2
if (sum(flag) > 0) {
LLR <- LLR[LLR[, 3] >= LD.r2, ]
colnames(LLR) <- c("Lead", "LD", "R2")
if (verbose) {
now <- Sys.time()
message(sprintf("Inclusion of LD SNPs is based on customised data (\%d
Lead SNPs and \%d LD SNPs) with R2>=\%f",
length(unique(LLR[, 1])), length(unique(LLR[,
2])), LD.r2), appendLF = T)
}
}
else {
LLR <- NULL
}
}
}
if (!is.null(LLR)) {
ld_list <- split(x = LLR[, -2], f = LLR[, 2])
res_list <- lapply(ld_list, function(x) {
ind <- match(x$Lead, leads)
min(sigs[ind]^x$R2)
})
vec <- unlist(res_list)
LD_Sig <- data.frame(SNP = names(vec), Sig = vec, row.names = NULL,
stringsAsFactors = F)
df <- rbind(Lead_Sig, as.matrix(LD_Sig))
res_list <- split(x = df$Sig, f = df$SNP)
vec <- unlist(lapply(res_list, min))
SNP_Sig <- data.frame(SNP = names(vec), FDR = vec, row.names = NULL,
stringsAsFactors = F)
}
else {
if (verbose) {
now <- Sys.time()
message(sprintf("Do not include any LD SNPs"), appendLF = T)
}
SNP_Sig <- Lead_Sig
}
if (verbose) {
now <- Sys.time()
message(sprintf("A total of \%d Lead/LD SNPs are considered",
nrow(SNP_Sig)), appendLF = T)
}
pval <- as.numeric(SNP_Sig[, 2])
names(pval) <- SNP_Sig[, 1]
if (is.null(significance.threshold)) {
scores <- log10((1 - pval)/pval)
}
else {
scores <- log10((1 - pval)/pval) - log10((1 -
significance.threshold)/significance.threshold)
}
tmp_max <- max(scores[!is.infinite(scores)])
tmp_min <- min(scores[!is.infinite(scores)])
scores[scores > tmp_max] <- tmp_max
scores[scores < tmp_min] <- tmp_min
ind_remained <- which(scores > 0)
seeds.snps <- scores[ind_remained]
pval <- pval[ind_remained]
if (verbose) {
now <- Sys.time()
message(sprintf("A total of \%d Lead/LD SNPs are scored positively",
sum(seeds.snps > 0)), appendLF = T)
}
flag <- rep("Lead", length(pval))
ind <- match(names(pval), Lead_Sig$SNP)
flag[is.na(ind)] <- "LD"
df_SNP <- data.frame(SNP = names(pval), Score = seeds.snps,
Pval = pval, Flag = flag, row.names = NULL, stringsAsFactors = F)
if (!is.null(score.cap)) {
score.cap <- as.numeric(score.cap)
if (score.cap <= max(df_SNP$Score)) {
df_SNP$Score[df_SNP$Score >= score.cap] <- score.cap
if (verbose) {
now <- Sys.time()
message(sprintf("SNP score capped to the maximum score \%d.",
score.cap), appendLF = T)
}
}
}
df_SNP <- df_SNP[order(df_SNP$Flag, df_SNP$Score, -df_SNP$Pval,
df_SNP$SNP, decreasing = TRUE), ]
invisible(df_SNP)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
