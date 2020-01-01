\name{xEnrichLadder}
\alias{xEnrichLadder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichLadder(eTerm, sortBy = c("fdr", "or", "adjp", "pvalue",
"zscore", "fc", "nAnno", "nOverlap", "none"), top_num = 5, FDR.cutoff =
0.05, CI.one = T, colormap = "lightgrey-grey-black", x.rotate = 90,
x.text.size = 6, y.text.size = 6, shape = 22, size = 2, label =
c("concise", "full"), verbose = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{eTerm}{
%% ~~Describe \code{eTerm} here~~
}
\item{sortBy}{
%% ~~Describe \code{sortBy} here~~
}
\item{top_num}{
%% ~~Describe \code{top_num} here~~
}
\item{FDR.cutoff}{
%% ~~Describe \code{FDR.cutoff} here~~
}
\item{CI.one}{
%% ~~Describe \code{CI.one} here~~
}
\item{colormap}{
%% ~~Describe \code{colormap} here~~
}
\item{x.rotate}{
%% ~~Describe \code{x.rotate} here~~
}
\item{x.text.size}{
%% ~~Describe \code{x.text.size} here~~
}
\item{y.text.size}{
%% ~~Describe \code{y.text.size} here~~
}
\item{shape}{
%% ~~Describe \code{shape} here~~
}
\item{size}{
%% ~~Describe \code{size} here~~
}
\item{label}{
%% ~~Describe \code{label} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
}
\item{\dots}{
%% ~~Describe \code{\dots} here~~
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
function (eTerm, sortBy = c("fdr", "or", "adjp", "pvalue", "zscore",
"fc", "nAnno", "nOverlap", "none"), top_num = 5, FDR.cutoff = 0.05,
CI.one = T, colormap = "lightgrey-grey-black", x.rotate = 90,
x.text.size = 6, y.text.size = 6, shape = 22, size = 2, label =
c("concise",
"full"), verbose = T, ...)
{
sortBy <- match.arg(sortBy)
label <- match.arg(label)
gp_heatmap <- NULL
if (class(eTerm) == "eTerm") {
df <- xEnrichViewer(eTerm, top_num = "all")
if (!CI.one) {
ind <- which(df$CIl > 1 | df$CIu < 1)
df <- df[ind, ]
}
if (FDR.cutoff == 1) {
FDR.cutoff <- FDR.cutoff + 0.01
}
if (top_num == "auto") {
top_num <- sum(df$adjp < FDR.cutoff)
if (top_num <= 1) {
top_num <- 10
}
}
df_enrichment <- xEnrichViewer(eTerm, top_num = top_num,
sortBy = sortBy, details = T)
}
else if (class(eTerm) == "data.frame") {
df_enrichment <- eTerm
}
if (!is.null(df_enrichment)) {
if (label == "concise") {
df_enrichment$label <- paste0(df_enrichment$name,
" [OR=", df_enrichment$or, ", FDR=", df_enrichment$adjp,
", n=", df_enrichment$nOverlap, "/", df_enrichment$nAnno,
"]")
}
else {
df_enrichment$label <- paste0(df_enrichment$name,
"\n[OR=", df_enrichment$or, ", P=", df_enrichment$pvalue,
", FDR=", df_enrichment$adjp, ", n=", df_enrichment$nOverlap,
"/", df_enrichment$nAnno, "]")
}
df_enrichment <- df_enrichment[!duplicated(df_enrichment$name),
]
ls_path <- lapply(1:nrow(df_enrichment), function(j) {
x <- df_enrichment$members_Overlap[j]
query <- unlist(strsplit(x, ", "))
})
names(ls_path) <- df_enrichment$name
if (length(ls_path) >= 1) {
all_genes <- unique(unlist(ls_path))
ls_vec <- lapply(1:length(ls_path), function(j) {
ind <- match(all_genes, ls_path[[j]])
vec <- rep(NA, length(all_genes))
vec[!is.na(ind)] <- names(ls_path)[j]
vec
})
df_res <- do.call(cbind, ls_vec)
colnames(df_res) <- names(ls_path)
rownames(df_res) <- all_genes
vec_sum <- apply(!is.na(df_res), 1, sum)
df_tmp <- data.frame(num = vec_sum, df_res, gene = rownames(df_res),
stringsAsFactors = F)
df_tmp <- df_tmp \%>\% dplyr::arrange_all()
df_tmp <- df_tmp[nrow(df_tmp):1, ]
colnames(df_tmp)[2:(ncol(df_tmp) - 1)] <- colnames(df_res)
vec_sum <- df_tmp$num
names(vec_sum) <- df_tmp$gene
df_res <- df_tmp \%>\% dplyr::select(2:(ncol(df_tmp) -
1))
rownames(df_res) <- df_tmp$gene
if (verbose) {
message(sprintf("heatmap of \%d rows X \%d columns (\%s) ...",
nrow(df_res), ncol(df_res), as.character(Sys.time())),
appendLF = TRUE)
}
if (1) {
df_heatmap <- 0 + !is.na(df_res)
df_heatmap[df_heatmap == 0] <- NA
for (i in 1:nrow(df_heatmap)) {
x <- df_heatmap[i, ]
df_heatmap[i, !is.na(x)] <- vec_sum[i]
}
mat_heatmap <- t(df_heatmap)
ind <- match(rownames(mat_heatmap), df_enrichment$name)
rownames(mat_heatmap) <- df_enrichment$label[ind]
gp_heatmap <- xHeatmap(mat_heatmap, reorder = "none",
colormap = colormap, zlim = c(0, max(mat_heatmap,
na.rm = T)), ncolors = 64, barwidth = 0.4,
x.rotate = x.rotate, x.text.size = x.text.size,
y.text.size = y.text.size, shape = shape, size = size,
na.color = "transparent", ...)
gp_heatmap <- gp_heatmap + theme(legend.title = element_text(size = 8),
legend.position = "none") + scale_y_discrete(position = "right")
colsep <- cumsum(table(vec_sum))
colsep <- length(vec_sum) - colsep[-length(colsep)]
gp_heatmap <- gp_heatmap + geom_vline(xintercept = colsep +
0.5, color = "grey90", size = 0.5)
gene <- sample <- val <- NULL
data_matrix <- gp_heatmap$data \%>\% dplyr::select(gene,
sample, val) \%>\% tidyr::spread(sample, val)
rownames(data_matrix) <- data_matrix$gene
data_matrix <- data_matrix[, -1]
gp_heatmap$matrix <- data_matrix
}
}
}
return(gp_heatmap)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
