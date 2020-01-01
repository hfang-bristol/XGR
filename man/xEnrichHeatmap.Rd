\name{xEnrichHeatmap}
\alias{xEnrichHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichHeatmap(list_eTerm, fdr.cutoff = 0.05, displayBy = c("zscore",
"fdr", "pvalue", "fc", "or"), colormap = NULL, zlim = NULL, reorder =
c("none", "row", "col", "both"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{list_eTerm}{
%% ~~Describe \code{list_eTerm} here~~
}
\item{fdr.cutoff}{
%% ~~Describe \code{fdr.cutoff} here~~
}
\item{displayBy}{
%% ~~Describe \code{displayBy} here~~
}
\item{colormap}{
%% ~~Describe \code{colormap} here~~
}
\item{zlim}{
%% ~~Describe \code{zlim} here~~
}
\item{reorder}{
%% ~~Describe \code{reorder} here~~
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
function (list_eTerm, fdr.cutoff = 0.05, displayBy = c("zscore",
"fdr", "pvalue", "fc", "or"), colormap = NULL, zlim = NULL,
reorder = c("none", "row", "col", "both"))
{
displayBy <- match.arg(displayBy)
reorder <- match.arg(reorder)
if (is.null(list_eTerm)) {
warnings("There is no enrichment in the 'eTerm' object.\n")
return(NULL)
}
if (class(list_eTerm) == "ls_eTerm" | class(list_eTerm) ==
"data.frame") {
if (class(list_eTerm) == "ls_eTerm") {
df <- list_eTerm$df
}
else if (class(list_eTerm) == "data.frame") {
if (displayBy == "fdr") {
if (is.null(list_eTerm$fdr)) {
list_eTerm$fdr <- list_eTerm$adjp
}
}
if (is.null(list_eTerm$id) & !is.null(list_eTerm$name)) {
list_eTerm$id <- list_eTerm$name
}
if (all(c("group", "ontology", "id", "name", "adjp",
displayBy) \%in\% colnames(list_eTerm))) {
df <- list_eTerm[, c("group", "ontology", "id",
"name", "adjp", displayBy)]
}
else if (all(c("group", "id", "name", "adjp", displayBy) \%in\%
colnames(list_eTerm))) {
df <- list_eTerm[, c("group", "id", "name", "adjp",
displayBy)]
df$ontology <- "ontology"
}
else if (all(c("ontology", "id", "name", "adjp",
displayBy) \%in\% colnames(list_eTerm))) {
df <- list_eTerm[, c("ontology", "id", "name",
"adjp", displayBy)]
df$group <- "group"
}
else if (all(c("id", "name", "adjp", displayBy) \%in\%
colnames(list_eTerm))) {
df <- list_eTerm[, c("id", "name", "adjp", displayBy)]
df$group <- "group"
df$ontology <- "ontology"
}
else {
warnings("The input data.frame does not contain required columns:
c('group','ontology','name','adjp').\n")
return(NULL)
}
}
df_all <- df
if (!is.null(df_all)) {
adjp <- NULL
gp <- NULL
mat <- NULL
ls_df <- split(x = df_all, f = df_all$ontology)
ls_df <- ls_df[unique(df_all$ontology)]
ls_mat <- lapply(1:length(ls_df), function(i) {
df <- ls_df[[i]]
ind <- which(df$adjp < fdr.cutoff)
if (length(ind) >= 1) {
df <- as.data.frame(df \%>\% dplyr::filter(adjp <
fdr.cutoff))
if (displayBy == "fdr") {
mat <- as.matrix(xSparseMatrix(df[, c("name",
"group", "adjp")], rows = unique(df$name),
columns = unique(df_all$group)))
mat[mat == 0] <- NA
mat <- -log10(mat)
}
else if (displayBy == "pvalue") {
mat <- as.matrix(xSparseMatrix(df[, c("name",
"group", "pvalue")], rows = unique(df$name),
columns = unique(df_all$group)))
mat[mat == 0] <- NA
mat <- -log10(mat)
}
else if (displayBy == "zscore") {
mat <- as.matrix(xSparseMatrix(df[, c("name",
"group", "zscore")], rows = unique(df$name),
columns = unique(df_all$group)))
mat[mat == 0] <- NA
}
else if (displayBy == "fc") {
mat <- as.matrix(xSparseMatrix(df[, c("name",
"group", "fc")], rows = unique(df$name),
columns = unique(df_all$group)))
mat[mat == 0] <- NA
mat <- log2(mat)
}
else if (displayBy == "or") {
mat <- as.matrix(xSparseMatrix(df[, c("name",
"group", "or")], rows = unique(df$name),
columns = unique(df_all$group)))
mat[mat == 0] <- NA
mat <- log2(mat)
}
if (nrow(mat) == 1) {
df_mat <- mat
}
else {
rname_ordered <- rownames(mat)[order(-nchar(rownames(mat)))]
if (names(ls_df)[i] == "PS2" || names(ls_df)[i] ==
"PSG") {
df_tmp <- unique(df[, c("id", "name")])
df_tmp <- df_tmp[with(df_tmp, order(as.numeric(df_tmp$id))),
]
rname_ordered <- df_tmp$name
}
ind <- match(rname_ordered, rownames(mat))
df_mat <- as.matrix(mat[ind, ], ncol = ncol(mat))
colnames(df_mat) <- colnames(mat)
colnames(df_mat) <- colnames(mat)
}
return(df_mat)
}
else {
return(NULL)
}
})
mat <- do.call(rbind, ls_mat)
if (!is.null(mat)) {
if (displayBy == "fdr" | displayBy == "pvalue") {
if (is.null(colormap)) {
colormap <- "grey100-darkorange"
}
if (is.null(zlim)) {
zlim <- c(0, ceiling(max(mat[!is.na(mat)])))
}
if (displayBy == "fdr") {
legend.title <- expression(-log[10]("FDR"))
}
else if (displayBy == "pvalue") {
legend.title <- expression(-log[10]("p-value"))
}
}
else if (displayBy == "fc" | displayBy == "zscore" |
displayBy == "or") {
tmp_max <- ceiling(max(mat[!is.na(mat)]))
tmp_min <- floor(min(mat[!is.na(mat)]))
if (tmp_max > 0 & tmp_min < 0) {
if (is.null(colormap)) {
colormap <- "deepskyblue-grey100-darkorange"
}
if (is.null(zlim)) {
tmp <- max(tmp_max, abs(tmp_min))
zlim <- c(-tmp, tmp)
}
}
else if (tmp_max <= 0) {
if (is.null(colormap)) {
colormap <- "deepskyblue-grey100"
}
if (is.null(zlim)) {
zlim <- c(tmp_min, 0)
}
}
else if (tmp_min >= 0) {
if (is.null(colormap)) {
colormap <- "grey100-darkorange"
}
if (is.null(zlim)) {
zlim <- c(0, tmp_max)
}
}
if (displayBy == "fc") {
legend.title <- expression(log[2]("FC"))
}
else if (displayBy == "zscore") {
legend.title <- ("Z-score")
}
else if (displayBy == "or") {
legend.title <- expression(log[2]("OR"))
}
}
gp <- xHeatmap(mat, reorder = reorder, colormap = colormap,
ncolors = 64, zlim = zlim, legend.title = legend.title,
barwidth = 0.4, x.rotate = 60, shape = 19,
size = 2, x.text.size = 6, y.text.size = 6,
na.color = "transparent", barheight = max(3,
min(5, nrow(mat))))
gp <- gp + theme(legend.title = element_text(size = 8))
gp$mat <- mat
}
}
else {
mat <- NULL
gp <- NULL
}
}
invisible(gp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
