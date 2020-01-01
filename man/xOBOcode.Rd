\name{xOBOcode}
\alias{xOBOcode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xOBOcode(g, node.level = "term_distance", node.level.value = 2,
node.label.size = 2, node.label.color = "darkblue", node.label.alpha =
0.8, node.label.padding = 0, node.label.arrow = 0.01, node.label.force
= 0, node.shape = 19, node.xcoord = NULL, node.ycoord = NULL,
node.color = NULL, node.color.title = NULL, colormap = "grey-grey",
ncolors = 64, zlim = NULL, node.size.range = 4, title = "", edge.size =
0.5, edge.color = "black", edge.color.alpha = 0.4, edge.curve = 0.1,
edge.arrow = 2, edge.arrow.gap = 0.02, node.table = "term_name",
node.table.wrap = 50, table.base.size = 7, table.row.space = 2,
table.nrow = 55, table.ncol = NULL, root.code = "RT")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{g}{
%% ~~Describe \code{g} here~~
}
\item{node.level}{
%% ~~Describe \code{node.level} here~~
}
\item{node.level.value}{
%% ~~Describe \code{node.level.value} here~~
}
\item{node.label.size}{
%% ~~Describe \code{node.label.size} here~~
}
\item{node.label.color}{
%% ~~Describe \code{node.label.color} here~~
}
\item{node.label.alpha}{
%% ~~Describe \code{node.label.alpha} here~~
}
\item{node.label.padding}{
%% ~~Describe \code{node.label.padding} here~~
}
\item{node.label.arrow}{
%% ~~Describe \code{node.label.arrow} here~~
}
\item{node.label.force}{
%% ~~Describe \code{node.label.force} here~~
}
\item{node.shape}{
%% ~~Describe \code{node.shape} here~~
}
\item{node.xcoord}{
%% ~~Describe \code{node.xcoord} here~~
}
\item{node.ycoord}{
%% ~~Describe \code{node.ycoord} here~~
}
\item{node.color}{
%% ~~Describe \code{node.color} here~~
}
\item{node.color.title}{
%% ~~Describe \code{node.color.title} here~~
}
\item{colormap}{
%% ~~Describe \code{colormap} here~~
}
\item{ncolors}{
%% ~~Describe \code{ncolors} here~~
}
\item{zlim}{
%% ~~Describe \code{zlim} here~~
}
\item{node.size.range}{
%% ~~Describe \code{node.size.range} here~~
}
\item{title}{
%% ~~Describe \code{title} here~~
}
\item{edge.size}{
%% ~~Describe \code{edge.size} here~~
}
\item{edge.color}{
%% ~~Describe \code{edge.color} here~~
}
\item{edge.color.alpha}{
%% ~~Describe \code{edge.color.alpha} here~~
}
\item{edge.curve}{
%% ~~Describe \code{edge.curve} here~~
}
\item{edge.arrow}{
%% ~~Describe \code{edge.arrow} here~~
}
\item{edge.arrow.gap}{
%% ~~Describe \code{edge.arrow.gap} here~~
}
\item{node.table}{
%% ~~Describe \code{node.table} here~~
}
\item{node.table.wrap}{
%% ~~Describe \code{node.table.wrap} here~~
}
\item{table.base.size}{
%% ~~Describe \code{table.base.size} here~~
}
\item{table.row.space}{
%% ~~Describe \code{table.row.space} here~~
}
\item{table.nrow}{
%% ~~Describe \code{table.nrow} here~~
}
\item{table.ncol}{
%% ~~Describe \code{table.ncol} here~~
}
\item{root.code}{
%% ~~Describe \code{root.code} here~~
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
function (g, node.level = "term_distance", node.level.value = 2,
node.label.size = 2, node.label.color = "darkblue", node.label.alpha =
0.8,
node.label.padding = 0, node.label.arrow = 0.01, node.label.force = 0,
node.shape = 19, node.xcoord = NULL, node.ycoord = NULL,
node.color = NULL, node.color.title = NULL, colormap = "grey-grey",
ncolors = 64, zlim = NULL, node.size.range = 4, title = "",
edge.size = 0.5, edge.color = "black", edge.color.alpha = 0.4,
edge.curve = 0.1, edge.arrow = 2, edge.arrow.gap = 0.02,
node.table = "term_name", node.table.wrap = 50, table.base.size = 7,
table.row.space = 2, table.nrow = 55, table.ncol = NULL,
root.code = "RT")
{
if (any(class(g) \%in\% c("igraph"))) {
ig <- g
}
else {
stop("The function must apply to a 'igraph' object.\n")
}
if (is.null(node.level)) {
stop("The node level must be provided!\n")
}
else {
node.level <- igraph::vertex_attr(ig, node.level)
if (is.null(node.level)) {
stop("The node level provided does not exist!\n")
}
}
ind <- which(node.level == node.level.value)
if (length(ind) == 0) {
stop("The node level value provided is wrong!\n")
}
vec.name <- V(ig)$name[ind]
V(ig)$node.code <- root.code
for (k in 1:length(vec.name)) {
x <- vec.name[k]
neighs.out <- igraph::neighborhood(ig, order = vcount(ig),
nodes = x, mode = "out")
neighbors <- names(unlist(neighs.out))
tmp <- igraph::distances(ig, v = x, to = V(ig), mode = "out")[1,
]
ind <- match(names(tmp), V(ig)$name[V(ig)$code != root.code])
tmp <- tmp[is.na(ind)]
tmp <- tmp[!is.infinite(tmp)]
if (k <= 26) {
code1 <- letters[k]
}
else {
code1 <- LETTERS[k - 26]
}
V(ig)[names(tmp)[tmp == 0]]$node.code <- code1
ind <- which(tmp != 0)
ttmp <- base::make.unique(as.character(tmp[ind]), sep = "-")
tmp_ls <- lapply(strsplit(ttmp, "-"), function(x) {
x <- as.numeric(x)
a <- x[1]
if (a >= 10 & a <= 35) {
a <- letters[a - 9]
}
else if (a > 35) {
a <- LETTERS[a - 35]
}
if (length(x) == 1) {
c(a, 1)
}
else {
b <- x[2] + 1
if (b >= 10 & b <= 35) {
b <- letters[b - 9]
}
else if (b > 35 & b <= 61) {
b <- LETTERS[b - 35]
}
c(a, b)
}
})
tmp_mat <- do.call(rbind, tmp_ls)
V(ig)[names(tmp)[ind]]$node.code <- paste0(code1, tmp_mat[,
1], tmp_mat[, 2])
}
gp_code <- xGGnetwork(g = ig, node.label = "node.code",
label.wrap.width = 30,
node.label.size = node.label.size, node.label.color = node.label.color,
node.label.alpha = node.label.alpha, node.label.padding =
node.label.padding,
node.label.arrow = node.label.arrow, node.label.force =
node.label.force,
node.shape = node.shape, node.xcoord = node.xcoord, node.ycoord =
node.ycoord,
node.color = node.color, node.color.title = node.color.title,
colormap = colormap, ncolors = ncolors, zlim = zlim,
node.size.range = node.size.range, title = title, edge.size =
edge.size,
edge.color = edge.color, edge.color.alpha = edge.color.alpha,
edge.curve = edge.curve, edge.arrow = edge.arrow, edge.arrow.gap =
edge.arrow.gap)
if (is.null(node.table)) {
stop("The node name must be provided for the table!\n")
}
else {
node.table <- igraph::vertex_attr(ig, node.table)
if (is.null(node.table)) {
stop("The node name for the table provided does not exist!\n")
}
}
res_list <- lapply(node.table, function(x) {
if (!is.na(x)) {
x <- gsub("_", " ", x)
y <- strwrap(x, width = node.table.wrap)
if (length(y) == 2) {
paste(y, collapse = "\n")
}
else if (length(y) > 2) {
paste0(paste(y[1:2], collapse = "\n"), "...")
}
else {
y
}
}
else {
x
}
})
V(ig)$node.table <- unlist(res_list)
df_code <- data.frame(Code = V(ig)$node.code, Name = V(ig)$node.table,
stringsAsFactors = FALSE)
Code <- NULL
df_code <- df_code \%>\% dplyr::arrange(Code)
tt <- gridExtra::ttheme_default(base_size = table.base.size,
padding = unit(c(1, 1), "mm"), core = list(fg_params = list(lineheight
= 0.8,
fontfamily = "sans"), bg_params = list(fill = c("snow1",
"snow2"))), colhead = list(bg_params = list(fill = "snow3"),
fg_params = list(cex = 1)), rowhead = list(fg_params = list(hjust = 0,
x = 0, fontface = "bold.italic"), bg_params = list(fill = c("snow2",
"snow1"))), )
if (0) {
ind <- base::nchar(df_code$Code) == 1 | df_code$Code ==
root.code
df_code_sub <- df_code[ind, ]
df_code_other <- df_code[!ind, ]
ls_gt <- lapply(df_code_sub$Code, function(x) {
ind <- grep(x, df_code_other$Code)
if (length(ind) > 0) {
gridExtra::tableGrob(df_code_other[ind, ], theme = tt,
rows = NULL)
}
else {
NULL
}
})
ls_gt <- base::Filter(base::Negate(is.null), ls_gt)
gridExtra::grid.arrange(grobs = ls_gt, ncol = 3, as.table = TRUE)
gx <- grid::textGrob("aa")
}
else {
ind_root <- df_code$Code == root.code
ind_sub <- base::nchar(df_code$Code) == 1
df_code_sub <- df_code[ind_sub, ]
df_code_other <- df_code[!(ind_sub | ind_root), ]
ls_df <- lapply(df_code_sub$Code, function(x) {
ind <- base::grep(paste0("^", x), df_code_other$Code)
df_code_other[ind, ]
})
df_code_order <- rbind(df_code[ind_root, ], df_code_sub,
do.call(rbind, ls_df))
if (0) {
vec <- ggplot2::cut_interval(1:nrow(df_code_order),
length = 20)
}
else {
vec_length <- unlist(lapply(strsplit(df_code_order$Name,
"\n"), length))
vec_length[vec_length == 2] <- table.row.space
vec <- ggplot2::cut_interval(cumsum(vec_length),
length = table.nrow)
}
ls_gt <- lapply(unique(vec), function(x) {
ind <- which(vec == x)
y <- df_code_order[ind, ]
diff <- table.nrow - sum(vec_length[ind])
if (diff > 0) {
z <- matrix(rep(c("", ""), diff), ncol = 2)
colnames(z) <- c("Code", "Name")
y <- rbind(y, z)
}
if (0) {
rownames(y) <- y$Code
y <- subset(y, select = "Name")
gridExtra::tableGrob(y, theme = tt)
}
else {
gridExtra::tableGrob(y, theme = tt, rows = NULL)
}
})
ls_gt <- base::Filter(base::Negate(is.null), ls_gt)
if (is.null(table.ncol)) {
table.ncol <- 3
}
if (length(unique(vec)) < table.ncol) {
table.ncol <- length(unique(vec))
}
top <- ""
if (table.ncol > 3) {
top <- quote(paste("page", g, "of", pages))
}
gp_table <- gridExtra::marrangeGrob(ls_gt, ncol = table.ncol,
nrow = 1, as.table = FALSE, top = top)
}
ls_gp <- list(code = gp_code, table = gp_table, ig = ig)
invisible(ls_gp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
