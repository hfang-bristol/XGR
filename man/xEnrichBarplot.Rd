\name{xEnrichBarplot}
\alias{xEnrichBarplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichBarplot(eTerm, top_num = 10, displayBy = c("fc", "adjp", "fdr",
"zscore", "pvalue"), FDR.cutoff = 0.05, bar.label = TRUE,
bar.label.size = 3, bar.color = "lightyellow-orange", bar.width = 0.8,
wrap.width = NULL, font.family = "sans", signature = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{eTerm}{
%% ~~Describe \code{eTerm} here~~
}
\item{top_num}{
%% ~~Describe \code{top_num} here~~
}
\item{displayBy}{
%% ~~Describe \code{displayBy} here~~
}
\item{FDR.cutoff}{
%% ~~Describe \code{FDR.cutoff} here~~
}
\item{bar.label}{
%% ~~Describe \code{bar.label} here~~
}
\item{bar.label.size}{
%% ~~Describe \code{bar.label.size} here~~
}
\item{bar.color}{
%% ~~Describe \code{bar.color} here~~
}
\item{bar.width}{
%% ~~Describe \code{bar.width} here~~
}
\item{wrap.width}{
%% ~~Describe \code{wrap.width} here~~
}
\item{font.family}{
%% ~~Describe \code{font.family} here~~
}
\item{signature}{
%% ~~Describe \code{signature} here~~
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
function (eTerm, top_num = 10, displayBy = c("fc", "adjp", "fdr",
"zscore", "pvalue"), FDR.cutoff = 0.05, bar.label = TRUE,
bar.label.size = 3, bar.color = "lightyellow-orange", bar.width = 0.8,
wrap.width = NULL, font.family = "sans", signature = TRUE)
{
displayBy <- match.arg(displayBy)
if (is.null(eTerm)) {
warnings("There is no enrichment in the 'eTerm' object.\n")
return(NULL)
}
df <- xEnrichViewer(eTerm, top_num = "all")
if (top_num == "auto") {
top_num <- sum(df$adjp < FDR.cutoff)
if (top_num <= 1) {
top_num <- 10
}
}
df <- xEnrichViewer(eTerm, top_num = top_num, sortBy = "adjp")
if (!is.null(wrap.width)) {
width <- as.integer(wrap.width)
res_list <- lapply(df$name, function(x) {
x <- gsub("_", " ", x)
y <- strwrap(x, width = width)
if (length(y) > 1) {
paste0(y[1], "...")
}
else {
y
}
})
df$name <- unlist(res_list)
}
name <- height <- direction <- hjust <- NULL
adjp <- zscore <- pvalue <- fc <- NULL
df <- df \%>\% dplyr::mutate(direction = ifelse(zscore > 0,
1, -1))
if (displayBy == "adjp" | displayBy == "fdr") {
df <- df \%>\% dplyr::arrange(direction, desc(adjp), zscore) \%>\%
dplyr::mutate(height = -1 * log10(adjp)) \%>\% dplyr::mutate(hjust = 1)
df$name <- factor(df$name, levels = df$name)
if (length(df$height[!is.infinite(df$height)]) == 0) {
df$height <- 10
}
else {
df$height[is.infinite(df$height)] <-
max(df$height[!is.infinite(df$height)])
}
p <- ggplot(df, aes(x = name, y = height))
p <- p + ylab(expression(paste("Enrichment significance: ",
-log[10]("FDR"))))
}
else if (displayBy == "fc") {
df <- df \%>\% dplyr::arrange(direction, fc, desc(adjp)) \%>\%
dplyr::mutate(height = log2(fc)) \%>\% dplyr::mutate(hjust =
ifelse(height >=
0, 1, 0))
df$name <- factor(df$name, levels = df$name)
p <- ggplot(df, aes(x = name, y = height))
p <- p + ylab(expression(paste("Enrichment changes: ",
log[2]("FC"))))
}
else if (displayBy == "pvalue") {
df <- df \%>\% dplyr::arrange(direction, desc(pvalue),
zscore) \%>\% dplyr::mutate(height = -1 * log10(pvalue)) \%>\%
dplyr::mutate(hjust = 1)
df$name <- factor(df$name, levels = df$name)
if (length(df$height[!is.infinite(df$height)]) == 0) {
df$height <- 10
}
else {
df$height[is.infinite(df$height)] <-
max(df$height[!is.infinite(df$height)])
}
p <- ggplot(df, aes(x = name, y = height))
p <- p + ylab(expression(paste("Enrichment significance: ",
-log[10]("p-value"))))
}
else if (displayBy == "zscore") {
df <- df \%>\% dplyr::arrange(direction, zscore, desc(adjp)) \%>\%
dplyr::mutate(height = zscore) \%>\% dplyr::mutate(hjust =
ifelse(height >=
0, 1, 0))
df$name <- factor(df$name, levels = df$name)
p <- ggplot(df, aes(x = name, y = height))
p <- p + ylab("Enrichment z-scores")
}
if (is.null(bar.color)) {
bp <- p + geom_col(color = "grey80", fill = "transparent",
width = bar.width)
}
else {
bar.color <- unlist(strsplit(bar.color, "-"))
if (length(bar.color) == 2) {
bp <- p + geom_col(aes(fill = height), width = bar.width) +
scale_fill_gradient(low = bar.color[1], high = bar.color[2])
}
else {
bp <- p + geom_col(color = "grey80", fill = "transparent",
width = bar.width)
}
}
bp <- bp + theme_bw() + theme(legend.position = "none", axis.title.y =
element_blank(),
axis.text.y = element_text(size = 10, color = "black"),
axis.title.x = element_text(size = 12, color = "black")) +
coord_flip()
bp <- bp + theme(panel.grid.major = element_blank(), panel.grid.minor =
element_blank())
if (bar.label) {
to_scientific_notation <- function(x) {
res <- format(x, scientific = T)
res <- sub("\\+0?", "", res)
sub("-0?", "-", res)
}
label <- to_scientific_notation(df$adjp)
df$label <- paste("FDR", as.character(label), sep = "=")
df_sub <- df \%>\% dplyr::filter(hjust == 1)
bp <- bp + geom_text(data = df_sub, aes(label = label),
hjust = 1, size = bar.label.size, family = font.family)
df_sub <- df \%>\% dplyr::filter(hjust == 0)
bp <- bp + geom_text(data = df_sub, aes(label = label),
hjust = 0, size = bar.label.size, family = font.family)
}
if (signature) {
caption <- paste("Created by xEnrichBarplot from XGR version",
utils::packageVersion("XGR"))
bp <- bp + labs(caption = caption) + theme(plot.caption =
element_text(hjust = 1,
face = "bold.italic", size = 8, colour = "#002147"))
}
bp <- bp + theme(text = element_text(family = font.family))
bp <- bp + theme(axis.line.x = element_line(arrow = arrow(angle = 30,
length = unit(0.25, "cm"), type = "open")))
bp <- bp + scale_y_continuous(position = "top")
invisible(bp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
