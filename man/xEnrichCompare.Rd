\name{xEnrichCompare}
\alias{xEnrichCompare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichCompare(list_eTerm, displayBy = c("fc", "adjp", "fdr", "zscore",
"pvalue"), FDR.cutoff = 0.05, bar.label = TRUE, bar.label.size = 2.5,
wrap.width = NULL, sharings = NULL, font.family = "sans", facet = TRUE,
signature = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{list_eTerm}{
%% ~~Describe \code{list_eTerm} here~~
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
\item{wrap.width}{
%% ~~Describe \code{wrap.width} here~~
}
\item{sharings}{
%% ~~Describe \code{sharings} here~~
}
\item{font.family}{
%% ~~Describe \code{font.family} here~~
}
\item{facet}{
%% ~~Describe \code{facet} here~~
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
function (list_eTerm, displayBy = c("fc", "adjp", "fdr", "zscore",
"pvalue"), FDR.cutoff = 0.05, bar.label = TRUE, bar.label.size = 2.5,
wrap.width = NULL, sharings = NULL, font.family = "sans",
facet = TRUE, signature = TRUE)
{
displayBy <- match.arg(displayBy)
list_eTerm <- base::Filter(base::Negate(is.null), list_eTerm)
if (length(list_eTerm) == 0) {
return(NULL)
}
list_names <- names(list_eTerm)
if (is.null(list_names)) {
list_names <- paste("Enrichment", 1:length(list_eTerm),
sep = " ")
}
res_ls <- lapply(1:length(list_eTerm), function(i) {
df <- xEnrichViewer(list_eTerm[[i]], top_num = "all",
sortBy = "none")
if (is.null(df)) {
return(NULL)
}
else {
cbind(group = rep(list_names[i], nrow(df)), id = rownames(df),
df, stringsAsFactors = F)
}
})
df_all <- do.call(rbind, res_ls)
rownames(df_all) <- NULL
ind <- which(df_all$adjp < FDR.cutoff)
d <- df_all[ind, c("id", "name", "fc", "adjp", "zscore",
"pvalue", "group")]
d$group <- factor(d$group, levels = rev(list_names))
if (!is.null(wrap.width)) {
width <- as.integer(wrap.width)
res_list <- lapply(d$name, function(x) {
x <- gsub("_", " ", x)
y <- strwrap(x, width = width)
if (length(y) > 1) {
paste0(y[1], "...")
}
else {
y
}
})
d$name <- unlist(res_list)
}
id_ls <- split(x = d$group, f = d$id)
ind <- match(d$id, names(id_ls))
id_full_ls <- id_ls[ind]
nSig <- unlist(lapply(id_full_ls, length))
d$nSig <- nSig
code <- lapply(id_full_ls, function(x) {
res <- rep(0, length(levels(x)))
ind <- match(x, levels(x))
res[ind] <- 1
paste(res, collapse = "-")
})
d$code <- unlist(code)
if (!is.null(sharings)) {
sharings <- as.numeric(sharings)
ind <- match(sharings, unique(d$nSig))
found <- sharings[!is.na(ind)]
if (length(found) > 0) {
flag <- match(d$nSig, found)
d <- d[!is.na(flag), ]
nSig <- nSig[!is.na(flag)]
}
}
name <- fc <- group <- adjp <- zscore <- pvalue <- label <- direction
<- height <- hjust <- NULL
d <- d \%>\% dplyr::mutate(direction = ifelse(zscore > 0, 1,
-1))
if (displayBy == "fc") {
d <- d \%>\% dplyr::arrange(nSig, group, direction, fc,
desc(adjp)) \%>\% dplyr::mutate(height = log2(fc)) \%>\%
dplyr::mutate(hjust = ifelse(height >= 0, 1, 0))
d$name <- factor(d$name, levels = unique(d$name))
ind <- match(unique(d$id), names(nSig))
xintercept <- which(!duplicated(nSig[ind]))[-1]
p <- ggplot(d, aes(x = name, y = height, fill = group))
p <- p + ylab(expression(paste("Enrichment changes: ",
log[2]("FC"))))
}
else if (displayBy == "adjp" | displayBy == "fdr") {
d <- d \%>\% dplyr::arrange(nSig, group, direction, desc(adjp),
zscore) \%>\% dplyr::mutate(height = -1 * log10(adjp)) \%>\%
dplyr::mutate(hjust = 1)
d$name <- factor(d$name, levels = unique(d$name))
ind <- match(unique(d$id), names(nSig))
xintercept <- which(!duplicated(nSig[ind]))[-1]
p <- ggplot(d, aes(x = name, y = height, fill = group))
p <- p + ylab(expression(paste("Enrichment significance: ",
-log[10]("FDR"))))
}
else if (displayBy == "zscore") {
d <- d \%>\% dplyr::arrange(nSig, group, direction, zscore,
desc(adjp)) \%>\% dplyr::mutate(height = zscore) \%>\%
dplyr::mutate(hjust = ifelse(height >= 0, 1, 0))
d$name <- factor(d$name, levels = unique(d$name))
ind <- match(unique(d$id), names(nSig))
xintercept <- which(!duplicated(nSig[ind]))[-1]
p <- ggplot(d, aes(x = name, y = height, fill = group))
p <- p + ylab("Enrichment z-scores")
}
else if (displayBy == "pvalue") {
d <- d \%>\% dplyr::arrange(nSig, group, direction, desc(pvalue),
zscore) \%>\% dplyr::mutate(height = -1 * log10(pvalue)) \%>\%
dplyr::mutate(hjust = 1)
d$name <- factor(d$name, levels = unique(d$name))
ind <- match(unique(d$id), names(nSig))
xintercept <- which(!duplicated(nSig[ind]))[-1]
p <- ggplot(d, aes(x = name, y = height, fill = group))
p <- p + ylab(expression(paste("Enrichment significance: ",
-log[10]("p-value"))))
}
bp <- p + geom_bar(stat = "identity") + theme_bw() +
theme(legend.position = "none",
legend.title = element_blank(), axis.title.y = element_blank(),
axis.text.y = element_text(size = 10, color = "black"),
axis.title.x = element_text(size = 14, color = "black")) +
geom_vline(xintercept = xintercept - 0.5, color = "black",
linetype = "dotdash") + coord_flip()
bp <- bp + theme(panel.grid.major = element_blank(), panel.grid.minor =
element_blank())
bp <- bp + theme(strip.background = element_rect(fill = "transparent",
color = "transparent"), strip.text = element_text(size = 12,
face = "italic"))
if (bar.label) {
to_scientific_notation <- function(x) {
res <- format(x, scientific = T)
res <- sub("\\+0?", "", res)
sub("-0?", "-", res)
}
label <- to_scientific_notation(d$adjp)
d$label <- paste("FDR", as.character(label), sep = "=")
df_sub <- d \%>\% dplyr::filter(hjust == 1)
bp <- bp + geom_text(data = df_sub, aes(label = label),
hjust = 1, size = bar.label.size, family = font.family)
df_sub <- d \%>\% dplyr::filter(hjust == 0)
bp <- bp + geom_text(data = df_sub, aes(label = label),
hjust = 0, size = bar.label.size, family = font.family)
}
title <- paste0("Enrichments (FDR < ", FDR.cutoff, ")")
bp <- bp + labs(title = title) + theme(plot.title = element_text(hjust
= 0.5))
if (signature) {
caption <- paste("Created by xEnrichCompare from XGR version",
utils::packageVersion("XGR"))
bp <- bp + labs(caption = caption) + theme(plot.caption =
element_text(hjust = 1,
face = "bold.italic", size = 8, colour = "#002147"))
}
bp <- bp + theme(text = element_text(family = font.family))
bp <- bp + theme(axis.line.x = element_line(arrow = arrow(angle = 30,
length = unit(0.25, "cm"), type = "open")))
if (facet) {
bp <- bp + facet_grid(~group)
}
else {
bp <- bp + theme(legend.position = "bottom", legend.title =
element_text(size = 10,
color = "black", face = "bold")) + guides(fill = guide_legend(title =
"",
title.position = "left", reverse = TRUE, nrow = 1,
byrow = TRUE))
}
flag <- sapply(list_eTerm, function(x) !is.null(x$g))
if (all(flag)) {
ls_edges <- lapply(list_eTerm, function(x) {
df_edge <- igraph::get.data.frame(x$g, what = "edges")
})
relations <- do.call(rbind, ls_edges)
relations <- relations[!duplicated(relations), ]
ls_nodes <- lapply(list_eTerm, function(x) {
df_nodes <- igraph::get.data.frame(x$g, what = "vertices")
})
nodes <- do.call(rbind, ls_nodes)[, 1:4]
nodes <- nodes[!duplicated(nodes), ]
ig <- igraph::graph.data.frame(d = relations, directed = T,
vertices = nodes)
bp$g <- ig
}
invisible(bp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
