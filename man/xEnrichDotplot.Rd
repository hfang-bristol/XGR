\name{xEnrichDotplot}
\alias{xEnrichDotplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnrichDotplot(eTerm, FDR.cutoff = 0.05, colors = c("pink", "red"),
y.scale = c("normal", "log"), slim = NULL, size.range = c(0.5, 3.5),
size.title = "Num of overlaps", label.top = "auto", label.direction.y =
c("left", "right", "none"), label.size = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{eTerm}{
%% ~~Describe \code{eTerm} here~~
}
\item{FDR.cutoff}{
%% ~~Describe \code{FDR.cutoff} here~~
}
\item{colors}{
%% ~~Describe \code{colors} here~~
}
\item{y.scale}{
%% ~~Describe \code{y.scale} here~~
}
\item{slim}{
%% ~~Describe \code{slim} here~~
}
\item{size.range}{
%% ~~Describe \code{size.range} here~~
}
\item{size.title}{
%% ~~Describe \code{size.title} here~~
}
\item{label.top}{
%% ~~Describe \code{label.top} here~~
}
\item{label.direction.y}{
%% ~~Describe \code{label.direction.y} here~~
}
\item{label.size}{
%% ~~Describe \code{label.size} here~~
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
function (eTerm, FDR.cutoff = 0.05, colors = c("pink", "red"),
y.scale = c("normal", "log"), slim = NULL, size.range = c(0.5,
3.5), size.title = "Num of overlaps", label.top = "auto",
label.direction.y = c("left", "right", "none"), label.size = 2,
...)
{
y.scale <- match.arg(y.scale)
label.direction.y <- match.arg(label.direction.y)
if (is.null(eTerm)) {
warnings("There is no enrichment in the 'eTerm' object.\n")
return(NULL)
}
if (class(eTerm) == "eTerm") {
df_enrichment_group <- xEnrichViewer(eTerm, top_num = "all",
sortBy = "adjp")
df_enrichment_group$group <- "group"
}
else if (class(eTerm) == "ls_eTerm" | class(eTerm) == "data.frame") {
if (class(eTerm) == "ls_eTerm") {
df_enrichment_group <- eTerm$df
}
else if (class(eTerm) == "data.frame") {
if (all(c("group", "name", "adjp", "nOverlap", "zscore") \%in\%
colnames(eTerm))) {
df_enrichment_group <- eTerm[, c("group", "name",
"adjp", "nOverlap", "zscore")]
}
else if (all(c("name", "adjp", "nOverlap", "zscore") \%in\%
colnames(eTerm))) {
df_enrichment_group <- eTerm[, c("name", "adjp",
"nOverlap", "zscore")]
df_enrichment_group$group <- "group"
}
}
}
if (class(df_enrichment_group$group) == "factor") {
if (length(unique(df_enrichment_group$group)) !=
length(levels(df_enrichment_group$group))) {
df_enrichment_group$group <- factor(df_enrichment_group$group,
levels = sort(unique(df_enrichment_group$group)))
}
}
ngroup <- length(unique(df_enrichment_group$group))
adjp <- zscore <- nOverlap <- flag <- name <- group <- rank <- NULL
df_enrichment_group <- df_enrichment_group \%>\% dplyr::mutate(flag =
ifelse(adjp >=
FDR.cutoff, "N", "Y"))
names(colors) <- c("N", "Y")
gp <- ggplot(df_enrichment_group, aes(x = zscore, y = -log10(adjp),
size = nOverlap))
gp <- gp + geom_point(aes(color = flag, size = nOverlap),
alpha = 0.6)
gp <- gp + scale_colour_manual(values = colors) + guides(color =
"none")
gp <- gp + xlab("Z-score") + ylab(expression(-log[10]("FDR")))
gp <- gp + theme_bw() + theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank())
if (is.null(slim)) {
slim <- range(df_enrichment_group$nOverlap)
}
if (any(is.na(slim))) {
gp <- gp + theme(legend.position = "none")
}
else {
gp <- gp + theme(legend.position = "bottom", legend.key =
element_rect(colour = "transparent")) +
scale_size_continuous(limits = slim, range = size.range,
guide = guide_legend(size.title, title.position = "left",
nrow = 1))
}
if (label.top == "auto") {
df <- subset(df_enrichment_group, adjp < FDR.cutoff)
}
else {
label.top <- as.integer(label.top)
df <- as.data.frame(df_enrichment_group \%>\% dplyr::group_by(group)
\%>\%
dplyr::group_by(rank = rank(adjp), add = TRUE) \%>\%
dplyr::filter(rank <= label.top & adjp < FDR.cutoff))
}
if (ngroup == 1 & label.direction.y != "none") {
offset <- (range(df_enrichment_group$zscore)[2] -
range(df_enrichment_group$zscore)[1]) *
0.1
if (label.direction.y == "right") {
df$nudge_x <- max(df_enrichment_group$zscore) - df$zscore +
offset
gp <- gp + ggrepel::geom_text_repel(data = df, aes(x = zscore,
y = -log10(adjp), label = name), size = label.size,
show.legend = F, segment.alpha = 0.5, segment.color = "grey50",
segment.size = 0.2, arrow = arrow(length = unit(0.01,
"npc")), direction = "y", hjust = 1, nudge_x = df$nudge_x,
...)
gp <- gp + scale_x_continuous(position = "bottom",
limits = c(min(df_enrichment_group$zscore),
max(df_enrichment_group$zscore) +
offset))
}
else if (label.direction.y == "left") {
df$nudge_x <- -1 * (df$zscore - min(df_enrichment_group$zscore)) -
offset
gp <- gp + ggrepel::geom_text_repel(data = df, aes(x = zscore,
y = -log10(adjp), label = name), size = label.size,
show.legend = F, segment.alpha = 0.5, segment.color = "grey50",
segment.size = 0.2, arrow = arrow(length = unit(0.01,
"npc")), direction = "y", hjust = 0, nudge_x = df$nudge_x,
...)
gp <- gp + scale_x_continuous(position = "bottom",
limits = c(min(df_enrichment_group$zscore) -
offset, max(df_enrichment_group$zscore)))
}
}
else {
gp <- gp + ggrepel::geom_text_repel(data = df, aes(x = zscore,
y = -log10(adjp), label = name), size = label.size,
show.legend = F, segment.alpha = 0.5, segment.color = "grey50",
segment.size = 0.2, arrow = arrow(length = unit(0.01,
"npc")), ...)
}
if (y.scale == "log") {
gp <- gp + scale_y_continuous(trans = "log1p")
}
if (ngroup != 1) {
scales <- "free_y"
space <- "free_y"
gp <- gp + facet_wrap(~group)
}
return(gp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
