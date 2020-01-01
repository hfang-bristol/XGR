\name{xPieplot}
\alias{xPieplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xPieplot(df, columns, colormap = "ggplot2", pie.radius = NULL,
pie.color = "transparent", pie.color.alpha = 1, pie.thick = 0.1,
legend.title = "", gp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{df}{
%% ~~Describe \code{df} here~~
}
\item{columns}{
%% ~~Describe \code{columns} here~~
}
\item{colormap}{
%% ~~Describe \code{colormap} here~~
}
\item{pie.radius}{
%% ~~Describe \code{pie.radius} here~~
}
\item{pie.color}{
%% ~~Describe \code{pie.color} here~~
}
\item{pie.color.alpha}{
%% ~~Describe \code{pie.color.alpha} here~~
}
\item{pie.thick}{
%% ~~Describe \code{pie.thick} here~~
}
\item{legend.title}{
%% ~~Describe \code{legend.title} here~~
}
\item{gp}{
%% ~~Describe \code{gp} here~~
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
function (df, columns, colormap = "ggplot2", pie.radius = NULL,
pie.color = "transparent", pie.color.alpha = 1, pie.thick = 0.1,
legend.title = "", gp = NULL)
{
if (class(df) != "data.frame") {
stop("The function must apply to a data frame.\n")
}
ind <- colnames(df) \%in\% c("x", "y", columns)
if (sum(ind) != (length(columns) + 2)) {
stop("The input data frame must contain columns 'x', 'y' and those in
columns provided.\n")
}
df_sub <- df[, ind]
df_sub <- df_sub[!duplicated(df_sub), ]
df_sub <- df_sub[apply(df_sub[, c(-1, -2)], 1, sum) != 0,
]
geom_pie <- function(mapping = NULL, data, columns, ...) {
if (is.null(mapping)) {
mapping <- ggplot2::aes_(x = ~x, y = ~y)
}
mapping <- utils::modifyList(mapping, ggplot2::aes_(r0 = 0,
fill = ~type, amount = ~value))
if (!("r" \%in\% names(mapping))) {
xvar <- "x"
tmp <- base::diff(range(data[, xvar]))
if (tmp == 0) {
tmp <- 1
}
size <- tmp/75
mapping <- utils::modifyList(mapping, aes_(r = size))
}
names(mapping)[match(c("x", "y"), names(mapping))] <- c("x0",
"y0")
df <- tidyr::gather(data, "type", "value", columns, factor_key = TRUE)
ggforce::geom_arc_bar(mapping, data = df, stat = "pie",
inherit.aes = FALSE, ...)
}
x <- y <- NULL
if (all(class(gp) \%in\% c("gg", "ggplot"))) {
gp <- gp
}
else {
gp <- ggplot()
}
if (is.null(pie.radius)) {
gp <- gp + geom_pie(data = df_sub, aes(x = x, y = y),
columns = columns, show.legend = TRUE, color = pie.color,
alpha = pie.color.alpha, size = pie.thick)
}
else {
df_sub$pie.radius <- rep(pie.radius, nrow(df_sub))
gp <- gp + geom_pie(data = df_sub, aes(x = x, y = y,
r = pie.radius), columns = columns, show.legend = TRUE,
color = pie.color, alpha = pie.color.alpha, size = pie.thick)
}
gp <- gp + coord_equal() + guides(fill = guide_legend(title =
legend.title,
title.position = "top", keywidth = 0.6, keyheight = 0.6,
override.aes = list(shape = 19)))
my_colors <- xColormap(colormap)(length(columns))
names(my_colors) <- columns
gp <- gp + scale_fill_manual(values = my_colors)
invisible(gp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
