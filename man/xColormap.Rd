\name{xColormap}
\alias{xColormap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xColormap(colormap = c("bwr", "jet", "gbr", "wyr", "br", "yr",
"rainbow", "wb", "heat", "terrain", "topo", "cm", "ggplot2"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{colormap}{
%% ~~Describe \code{colormap} here~~
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
function (colormap = c("bwr", "jet", "gbr", "wyr", "br", "yr",
"rainbow", "wb", "heat", "terrain", "topo", "cm", "ggplot2"))
{
if (length(colormap) > 1) {
colormap <- colormap[1]
}
if (colormap == "ggplot2") {
my_hue_pal <- function(h = c(0, 360) + 15, c = 100, l = 65,
h.start = 0, direction = 1) {
function(n) {
if ((diff(h)\%\%360) < 1) {
h[2] <- h[2] - 360/n
}
rotate <- function(x) (x + h.start)\%\%360 * direction
hues <- rotate(seq(h[1], h[2], length.out = n))
grDevices::hcl(hues, c, l)
}
}
palette.name <- my_hue_pal(h = c(0, 360) + 15, c = 100,
l = 65, h.start = 0, direction = 1)
}
else if (colormap == "heat") {
palette.name <- grDevices::heat.colors
}
else if (colormap == "terrain") {
palette.name <- grDevices::terrain.colors
}
else if (colormap == "topo") {
palette.name <- grDevices::topo.colors
}
else if (colormap == "cm") {
palette.name <- grDevices::cm.colors
}
else {
palette.name <- supraHex::visColormap(colormap = colormap)
}
invisible(palette.name)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
