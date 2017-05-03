\name{xVisInterp}
\alias{xVisInterp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xVisInterp(ls_xyz, interpolation = c("spline", "linear"), nx = 100, ny
= 100, zlim = NULL, nD = c("auto", "2D", "3D"), colkey = TRUE, contour
= FALSE, image = FALSE, clab = c("Value", ""), nlevels = 20, colormap =
"terrain", theta.3D = 40, phi.3D = 20, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ls_xyz}{
%% ~~Describe \code{ls_xyz} here~~
}
\item{interpolation}{
%% ~~Describe \code{interpolation} here~~
}
\item{nx}{
%% ~~Describe \code{nx} here~~
}
\item{ny}{
%% ~~Describe \code{ny} here~~
}
\item{zlim}{
%% ~~Describe \code{zlim} here~~
}
\item{nD}{
%% ~~Describe \code{nD} here~~
}
\item{colkey}{
%% ~~Describe \code{colkey} here~~
}
\item{contour}{
%% ~~Describe \code{contour} here~~
}
\item{image}{
%% ~~Describe \code{image} here~~
}
\item{clab}{
%% ~~Describe \code{clab} here~~
}
\item{nlevels}{
%% ~~Describe \code{nlevels} here~~
}
\item{colormap}{
%% ~~Describe \code{colormap} here~~
}
\item{theta.3D}{
%% ~~Describe \code{theta.3D} here~~
}
\item{phi.3D}{
%% ~~Describe \code{phi.3D} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
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
function (ls_xyz, interpolation = c("spline", "linear"), nx = 100,
ny = 100, zlim = NULL, nD = c("auto", "2D", "3D"), colkey = TRUE,
contour = FALSE, image = FALSE, clab = c("Value", ""), nlevels = 20,
colormap = "terrain", theta.3D = 40, phi.3D = 20, verbose = TRUE)
{
interpolation <- match.arg(interpolation)
nD <- match.arg(nD)
if (class(ls_xyz) == "list") {
ind <- names(ls_xyz) \%in\% c("x", "y", "z")
ls_xyz <- ls_xyz[ind]
if (length(ls_xyz) != 3) {
return(NULL)
}
}
else if (class(ls_xyz) == "data.frame" | class(ls_xyz) ==
"matrix") {
if (ncol(ls_xyz) == 3) {
ls_xyz <- ls_xyz[, 1:3]
colnames(ls_xyz) <- c("x", "y", "z")
}
else if (ncol(ls_xyz) == 4) {
ls_xyz <- ls_xyz[, 1:4]
colnames(ls_xyz) <- c("x", "y", "z", "label")
}
ls_xyz <- as.list(ls_xyz)
}
else {
stop("The function must apply to a 'list' object, or a
'data.frame'/'matrix' object.\n")
}
if (interpolation == "linear") {
linear <- TRUE
}
else {
linear <- FALSE
}
x <- y <- z <- label <- NULL
ls_xyz.smooth <- with(ls_xyz, akima::interp(x, y, z, nx = nx,
ny = ny, linear = linear))
col <- xColormap(colormap)(nlevels)
if (class(colkey) != "list") {
if (colkey) {
colkey <- list(side = 4, length = 0.15, width = 0.5,
shift = 0.35, dist = -0.15, cex.axis = 0.6, cex.clab = 0.8,
side.clab = 3)
}
}
if (nD == "auto") {
si.zmin <- min(ls_xyz.smooth$z, na.rm = TRUE)
si.zmax <- max(ls_xyz.smooth$z, na.rm = TRUE)
breaks <- pretty(c(si.zmin, si.zmax), nlevels + 1)
colors <- xColormap(colormap)(length(breaks) - 1)
if (image) {
graphics::image(ls_xyz.smooth, axes = FALSE, breaks = breaks,
col = colors, xlim = grDevices::extendrange(ls_xyz.smooth$x,
f = 0.1), ylim = grDevices::extendrange(ls_xyz.smooth$y,
f = 0.1))
}
if (contour) {
if (image) {
graphics::contour(ls_xyz.smooth, add = TRUE,
labcex = 1, levels = breaks, col = "thistle")
}
else {
graphics::contour(ls_xyz.smooth, labcex = 1,
levels = breaks, col = colors)
}
graphics::points(ls_xyz, pch = 20, cex = 1, col = "black")
}
else {
if (!image) {
plot(y ~ x, data = ls_xyz, pch = 20, col = "blue",
axes = FALSE, ann = FALSE)
}
else {
graphics::points(ls_xyz, pch = 20, cex = 1)
}
if (is.null(ls_xyz$label)) {
with(ls_xyz, graphics::text(x, y, formatC(z,
dig = 2), adj = -0.2, cex = 0.6))
}
else {
with(ls_xyz, graphics::text(x, y, label, adj = -0.2,
cex = 0.6))
}
}
}
else if (nD == "2D") {
if (class(contour) != "list") {
if (contour) {
contour <- list(col = "black", labcex = 1, lwd = 2,
alpha = 0.5, addbox = FALSE, nlevels = nlevels)
}
}
graphics::par(mar = c(0, 0, 0, 0))
plot3D::image2D(z = ls_xyz.smooth$z, axes = FALSE, cex.axis = 0.8,
cex.lab = 1.2, col = col, colkey = colkey, clab = clab,
plot = TRUE, lighting = FALSE, lphi = 90, contour = contour)
}
else {
if (is.null(zlim)) {
zlim <- c(floor(min(ls_xyz.smooth$z, na.rm = TRUE) *
100)/100, ceiling(max(ls_xyz.smooth$z, na.rm = TRUE) *
100)/100)
if (verbose) {
now <- Sys.time()
message(sprintf("The range of interpolated values: [\%.3f, \%.3f]",
zlim[1], zlim[2]), appendLF = TRUE)
}
if (class(contour) != "list" & class(image) != "list") {
if (contour | image) {
zlim <- c(1.5 * zlim[1] - 0.5 * zlim[2], zlim[2])
}
}
}
if (class(contour) == "list" | class(image) == "list") {
graphics::par(mar = c(0, 0, 0, 0))
plot3D::persp3D(z = ls_xyz.smooth$z, axes = FALSE,
box = FALSE, zlim = zlim, cex.axis = 0.8, cex.lab = 1.2,
ticktype = c("simple", "detailed")[1], col = col,
colkey = colkey, clab = clab, bty = "b", facets = TRUE,
curtain = FALSE, plot = TRUE, lighting = TRUE,
lphi = 90, theta = theta.3D, phi = phi.3D, d = 1,
image = image, contour = contour)
}
else {
if (contour | image) {
plot <- FALSE
}
else {
plot <- TRUE
}
graphics::par(mar = c(0, 0, 0, 0))
plot3D::persp3D(z = ls_xyz.smooth$z, axes = FALSE,
box = FALSE, zlim = zlim, cex.axis = 0.8, cex.lab = 1.2,
ticktype = c("simple", "detailed")[1], col = col,
colkey = colkey, clab = clab, bty = "b", facets = TRUE,
curtain = FALSE, plot = plot, lighting = TRUE,
lphi = 90, theta = theta.3D, phi = phi.3D, d = 1)
if (image) {
if (contour) {
plot2 <- FALSE
}
else {
plot2 <- TRUE
}
plot3D::image3D(z = zlim[1], colvar = ls_xyz.smooth$z,
col = col, box = FALSE, colkey = FALSE, add = TRUE,
plot = plot2)
}
if (contour) {
if (image) {
col <- c("thistle", "grey")[2]
}
plot3D::contour3D(z = zlim[1], colvar = ls_xyz.smooth$z,
labcex = 1, col = col, nlevels = nlevels, lwd = 0.8,
box = FALSE, addbox = FALSE, colkey = FALSE,
add = TRUE, plot = TRUE)
}
}
}
invisible(ls_xyz.smooth)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
