\name{xVisKernels}
\alias{xVisKernels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xVisKernels(exponent = 2, newpage = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{exponent}{
%% ~~Describe \code{exponent} here~~
}
\item{newpage}{
%% ~~Describe \code{newpage} here~~
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
function (exponent = 2, newpage = T)
{
if (newpage) {
dev.new(width = 7, height = 7)
}
cl <- rainbow(3)
ph <- c(21, 22, 17)
x <- seq(0, 1, by = 0.01)
y1 <- 1 - (x)^exponent
y2 <- 1 - x
y3 <- (1 - x)^exponent
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), type = "n", xlab =
expression(paste("Relative distance of ",
d[gs], " (within window D) between gene g and SNP s",
collapse = " ")), ylab = expression(paste("Distance kernel function ",
h[gs], collapse = " ")), main = bquote("Decay exponent is" ~
lambda == .(exponent)))
lines(x, y1, type = "b", pch = ph[1], col = cl[1])
lines(x, y2, type = "b", pch = ph[2], col = cl[2])
lines(x, y3, type = "b", pch = ph[3], col = cl[3])
leg.txt <- c(expression(paste("slow decay:", h[gs], "=",
(1 - (frac(d[gs], D))^lambda) * (d[gs] <= D), collapse = " ")),
expression(paste("linear decay:", h[gs], "=", (1 - frac(d[gs],
D)) * (d[gs] <= D), collapse = " ")), expression(paste("rapid decay: ",
, h[gs], "=", ((1 - frac(d[gs], D))^lambda) * (d[gs] <=
D), collapse = " ")))
legend("topright", legend = leg.txt, pch = ph, col = cl,
border = "transparent", box.col = "transparent", cex = 0.6)
invisible()
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
