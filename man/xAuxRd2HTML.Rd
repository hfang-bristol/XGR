\name{xAuxRd2HTML}
\alias{xAuxRd2HTML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xAuxRd2HTML(path.from = "./XGR/man", path.to = "./XGR/vignettes")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{path.from}{
%% ~~Describe \code{path.from} here~~
}
\item{path.to}{
%% ~~Describe \code{path.to} here~~
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
function (path.from = "./XGR/man", path.to = "./XGR/vignettes")
{
for (nm in list.files(path.from, pattern = "\\.Rd$")) {
nm.out <- gsub(pattern = "\\.Rd$", "\\.html", nm, ignore.case = T,
perl = T)
input_file <- file.path(path.from, nm)
output_file <- file.path(path.to, nm.out)
message(sprintf("\%s => \%s", input_file, output_file),
appendLF = T)
tools::Rd2HTML(Rd = input_file, out = output_file, package = "",
stages = c("install", "render"))
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
