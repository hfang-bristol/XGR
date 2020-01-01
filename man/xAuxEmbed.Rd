\name{xAuxEmbed}
\alias{xAuxEmbed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xAuxEmbed(file, download.attribute = basename(file), link.text =
paste0("Download ", download.attribute))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{file}{
%% ~~Describe \code{file} here~~
}
\item{download.attribute}{
%% ~~Describe \code{download.attribute} here~~
}
\item{link.text}{
%% ~~Describe \code{link.text} here~~
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
function (file, download.attribute = basename(file), link.text =
paste0("Download ",
download.attribute))
{
ext <- tools::file_ext(file)
if (ext \%in\% c("txt")) {
type <- "text/plain"
}
else if (ext \%in\% c("csv")) {
type <- "text/csv"
}
else if (ext \%in\% c("html")) {
type <- "text/html"
}
else if (ext \%in\% c("pdf")) {
type <- "application/pdf"
}
else if (ext \%in\% c("Rmd", "rmd")) {
type <- "text/x-markdown"
}
else {
type <- "text/plain"
}
encoded <- sprintf("data:\%s;base64,\%s", type,
base64enc::base64encode(file))
res <- sprintf("<a href='\%s' download='\%s'>\%s</a>", encoded,
download.attribute, link.text)
invisible(res)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
