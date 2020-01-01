\name{xAuxFunArgs}
\alias{xAuxFunArgs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xAuxFunArgs(fun, action = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{fun}{
%% ~~Describe \code{fun} here~~
}
\item{action}{
%% ~~Describe \code{action} here~~
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
function (fun, action = NULL, verbose = TRUE)
{
args_list <- base::formals(fun)
args_names <- names(args_list)
if (is.null(action)) {
ls_vec <- lapply(1:length(args_list), function(i) {
lft <- args_names[[i]]
rgt <- paste(base::deparse(args_list[[i]]), collapse = "")
if (rgt != "") {
paste0(lft, "=", rgt, "[1];")
}
else {
NULL
}
})
res <- noquote(paste(unlist(ls_vec), collapse = ""))
return(res)
}
else {
for (i in 1:length(args_list)) {
lft <- args_names[[i]]
rgt <- paste(base::deparse(args_list[[i]]), collapse = "")
if (rgt != "") {
tmp <- paste(lft, "<-", rgt, sep = " ")
if (action == T) {
base::eval(base::parse(text = tmp), envir = parent.frame())
}
else if (action == F) {
base::eval(base::parse(text = tmp))
}
}
}
}
if (verbose) {
if (action == T) {
message(sprintf("For the function '\%s', \%d arguments have been
assigned with default values in the current environment.",
fun, length(args_names)), appendLF = T)
}
else {
message(sprintf("For the function '\%s', there are \%d arguments.",
fun, length(args_names)), appendLF = T)
}
}
invisible(args_list)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
