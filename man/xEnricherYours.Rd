\name{xEnricherYours}
\alias{xEnricherYours}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xEnricherYours(data.file, annotation.file, background.file = NULL,
size.range = c(10, 2000), min.overlap = 5, test = c("fisher",
"hypergeo", "binomial"), background.annotatable.only = NULL, p.tail =
c("one-tail", "two-tails"), p.adjust.method = c("BH", "BY",
"bonferroni", "holm", "hochberg", "hommel"), verbose = T, silent =
FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data.file}{
%% ~~Describe \code{data.file} here~~
}
\item{annotation.file}{
%% ~~Describe \code{annotation.file} here~~
}
\item{background.file}{
%% ~~Describe \code{background.file} here~~
}
\item{size.range}{
%% ~~Describe \code{size.range} here~~
}
\item{min.overlap}{
%% ~~Describe \code{min.overlap} here~~
}
\item{test}{
%% ~~Describe \code{test} here~~
}
\item{background.annotatable.only}{
%% ~~Describe \code{background.annotatable.only} here~~
}
\item{p.tail}{
%% ~~Describe \code{p.tail} here~~
}
\item{p.adjust.method}{
%% ~~Describe \code{p.adjust.method} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
}
\item{silent}{
%% ~~Describe \code{silent} here~~
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
function (data.file, annotation.file, background.file = NULL,
size.range = c(10, 2000), min.overlap = 5, test = c("fisher",
"hypergeo", "binomial"), background.annotatable.only = NULL,
p.tail = c("one-tail", "two-tails"), p.adjust.method = c("BH",
"BY", "bonferroni", "holm", "hochberg", "hommel"), verbose = T,
silent = FALSE)
{
startT <- Sys.time()
if (!silent) {
message(paste(c("Start at ", as.character(startT)), collapse = ""),
appendLF = TRUE)
message("", appendLF = TRUE)
}
else {
verbose <- FALSE
}
test <- match.arg(test)
p.tail <- match.arg(p.tail)
p.adjust.method <- match.arg(p.adjust.method)
p.tail <- match.arg(p.tail)
if (length(data.file) == 0) {
return(NULL)
}
if (is.matrix(data.file) | is.data.frame(data.file)) {
data <- unique(data.file[, 1])
}
else if (!is.null(data.file) & any(!is.na(data.file))) {
if (length(data.file) == 1) {
if (class(suppressWarnings(try(data <- utils::read.delim(file =
data.file,
header = F, row.names = NULL, stringsAsFactors = F),
T))) == "try-error") {
data <- data.file
}
else {
data <- unique(data[, 1])
}
}
else {
data <- data.file
}
}
else {
stop("The file 'data.file' must be provided!\n")
}
if (is.matrix(background.file) | is.data.frame(background.file)) {
background <- unique(background.file[, 1])
}
else if (!is.null(background.file)) {
if (length(background.file) == 1) {
background <- utils::read.delim(file = background.file,
header = F, row.names = NULL, stringsAsFactors = F)
background <- unique(background[, 1])
}
else {
background <- background.file
}
}
else {
background <- background.file
}
if (is.matrix(annotation.file) | is.data.frame(annotation.file)) {
input <- cbind(annotation.file[, 1], annotation.file[,
2])
}
else if (!is.null(annotation.file) & any(!is.na(annotation.file))) {
input <- utils::read.delim(file = annotation.file, header = F,
row.names = NULL, stringsAsFactors = F)
}
else {
stop("The file 'annotation.file' must be provided!\n")
}
if (0) {
if (background.annotatable.only == FALSE) {
if (!is.null(background.file)) {
tmp <- cbind(background.file, rep("BG", length(background.file)))
input <- rbind(input, tmp)
}
}
}
anno <- split(x = input[, 1], f = input[, 2])
terms <- names(anno)
nodes <- data.frame(name = terms, term_id = terms, term_name = terms,
term_distance = rep(1, length(terms)), term_namespace =
rep("Customised",
length(terms)), stringsAsFactors = F)
root <- c("Root", "Root", "Root", 0, "Customised")
nodes <- rbind(nodes, root)
relations <- data.frame(from = "Root", to = nodes$name)
g <- igraph::graph.data.frame(d = relations, directed = T,
vertices = nodes)
if (verbose) {
now <- Sys.time()
message(sprintf("\n#######################################################",
appendLF = T))
message(sprintf("'xEnricher' is being called (\%s):",
as.character(now)), appendLF = T)
message(sprintf("#######################################################",
appendLF = T))
}
eTerm <- xEnricher(data = data, annotation = anno, g = g,
background = background, size.range = size.range, min.overlap =
min.overlap,
test = test, background.annotatable.only = background.annotatable.only,
p.tail = p.tail, p.adjust.method = p.adjust.method, ontology.algorithm
= "none",
true.path.rule = F, verbose = verbose)
if (verbose) {
now <- Sys.time()
message(sprintf("#######################################################",
appendLF = T))
message(sprintf("'xEnricher' has been finished (\%s)!",
as.character(now)), appendLF = T)
message(sprintf("#######################################################\n",
appendLF = T))
}
endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, "\%Y-\%m-\%d
\%H:\%M:\%S"),
strptime(startT, "\%Y-\%m-\%d \%H:\%M:\%S"), units = "secs"))
if (!silent) {
message(paste(c("\nEnd at ", as.character(endT)), collapse = ""),
appendLF = TRUE)
message(paste(c("Runtime in total (xEnricherYours): ",
runTime, " secs\n"), collapse = ""), appendLF = TRUE)
}
invisible(eTerm)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
