\name{xSymbol2GeneID}
\alias{xSymbol2GeneID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xSymbol2GeneID(data, org = c("human", "mouse"), check.symbol.identity =
F, details = F, verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{org}{
%% ~~Describe \code{org} here~~
}
\item{check.symbol.identity}{
%% ~~Describe \code{check.symbol.identity} here~~
}
\item{details}{
%% ~~Describe \code{details} here~~
}
\item{verbose}{
%% ~~Describe \code{verbose} here~~
}
\item{RData.location}{
%% ~~Describe \code{RData.location} here~~
}
\item{guid}{
%% ~~Describe \code{guid} here~~
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
function (data, org = c("human", "mouse"), check.symbol.identity = F,
details = F, verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
if (!is.vector(data)) {
stop("The input data must be a vector.\n")
}
Symbol <- as.character(data)
if (class(org) == "EG") {
df_eg <- org$gene_info
if (verbose) {
message(sprintf("Customised organism (\%s)", as.character(Sys.time())),
appendLF = T)
}
}
else {
org <- org[1]
if (org == "human") {
df_eg <- xRDataLoader("org.Hs.eg", RData.location = RData.location,
guid = guid, verbose = verbose)$gene_info
}
else if (org == "mouse") {
df_eg <- xRDataLoader("org.Mm.eg", RData.location = RData.location,
guid = guid, verbose = verbose)$gene_info
}
if (verbose) {
message(sprintf("\%s organism (\%s)", org, as.character(Sys.time())),
appendLF = T)
}
}
type_of_gene <- ""
df_eg <- rbind(subset(df_eg, type_of_gene == "protein-coding"),
subset(df_eg, type_of_gene != "protein-coding"))
allGeneID <- df_eg$GeneID
allSymbol <- as.vector(df_eg$Symbol)
allSynonyms <- as.vector(df_eg$Synonyms)
allDescription <- as.vector(df_eg$description)
if (0) {
a <- Symbol
flag <- grep("-Mar$", a, ignore.case = T, perl = T, value = F)
if (length(flag) >= 1) {
b <- a[flag]
c <- sub("-Mar$", "", b, ignore.case = T, perl = T)
d <- sub("^0", "", c, ignore.case = T, perl = T)
e <- sapply(d, function(x) paste(c("March", x), collapse = ""))
a[flag] <- e
Symbol <- a
}
a <- Symbol
flag <- grep("-Sep$", a, ignore.case = T, perl = T, value = F)
if (length(flag) >= 1) {
b <- a[flag]
c <- sub("-Sep$", "", b, ignore.case = T, perl = T)
d <- sub("^0", "", c, ignore.case = T, perl = T)
e <- sapply(d, function(x) paste(c("Sept", x), collapse = ""))
a[flag] <- e
Symbol <- a
}
}
match_flag <- match(Symbol, allSymbol)
if (check.symbol.identity) {
na_flag <- is.na(match_flag)
a <- Symbol[na_flag]
tmp_flag <- is.na(match(tolower(allSymbol), tolower(Symbol)))
tmp_Synonyms <- allSynonyms[tmp_flag]
Orig.index <- seq(1, length(allSynonyms))
Orig.index <- Orig.index[tmp_flag]
b <- sapply(1:length(a), function(x) {
tmp_pattern1 <- paste("^", a[x], "\\|", sep = "")
tmp_pattern2 <- paste("\\|", a[x], "\\|", sep = "")
tmp_pattern3 <- paste("\\|", a[x], "$", sep = "")
tmp_pattern <- paste(tmp_pattern1, "|", tmp_pattern2,
"|", tmp_pattern3, sep = "")
tmp_result <- grep(tmp_pattern, tmp_Synonyms, ignore.case = T,
perl = T, value = F)
ifelse(length(tmp_result) == 1, Orig.index[tmp_result[1]],
NA)
})
match_flag[na_flag] <- b
if (verbose) {
now <- Sys.time()
message(sprintf("Among \%d symbols of input data, there are \%d
mappable via official gene symbols, \%d mappable via gene alias but \%d
left unmappable",
length(Symbol), (length(Symbol) - length(a)),
sum(!is.na(b)), sum(is.na(b))), appendLF = T)
}
}
else {
if (verbose) {
now <- Sys.time()
message(sprintf("Among \%d symbols of input data, there are \%d
mappable via official gene symbols but \%d left unmappable",
length(Symbol), (sum(!is.na(match_flag))), (sum(is.na(match_flag)))),
appendLF = T)
}
}
df_res <- df_eg[match_flag, ]
if (details) {
df_res <- data.frame(Input = Symbol, df_res, stringsAsFactors = F)
return(df_res)
}
else {
return(df_res$GeneID)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
