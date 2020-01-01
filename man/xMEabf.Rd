\name{xMEabf}
\alias{xMEabf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xMEabf(eqtl.summary, gwas.summary, prior.eqtl = 1e-04, prior.gwas =
1e-04, prior.both = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{eqtl.summary}{
%% ~~Describe \code{eqtl.summary} here~~
}
\item{gwas.summary}{
%% ~~Describe \code{gwas.summary} here~~
}
\item{prior.eqtl}{
%% ~~Describe \code{prior.eqtl} here~~
}
\item{prior.gwas}{
%% ~~Describe \code{prior.gwas} here~~
}
\item{prior.both}{
%% ~~Describe \code{prior.both} here~~
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
function (eqtl.summary, gwas.summary, prior.eqtl = 1e-04, prior.gwas =
1e-04,
prior.both = 1e-05)
{
if (all(c("beta", "varbeta", "N", "MAF", "snp") \%in\%
names(eqtl.summary))) {
oneover <- 1/eqtl.summary$MAF
nvx <- 2 * eqtl.summary$N * eqtl.summary$MAF * (1 - eqtl.summary$MAF)
m <- stats::lm(nvx ~ oneover - 1)
cf <- stats::coef(m)[["oneover"]]
if (cf > 0) {
eqtl.summary$sdY <- sqrt(cf)
z <- eqtl.summary$beta/sqrt(eqtl.summary$varbeta)
V <- eqtl.summary$varbeta
sd.prior <- 0.15 * eqtl.summary$sdY
r <- sd.prior^2/(sd.prior^2 + V)
lABF <- 0.5 * (log(1 - r) + (r * z^2))
df_eQTL <- data.frame(eqtl.V = V, eqtl.z = z, eqtl.r = r,
eqtl.lABF = lABF, snp = eqtl.summary$snp, stringsAsFactors = F)
}
else {
return(NULL)
}
}
else {
return(NULL)
}
if (all(c("beta", "varbeta", "snp") \%in\% names(gwas.summary))) {
z <- gwas.summary$beta/sqrt(gwas.summary$varbeta)
V <- gwas.summary$varbeta
sd.prior <- 0.2
r <- sd.prior^2/(sd.prior^2 + V)
lABF <- 0.5 * (log(1 - r) + (r * z^2))
df_gwas <- data.frame(gwas.V = V, gwas.z = z, gwas.r = r,
gwas.lABF = lABF, snp = gwas.summary$snp, stringsAsFactors = F)
}
else {
return(NULL)
}
merged.df <- merge(df_eQTL, df_gwas)
if (nrow(merged.df) == 0) {
warnings("fail to merge df_eQTL and df_gwas. Both should contain the
common snps")
return(NULL)
}
else {
merged.df$both.sum.lABF <- merged.df$eqtl.lABF + merged.df$gwas.lABF
logsum <- function(x) {
max(x) + log(sum(exp(x - max(x))))
}
logdiff <- function(x, y) {
max(x, y) + log(exp(x - max(x, y)) - exp(y - max(x,
y)))
}
x <- merged.df$both.sum.lABF
my.denom.log.abf <- logsum(x)
merged.df$SNP.PP.H4 <- exp(x - my.denom.log.abf)
SNP.PP.H4 <- NULL
merged.df <- merged.df \%>\% dplyr::arrange(-SNP.PP.H4)
p1 <- prior.eqtl
p2 <- prior.gwas
p12 <- prior.both
l1 <- merged.df$eqtl.lABF
l2 <- merged.df$gwas.lABF
lsum <- l1 + l2
lH0.abf <- 0
lH1.abf <- log(p1) + logsum(l1)
lH2.abf <- log(p2) + logsum(l2)
lH3.abf <- log(p1) + log(p2) + logdiff(logsum(l1) + logsum(l2),
logsum(lsum))
lH4.abf <- log(p12) + logsum(lsum)
all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
pp.abf <- exp(all.abf - logsum(all.abf))
names(pp.abf) <- paste0("PP.H", (1:length(pp.abf)) -
1, ".abf")
output <- list(summary = c(nsnps = nrow(merged.df), pp.abf),
results = merged.df)
invisible(output)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
