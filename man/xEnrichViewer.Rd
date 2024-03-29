% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xEnrichViewer.r
\name{xEnrichViewer}
\alias{xEnrichViewer}
\title{Function to view enrichment results}
\usage{
xEnrichViewer(
eTerm,
top_num = 10,
sortBy = c("adjp", "fdr", "pvalue", "zscore", "fc", "nAnno",
"nOverlap", "or", "none"),
decreasing = NULL,
details = FALSE
)
}
\arguments{
\item{eTerm}{an object of class "eTerm"}

\item{top_num}{the number of the top terms (sorted according to
'sortBy' below) will be viewed}

\item{sortBy}{which statistics will be used for sorting and viewing
gene sets (terms). It can be "adjp" or "fdr" for adjusted p value
(FDR), "pvalue" for p value, "zscore" for enrichment z-score, "fc" for
enrichment fold change, "nAnno" for the number of sets (terms),
"nOverlap" for the number in overlaps, "or" for the odds ratio, and
"none" for ordering according to ID of terms}

\item{decreasing}{logical to indicate whether to sort in a decreasing
order. If it is null, it would be true for "zscore", "nAnno" or
"nOverlap"; otherwise it would be false}

\item{details}{logical to indicate whether the detailed information of
gene sets (terms) is also viewed. By default, it sets to false for no
inclusion}
}
\value{
a data frame with following components:
\itemize{
\item{\code{id}: term ID; as rownames}
\item{\code{name}: term name}
\item{\code{nAnno}: number in members annotated by a term}
\item{\code{nOverlap}: number in overlaps}
\item{\code{fc}: enrichment fold changes}
\item{\code{zscore}: enrichment z-score}
\item{\code{pvalue}: nominal p value}
\item{\code{adjp}: adjusted p value (FDR)}
\item{\code{or}: a vector containing odds ratio}
\item{\code{CIl}: a vector containing lower bound confidence interval
for the odds ratio}
\item{\code{CIu}: a vector containing upper bound confidence interval
for the odds ratio}
\item{\code{distance}: term distance or other information; optional, it
is only appended when "details" is true}
\item{\code{members_Overlap}: members (represented as Gene Symbols) in
overlaps; optional, it is only appended when "details" is true}
\item{\code{members_Anno}: members (represented as Gene Symbols) in
annotations; optional, it is only appended when "details" is true}
}
}
\description{
\code{xEnrichViewer} is supposed to view results of enrichment
analysis.
}
\note{
none
}
\examples{
RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
\dontrun{
xEnrichViewer(eTerm)
}
}
\seealso{
\code{\link{xEnrichViewer}}
}
