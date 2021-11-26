\name{ImmunoBase}
\alias{ImmunoBase}

\title{Immune-disease associated variants, regions and genes from ImmunoBase (hg19)}

\usage{
data(ImmunoBase)
}

\description{
This dataset contains data obtained from ImmunoBase. For each of 20
immune-diseases, its associated variants, regions, and nearby genes
(within 500kb) are stored.
}

\value{
a list with 5 components:
\itemize{
\item{\code{disease}: a character of disease name}
\item{\code{variants}: an object of class "GRanges", storing genomic
locations of associated variants plus their significance and odd
ratios}
\item{\code{regions}: an object of class "GRanges", storing genomic
locations of associated regions}
\item{\code{genes_variants}: a named vector for nearby genes within
500kb of associated variants; the element names are gene symbols, the
element values for the shortest distance to all associated variants}
\item{\code{genes_regions}: a named vector for nearby genes within
500kb of associated regions; the element names are gene symbols, the
element values for the shortest distance to all associated regions}
}
}

\keyword{datasets}

\examples{
\dontrun{
ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase')
names(ImmunoBase)
ImmunoBase$AS$disease
ImmunoBase$AS$variants
head(ImmunoBase$AS$genes_variants)
head(ImmunoBase$AS$genes_regions)
}
}
