\name{Haploid_regulators}
\alias{Haploid_regulators}

\title{Haploid mutagenesis screens for regulators of protein phenotypes }

\usage{
data(Haploid_regulators)
}

\description{
This dataset contains regulators of 11 protein phenotypes identified by
haploid mutagenesis screens.
}

\value{
a data frame. It has the following columns: "Gene" (gene symbol;
regulators), "MI" (mutation index; negative values for positive
regulators and positive values for negative regulators), "FDR" and
"Phenotype" (one of 11 protein phenotypes).
}

\references{
Brockmann et al. (2017). Genetic wiring maps of single-cell protein
states reveal an off-switch for GPCR signalling. \emph{Nature},
546:307-11.
\cr
Mezzadra et al. (2017). Identification of CMTM6 and CMTM4 as PD-L1
protein regulators. \emph{Nature}, 549:106-10.
}

\keyword{datasets}

\examples{
\dontrun{
Haploid_regulators <- xRDataLoader('Haploid_regulators')
Haploid_regulators[1:5,]
## for 'PDL1'
ind <- grepl('PDL1', Haploid_regulators$Phenotype)
df <- Haploid_regulators[ind,c('Gene','MI','FDR')]
}
}
