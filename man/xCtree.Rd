% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xCtree.r
\name{xCtree}
\alias{xCtree}
\title{Function to draw a tree-like circular plot}
\usage{
xCtree(
ig,
leave.label.orientation = c("outwards", "inwards"),
leave.label.size = 2,
leave.label.color = "steelblue",
leave.label.alpha = 0.7,
leave.label.wrap = NULL,
leave.label.expansion = NULL,
leave.size = 0,
limit.expansion = 1.1,
edge.color = "grey",
edge.alpha = 0.5,
edge.width = 0.5
)
}
\arguments{
\item{ig}{an object of class "igraph" with node attribute 'name'. It
could be a 'phylo' object converted to. Note: the leave labels would be
the node attribute 'name' unless the node attribute 'label' is
explicitely provided}

\item{leave.label.orientation}{the leave label orientation. It can be
"outwards" and "inwards"}

\item{leave.label.size}{the text size of the leave labelings. By
default, it is 2}

\item{leave.label.color}{the color of the leave labelings}

\item{leave.label.alpha}{the alpha of the leave labelings}

\item{leave.label.wrap}{the wrap width of the leave labelings}

\item{leave.label.expansion}{the x- and y-expansion of the leave
labelings. The value of 1 for the exact location of the leave, and the
outwards (>1; by default 1.05 if NULL) and the inwards (<1; by default
0.98 if NULL)}

\item{leave.size}{the size of the leave nodes. By default, it is 0}

\item{limit.expansion}{the x- and y-limit expansion. By default, it is
1.1. Beware the orignial limit is [-1,1]}

\item{edge.color}{the color of edges}

\item{edge.alpha}{the alpha of edges}

\item{edge.width}{the width of edges}
}
\value{
a ggplot2 object appended with 'ig' and 'data' which should contain
columns 'x','y', 'leaf' (T/F), 'name' (the same as V(ig)$name), 'tipid'
(tip id), 'label' (if not given in ig, a 'name' varient), 'angle' and
'hjust' (assist in leave label orientation).
}
\description{
\code{xCtree} is supposed to draw a tree-like circular plot (dendrogram
circular layout), with tips labelled (outwards or inwards). The tree is
provided as an object of class "igraph".
}
\note{
none
}
\examples{
\dontrun{
library(XGR)
RData.location <- "http://galahad.well.ox.ac.uk/bigdata"

AA.template <- xRDataLoader("AA.template",
RData.location=RData.location)
# consensus tree
ig <- AA.template$consensus$ig

# outwards
gp <- xCtree(ig, leave.label.orientation="outwards",
leave.label.wrap=50, limit.expansion=1.5, leave.size=2)
head(gp$data \%>\% dplyr::arrange(tipid))

# inwards
gp <- xCtree(ig, leave.label.orientation="inwards",
leave.label.wrap=30)

# obtain 'xcoord' and 'ycoord'
gp <- ggraph::ggraph(ig, layout='dendrogram', circular=TRUE)
data <- gp$data \%>\% dplyr::arrange(.ggraph.orig_index)
V(ig)$xcoord <- data[,'x']
V(ig)$ycoord <- data[,'y']
}
}
\seealso{
\code{\link{xCtree}}
}
