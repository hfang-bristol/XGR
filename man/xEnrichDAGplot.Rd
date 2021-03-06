% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xEnrichDAGplot.r
\name{xEnrichDAGplot}
\alias{xEnrichDAGplot}
\title{Function to visualise enrichment results using a direct acyclic graph (DAG)}
\usage{
xEnrichDAGplot(
eTerm,
top_num = 10,
ig = NULL,
displayBy = c("fc", "adjp", "fdr", "zscore", "pvalue"),
path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"),
height = 7,
width = 7,
margin = rep(0.1, 4),
colormap = c("yr", "bwr", "jet", "gbr", "wyr", "br", "rainbow", "wb",
"lightyellow-orange"),
ncolors = 40,
zlim = NULL,
colorbar = T,
colorbar.fraction = 0.1,
newpage = T,
layout.orientation = c("top_bottom", "left_right", "bottom_top",
"right_left"),
node.info = c("none", "term_id", "term_name", "both",
"full_term_name"),
wrap.width = NULL,
graph.node.attrs = NULL,
graph.edge.attrs = NULL,
node.attrs = NULL
)
}
\arguments{
\item{eTerm}{an object of class "eTerm"}

\item{top_num}{the number of the top terms (sorted according to FDR or
adjusted p-values). If it is 'auto', only the significant terms (FDR <
0.05) will be displayed}

\item{ig}{the igraph object. If provided, only those terms within it
will be visualised. By default, it is NULL meaning no surch
restriction}

\item{displayBy}{which statistics will be used for displaying. It can
be "fc" for enrichment fold change (by default), "adjp" or "fdr" for
adjusted p value (or FDR), "pvalue" for p value, "zscore" for
enrichment z-score}

\item{path.mode}{the mode of paths induced by nodes in query. It can be
"all_paths" for all possible paths to the root, "shortest_paths" for
only one path to the root (for each node in query),
"all_shortest_paths" for all shortest paths to the root (i.e. for each
node, find all shortest paths with the equal lengths)}

\item{height}{a numeric value specifying the height of device}

\item{width}{a numeric value specifying the width of device}

\item{margin}{margins as units of length 4 or 1}

\item{colormap}{short name for the colormap. It can be one of "jet"
(jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red
colormap), "wyr" (white-yellow-red colormap), "br" (black-red
colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and
"rainbow" (rainbow colormap, that is,
red-yellow-green-cyan-blue-magenta). Alternatively, any
hyphen-separated HTML color names, e.g. "lightyellow-orange" (by
default), "blue-black-yellow", "royalblue-white-sandybrown",
"darkgreen-white-darkviolet". A list of standard color names can be
found in \url{http://html-color-codes.info/color-names}}

\item{ncolors}{the number of colors specified over the colormap}

\item{zlim}{the minimum and maximum z/data values for which colors
should be plotted, defaulting to the range of the finite values of z.
Each of the given colors will be used to color an equispaced interval
of this range. The midpoints of the intervals cover the range, so that
values just outside the range will be plotted}

\item{colorbar}{logical to indicate whether to append a colorbar. If
data is null, it always sets to false}

\item{colorbar.fraction}{the relative fraction of colorbar block
against the device size}

\item{newpage}{logical to indicate whether to open a new page. By
default, it sets to true for opening a new page}

\item{layout.orientation}{the orientation of the DAG layout. It can be
one of "left_right" for the left-right layout (viewed from the DAG root
point), "top_bottom" for the top-bottom layout, "bottom_top" for the
bottom-top layout, and "right_left" for the right-left layout}

\item{node.info}{tells the ontology term information used to label
nodes. It can be one of "none" for no node labeling, "term_id" for
using Term ID, "term_name" for using Term Name (the first 15
characters), "both" for using both of Term ID and Name (the first 15
characters), and "full_term_name" for using the full Term Name}

\item{wrap.width}{a positive integer specifying wrap width of Term
Name}

\item{graph.node.attrs}{a list of global node attributes. These node
attributes will be changed globally. See 'Note' below for details on
the attributes}

\item{graph.edge.attrs}{a list of global edge attributes. These edge
attributes will be changed globally. See 'Note' below for details on
the attributes}

\item{node.attrs}{a list of local edge attributes. These node
attributes will be changed locally; as such, for each attribute, the
input value must be a named vector (i.e. using Term ID as names). See
'Note' below for details on the attributes}
}
\value{
An object of class 'Ragraph'
}
\description{
\code{xEnrichDAGplot} is supposed to visualise enrichment results using
a direct acyclic graph (DAG) with node colorings. By default,
significant terms (of interest) are highlighted by box-shaped nodes,
the others by ellipse nodes. It returns an object of class 'Ragraph'.
}
\note{
A list of global node attributes used in "graph.node.attrs":
\itemize{
\item{"shape": the shape of the node: "circle", "rectangle", "rect",
"box" and "ellipse"}
\item{"fixedsize": the logical to use only width and height attributes.
By default, it sets to true for not expanding for the width of the
label}
\item{"fillcolor": the background color of the node}
\item{"color": the color for the node, corresponding to the outside
edge of the node}
\item{"fontcolor": the color for the node text/labelings}
\item{"fontsize": the font size for the node text/labelings}
\item{"height": the height (in inches) of the node: 0.5 by default}
\item{"width": the width (in inches) of the node: 0.75 by default}
\item{"style": the line style for the node: "solid", "dashed",
"dotted", "invis" and "bold"}
}
A list of global edge attributes used in "graph.edge.attrs":
\itemize{
\item{"color": the color of the edge: gray by default}
\item{"weight": the weight of the edge: 1 by default}
\item{"style": the line style for the edge: "solid", "dashed",
"dotted", "invis" and "bold"}
}
A list of local node attributes used in "node.attrs" (only those named
Term IDs will be changed locally!):
\itemize{
\item{"label": a named vector specifying the node text/labelings}
\item{"shape": a named vector specifying the shape of the node:
"circle", "rectangle", "rect", "box" and "ellipse"}
\item{"fixedsize": a named vector specifying whether it sets to true
for not expanding for the width of the label}
\item{"fillcolor": a named vector specifying the background color of
the node}
\item{"color": a named vector specifying the color for the node,
corresponding to the outside edge of the node}
\item{"fontcolor": a named vector specifying the color for the node
text/labelings}
\item{"fontsize": a named vector specifying the font size for the node
text/labelings}
\item{"height": a named vector specifying the height (in inches) of the
node: 0.5 by default}
\item{"width": a named vector specifying the width (in inches) of the
node: 0.75 by default}
\item{"style": a named vector specifying the line style for the node:
"solid", "dashed", "dotted", "invis" and "bold"}
}
}
\examples{
\dontrun{
# Load the XGR package and specify the location of built-in data
library(XGR)
RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"

# 1) load eQTL mapping results: cis-eQTLs significantly induced by IFN
cis <- xRDataLoader(RData.customised='JKscience_TS2A',
RData.location=RData.location)
ind <- which(cis$IFN_t > 0 & cis$IFN_fdr < 0.05)
df_cis <- cis[ind, c('variant','Symbol','IFN_t','IFN_fdr')]
data <- df_cis$variant

# 2) Enrichment analysis using Experimental Factor Ontology (EFO)
# Considering LD SNPs and respecting ontology tree
eTerm <- xEnricherSNPs(data, ontology="EF", include.LD="EUR",
LD.r2=0.8, ontology.algorithm="lea", RData.location=RData.location)

# 3) DAG plot of enrichment results
agDAG <- xEnrichDAGplot(eTerm, top_num="auto", displayBy="fc",
node.info=c("full_term_name"))
## modify node labels
xEnrichDAGplot(eTerm, top_num="auto", displayBy="fc",
node.info=c("full_term_name"),
graph.node.attrs=list(fontsize=25,fontcolor="blue",color="transparent"))
## modify node shapes
xEnrichDAGplot(eTerm, top_num="auto", displayBy="fc",
node.info=c("full_term_name"),
graph.node.attrs=list(fixedsize=FALSE,shape=c("ellipse","box","circle","plaintext")[2]))
## further modify edge color
xEnrichDAGplot(eTerm, top_num="auto", displayBy="fc",
node.info=c("full_term_name"), graph.node.attrs=list(fontsize=25),
graph.edge.attrs=list(color="black"))

# 4) hide labels for ellipse nodes
library(Rgraphviz)
name_nodes <- sapply(AgNode(agDAG), name)
shape_nodes <- sapply(AgNode(agDAG), shape)
names(shape_nodes) <- name_nodes
ind <- which(shape_nodes=='ellipse')
label_nodes <- rep('', length(ind))
names(label_nodes) <- name_nodes[ind]
xEnrichDAGplot(eTerm, top_num="auto", displayBy="fc",
node.info=c("full_term_name"), node.attrs=list(label=label_nodes,
shape=shape_nodes))

}
}
\seealso{
\code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}},
\code{\link{xEnrichViewer}}
}
