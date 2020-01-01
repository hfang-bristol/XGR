\name{xDefineOntology}
\alias{xDefineOntology}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xDefineOntology(ontology = c(NA, "GOBP", "GOMF", "GOCC", "PSG", "PS",
"PS2", "SF", "Pfam", "DO", "HPPA", "HPMI", "HPCM", "HPMA", "MP", "EF",
"MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP",
"MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT",
"MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF",
"MsigdbC5CC", "MsigdbC6", "MsigdbC7", "DGIdb", "GTExV4", "GTExV6p",
"GTExV7", "CreedsDisease", "CreedsDiseaseUP", "CreedsDiseaseDN",
"CreedsDrug", "CreedsDrugUP", "CreedsDrugDN",
"CreedsGene", "CreedsGeneUP", "CreedsGeneDN", "KEGG", "KEGGmetabolism",
"KEGGgenetic", "KEGGenvironmental", "KEGGcellular", "KEGGorganismal",
"KEGGdisease", "REACTOME", "REACTOME_ImmuneSystem",
"REACTOME_SignalTransduction", "CGL", "SIFTS2GOBP", "SIFTS2GOMF",
"SIFTS2GOCC", "EnrichrARCHS4Cells", "EnrichrARCHS4Tissues",
"EnrichrHumanGeneAtlas", "EnrichrTissueHumanProteomeMap",
"EnrichrTissueProteomicsDB", "EnrichrAchillesFitnessD",
"EnrichrAchillesFitnessI", "EnrichrDSigDB", "EnrichrOMIM",
"EnrichrOMIMexpanded",
"EnrichrdbGaP", "EnrichrJensenDiseases", "EnrichrJensenTissues",
"EnrichrBioCarta", "EnrichrKEGG", "EnrichrNCIpathway",
"EnrichrPanther", "EnrichrReactome", "EnrichrWikiPathways",
"EnrichrhuMAP", "EnrichrChEA", "EnrichrConsensusTFs",
"EnrichrEncodeTF", "EnrichrTFlof", "EnrichrTFpert"),
ontology.customised = NULL, anno.identity = c("GeneID", "Symbol"),
verbose = T, RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ontology}{
%% ~~Describe \code{ontology} here~~
}
\item{ontology.customised}{
%% ~~Describe \code{ontology.customised} here~~
}
\item{anno.identity}{
%% ~~Describe \code{anno.identity} here~~
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
function (ontology = c(NA, "GOBP", "GOMF", "GOCC", "PSG", "PS",
"PS2", "SF", "Pfam", "DO", "HPPA", "HPMI", "HPCM", "HPMA",
"MP", "EF", "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall",
"MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA",
"MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM",
"MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7",
"DGIdb", "GTExV4", "GTExV6p", "GTExV7", "CreedsDisease",
"CreedsDiseaseUP", "CreedsDiseaseDN", "CreedsDrug", "CreedsDrugUP",
"CreedsDrugDN", "CreedsGene", "CreedsGeneUP", "CreedsGeneDN",
"KEGG", "KEGGmetabolism", "KEGGgenetic", "KEGGenvironmental",
"KEGGcellular", "KEGGorganismal", "KEGGdisease", "REACTOME",
"REACTOME_ImmuneSystem", "REACTOME_SignalTransduction", "CGL",
"SIFTS2GOBP", "SIFTS2GOMF", "SIFTS2GOCC", "EnrichrARCHS4Cells",
"EnrichrARCHS4Tissues", "EnrichrHumanGeneAtlas",
"EnrichrTissueHumanProteomeMap",
"EnrichrTissueProteomicsDB", "EnrichrAchillesFitnessD",
"EnrichrAchillesFitnessI",
"EnrichrDSigDB", "EnrichrOMIM", "EnrichrOMIMexpanded", "EnrichrdbGaP",
"EnrichrJensenDiseases", "EnrichrJensenTissues", "EnrichrBioCarta",
"EnrichrKEGG", "EnrichrNCIpathway", "EnrichrPanther",
"EnrichrReactome",
"EnrichrWikiPathways", "EnrichrhuMAP", "EnrichrChEA",
"EnrichrConsensusTFs",
"EnrichrEncodeTF", "EnrichrTFlof", "EnrichrTFpert"),
ontology.customised = NULL,
anno.identity = c("GeneID", "Symbol"), verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
anno.identity <- match.arg(anno.identity)
ontology <- ontology[1]
ontology.input <- ontology
g <- NULL
anno <- NULL
if (!is.null(ontology.customised)) {
if (class(ontology.customised) == "GS") {
GS <- ontology.customised
anno <- GS$gs
nodes <- data.frame(name = as.character(GS$set_info$setID),
term_id = as.character(GS$set_info$setID), term_name =
as.character(GS$set_info$name),
term_distance = as.character(GS$set_info$distance),
term_namespace = as.character(GS$set_info$namespace),
stringsAsFactors = F)
nodes <- rbind(nodes, c("root", "root", "root", "root",
"root"))
relations <- data.frame(from = "root", to = nodes$name)
g <- igraph::graph.data.frame(d = relations, directed = T,
vertices = nodes)
g$type <- "iso"
ontology.input <- "customised"
}
}
else {
if (!is.na(ontology)) {
if (verbose) {
message(sprintf("Load the ontology \%s and its gene annotations (\%s)
...",
ontology, as.character(Sys.time())), appendLF = T)
}
flag_PS2 <- FALSE
if (ontology == "PS2") {
flag_PS2 <- TRUE
ontology <- "PS"
}
flag_REACTOME <- FALSE
if (grepl("REACTOME_", ontology)) {
flag_REACTOME <- TRUE
ontology_REACTOME <- ontology
ontology <- "REACTOME"
}
GS <- xRDataLoader(paste("org.Hs.eg", ontology, sep = ""),
RData.location = RData.location, guid = guid,
verbose = verbose)
if (flag_PS2) {
tmp <- as.character(unique(GS$set_info$name))
inds <- sapply(tmp, function(x) which(GS$set_info$name ==
x))
set_info <- data.frame()
for (i in 1:length(inds)) {
set_info <- rbind(set_info, as.matrix(GS$set_info[max(inds[[i]]),
]))
}
gs <- list()
for (i in 1:length(inds)) {
gs[[i]] <- unlist(GS$gs[inds[[i]]], use.names = F)
}
names(gs) <- rownames(set_info)
GS$set_info <- set_info
GS$gs <- gs
}
if (flag_REACTOME) {
flag <- unlist(strsplit(ontology_REACTOME, "_"))[2]
if (flag \%in\% GS$set_info$namespace) {
GS$set_info <- GS$set_info[GS$set_info$namespace ==
flag, ]
ind <- match(names(GS$gs), GS$set_info$setID)
GS$gs <- GS$gs[!is.na(ind)]
}
}
anno <- GS$gs
if (anno.identity == "Symbol") {
EG <- xRDataLoader(RData.customised = paste("org.Hs.eg",
sep = ""), RData.location = RData.location,
guid = guid, verbose = verbose)
anno_symbols <- lapply(anno, function(x) {
xGeneID2Symbol(x, org = EG, details = F, verbose = F)
})
anno <- anno_symbols
}
all.ontologies <- c("GOBP", "GOMF", "GOCC", "DO",
"HPPA", "HPMI", "HPCM", "HPMA", "MP", "EF", "SIFTS2GOBP",
"SIFTS2GOMF", "SIFTS2GOCC", "REACTOME")
flag_ontology <- ontology \%in\% all.ontologies
if (flag_ontology) {
ontology <- gsub("^SIFTS2", "", ontology)
g <- xRDataLoader(RData.customised = paste("ig.",
ontology, sep = ""), RData.location = RData.location,
guid = guid, verbose = verbose)
if (is.null(V(g)$term_namespace)) {
V(g)$term_namespace <- ontology
}
if (flag_REACTOME) {
if (0) {
g <- dnet::dDAGinduce(g, nodes_query = GS$set_info$setID,
path.mode = "all_paths")
}
else {
ind <- match(V(g)$name, GS$set_info$setID)
vids <- V(g)$name[!is.na(ind)]
neighs.in <- igraph::neighborhood(g, order = vcount(g),
nodes = vids, mode = "in")
neighbors <- unique(names(unlist(neighs.in)))
g <- dnet::dNetInduce(g, nodes_query = neighbors,
knn = 0, remove.loops = TRUE, largest.comp = T)
}
}
g$type <- "dag"
}
else {
nodes <- data.frame(name = as.character(GS$set_info$setID),
term_id = as.character(GS$set_info$setID),
term_name = as.character(GS$set_info$name),
term_distance = as.character(GS$set_info$distance),
term_namespace = as.character(GS$set_info$namespace),
stringsAsFactors = F)
nodes <- rbind(nodes, c("root", "root", "root",
"root", "root"))
relations <- data.frame(from = "root", to = nodes$name)
g <- igraph::graph.data.frame(d = relations,
directed = T, vertices = nodes)
g$type <- "iso"
}
}
}
g$ontology <- ontology.input
aOnto <- list(g = g, anno = anno)
class(aOnto) <- "aOnto"
invisible(aOnto)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
