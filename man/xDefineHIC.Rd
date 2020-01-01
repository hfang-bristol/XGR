\name{xDefineHIC}
\alias{xDefineHIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xDefineHIC(data = NULL, entity = c("SNP", "chr:start-end",
"data.frame", "bed", "GRanges"), include.HiC = c(NA, "Monocytes",
"Macrophages_M0", "Macrophages_M1", "Macrophages_M2", "Neutrophils",
"Megakaryocytes", "Endothelial_precursors", "Erythroblasts",
"Fetal_thymus", "Naive_CD4_T_cells", "Total_CD4_T_cells",
"Activated_total_CD4_T_cells", "Nonactivated_total_CD4_T_cells",
"Naive_CD8_T_cells", "Total_CD8_T_cells", "Naive_B_cells",
"Total_B_cells", "PE.Monocytes", "PE.Macrophages_M0",
"PE.Macrophages_M1", "PE.Macrophages_M2", "PE.Neutrophils",
"PE.Megakaryocytes", "PE.Erythroblasts", "PE.Naive_CD4_T_cells",
"PE.Naive_CD8_T_cells", "Combined", "Combined_PE"), GR.SNP =
c("dbSNP_GWAS", "dbSNP_Common", "dbSNP_Single"), verbose = TRUE,
RData.location = "http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{entity}{
%% ~~Describe \code{entity} here~~
}
\item{include.HiC}{
%% ~~Describe \code{include.HiC} here~~
}
\item{GR.SNP}{
%% ~~Describe \code{GR.SNP} here~~
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
function (data = NULL, entity = c("SNP", "chr:start-end", "data.frame",
"bed", "GRanges"), include.HiC = c(NA, "Monocytes", "Macrophages_M0",
"Macrophages_M1", "Macrophages_M2", "Neutrophils", "Megakaryocytes",
"Endothelial_precursors", "Erythroblasts", "Fetal_thymus",
"Naive_CD4_T_cells", "Total_CD4_T_cells",
"Activated_total_CD4_T_cells",
"Nonactivated_total_CD4_T_cells", "Naive_CD8_T_cells",
"Total_CD8_T_cells",
"Naive_B_cells", "Total_B_cells", "PE.Monocytes", "PE.Macrophages_M0",
"PE.Macrophages_M1", "PE.Macrophages_M2", "PE.Neutrophils",
"PE.Megakaryocytes", "PE.Erythroblasts", "PE.Naive_CD4_T_cells",
"PE.Naive_CD8_T_cells", "Combined", "Combined_PE"), GR.SNP =
c("dbSNP_GWAS",
"dbSNP_Common", "dbSNP_Single"), verbose = TRUE, RData.location =
"http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
entity <- match.arg(entity)
default.include.HiC <- c("Monocytes", "Macrophages_M0",
"Macrophages_M1",
"Macrophages_M2", "Neutrophils", "Megakaryocytes",
"Endothelial_precursors",
"Erythroblasts", "Fetal_thymus", "Naive_CD4_T_cells",
"Total_CD4_T_cells", "Activated_total_CD4_T_cells",
"Nonactivated_total_CD4_T_cells",
"Naive_CD8_T_cells", "Total_CD8_T_cells", "Naive_B_cells",
"Total_B_cells", "PE.Monocytes", "PE.Macrophages_M0",
"PE.Macrophages_M1", "PE.Macrophages_M2", "PE.Neutrophils",
"PE.Megakaryocytes", "PE.Erythroblasts", "PE.Naive_CD4_T_cells",
"PE.Naive_CD8_T_cells", "Combined", "Combined_PE")
ind <- match(default.include.HiC, include.HiC)
include.HiC <- default.include.HiC[!is.na(ind)]
if (!is.null(data)) {
if (entity == "SNP") {
data_gr <- xSNPlocations(data, GR.SNP = GR.SNP, verbose = verbose,
RData.location = RData.location, guid = guid)
}
else {
data_gr <- xGR(data, format = entity, verbose = verbose,
RData.location = RData.location, guid = guid)
}
if (is.null(data_gr)) {
return(NULL)
}
}
df_returned <- NULL
if (length(include.HiC) > 0) {
res_list <- lapply(include.HiC, function(x) {
if (verbose) {
now <- Sys.time()
message(sprintf("Processing \%s ...", x), appendLF = TRUE)
}
if (x == "Combined") {
g <- xRDataLoader(RData.customised = "ig.PCHiC",
RData.location = RData.location, guid = guid,
verbose = verbose)
df <- do.call(cbind, igraph::edge_attr(g))
if (1) {
df[df < 5] <- NA
res <- xAggregate(log(df), verbose = FALSE)
vec_score <- res$Aggregate
}
else {
num <- apply(df >= 5, 1, sum)
dff <- df
dff[dff < 5] <- 0
total <- apply(dff, 1, sum)
ave <- log(total/num)
ave_scale <- (ave - min(ave))/(max(ave) - min(ave)) *
0.9999999
vec_score <- num + ave_scale
}
ig <- g
E(ig)$score <- vec_score
for (i in igraph::edge_attr_names(g)) {
ig <- igraph::delete_edge_attr(ig, i)
}
}
else if (x == "Combined_PE") {
g <- xRDataLoader(RData.customised = "ig.PCHiC_PE",
RData.location = RData.location, guid = guid,
verbose = verbose)
df <- do.call(cbind, igraph::edge_attr(g))
if (1) {
df[df < 5] <- NA
res <- xAggregate(log(df), verbose = FALSE)
vec_score <- res$Aggregate
}
else {
num <- apply(df >= 5, 1, sum)
dff <- df
dff[dff < 5] <- 0
total <- apply(dff, 1, sum)
ave <- log(total/num)
ave_scale <- (ave - min(ave))/(max(ave) - min(ave)) *
0.9999999
vec_score <- num + ave_scale
}
ig <- g
E(ig)$score <- vec_score
for (i in igraph::edge_attr_names(g)) {
ig <- igraph::delete_edge_attr(ig, i)
}
}
else {
if (sum(grep("^PE.", x, perl = TRUE)) > 0) {
RData.customised <- paste("ig.PCHiC_", x, sep = "")
}
else {
RData.customised <- paste("ig.PCHiC.", x, sep = "")
}
ig <- xRDataLoader(RData.customised = RData.customised,
RData.location = RData.location, guid = guid,
verbose = verbose)
}
df_nodes <- igraph::get.data.frame(ig, what = "vertices")
df_edges <- igraph::get.data.frame(ig, what = "edges")
if (!is.null(data)) {
nodes_gr <- xGR(data = df_nodes[, 1], format = "chr:start-end",
verbose = verbose, RData.location = RData.location,
guid = guid)
maxgap <- -1L
minoverlap <- 0L
subject <- nodes_gr
query <- data_gr
q2r <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query =
query,
subject = subject, maxgap = maxgap, minoverlap = minoverlap,
type = "any", select = "all", ignore.strand = TRUE)))
res_df <- data.frame(SNP = names(data_gr)[q2r[,
1]], nodes = names(nodes_gr)[q2r[, 2]], stringsAsFactors = FALSE)
ind_from <- match(res_df$nodes, df_edges[, "from"])
nodes_from <- res_df[!is.na(ind_from), ]
partner_from <- df_edges[ind_from[!is.na(ind_from)],
c("to", "score")]
partner_from_gene <- df_nodes[match(partner_from[,
1], df_nodes[, 1]), 2]
nodes_from_gene <- df_nodes[match(nodes_from[,
2], df_nodes[, 1]), 2]
df_from <- cbind(nodes_from, partner_from, nodes_from_gene,
partner_from_gene, stringsAsFactors = FALSE)
colnames(df_from) <- c("SNP", "harbor", "partner",
"score", "harbor_genes", "partner_genes")
df_from$harbor_end <- rep("bait/from", nrow(df_from))
ind_to <- match(res_df$nodes, df_edges[, "to"])
nodes_to <- res_df[!is.na(ind_to), ]
partner_to <- df_edges[ind_to[!is.na(ind_to)],
c("from", "score")]
partner_to_gene <- df_nodes[match(partner_to[,
1], df_nodes[, 1]), 2]
nodes_to_gene <- df_nodes[match(nodes_to[, 2],
df_nodes[, 1]), 2]
df_to <- cbind(nodes_to, partner_to, nodes_to_gene,
partner_to_gene, stringsAsFactors = FALSE)
colnames(df_to) <- c("SNP", "harbor", "partner",
"score", "harbor_genes", "partner_genes")
df_to$harbor_end <- rep("prey/to", nrow(df_to))
df <- rbind(df_from, df_to)
df$Context <- rep(x, nrow(df))
y <- df
ind <- which(y$harbor_end == "bait/from")
df_bait <- y[ind, c("harbor", "partner", "score",
"harbor_genes", "partner_genes", "SNP", "harbor_end",
"harbor", "Context")]
colnames(df_bait) <- c("from", "to", "score",
"from_genes", "to_genes", "SNP", "SNP_end",
"SNP_harbor", "Context")
ind <- which(y$harbor_end == "prey/to")
df_prey <- y[ind, c("partner", "harbor", "score",
"partner_genes", "harbor_genes", "SNP", "harbor_end",
"harbor", "Context")]
colnames(df_prey) <- c("from", "to", "score",
"from_genes", "to_genes", "SNP", "SNP_end",
"SNP_harbor", "Context")
df <- rbind(df_bait, df_prey)
rownames(df) <- NULL
}
else {
df <- df_edges
colnames(df) <- c("from", "to", "score")
df$Context <- rep(x, nrow(df))
}
return(df)
})
df_returned <- do.call(rbind, res_list)
}
if (!is.null(data) & !is.null(df_returned)) {
if (verbose) {
now <- Sys.time()
message(sprintf("Amongst \%d SNPs, \%d SNPs are falling into \%d
interacton regions",
length(unique(data)), length(unique(df_returned$SNP)),
length(unique(df_returned$SNP_harbor))), appendLF = TRUE)
}
context_ls <- split(x = df_returned, f = df_returned$Context)
ls_ig <- lapply(1:length(context_ls), function(i) {
df <- context_ls[[i]]
if (1) {
res_ls <- lapply(split(x = df$SNP, f = df$SNP_harbor),
function(x) {
paste(x, collapse = ";")
})
vec_harbor_SNPs <- unlist(res_ls)
relations <- df[, 1:3]
relations <- relations[!duplicated(relations),
]
nodes <- base::as.data.frame(rbind(as.matrix(df[,
c("from", "from_genes")]), as.matrix(df[, c("to",
"to_genes")])), stringsAsFactors = FALSE)
nodes <- nodes[!duplicated(nodes), ]
colnames(nodes) <- c("name", "target")
nodes$SNP <- NA
ind <- match(nodes$name, names(vec_harbor_SNPs))
nodes$SNP[!is.na(ind)] <- vec_harbor_SNPs[ind[!is.na(ind)]]
ig <- graph.data.frame(d = relations, directed = TRUE,
vertices = nodes)
class(ig) <- c("PCHiC", "igraph")
return(ig)
}
})
names(ls_ig) <- names(context_ls)
if (length(ls_ig) == 1) {
ls_ig <- ls_ig[[1]]
}
output <- list(df = df_returned, ig = ls_ig)
}
else {
output <- df_returned
}
invisible(output)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
