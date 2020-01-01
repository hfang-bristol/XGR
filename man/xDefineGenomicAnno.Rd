\name{xDefineGenomicAnno}
\alias{xDefineGenomicAnno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xDefineGenomicAnno(GR.annotation = c(NA, "Uniform_TFBS",
"ENCODE_TFBS_ClusteredV3", "ENCODE_TFBS_ClusteredV3_CellTypes",
"Uniform_DNaseI_HS", "ENCODE_DNaseI_ClusteredV3",
"ENCODE_DNaseI_ClusteredV3_CellTypes", "Broad_Histone", "SYDH_Histone",
"UW_Histone", "FANTOM5_Enhancer_Cell", "FANTOM5_Enhancer_Tissue",
"FANTOM5_Enhancer_Extensive", "FANTOM5_Enhancer",
"Segment_Combined_Gm12878", "Segment_Combined_H1hesc",
"Segment_Combined_Helas3", "Segment_Combined_Hepg2",
"Segment_Combined_Huvec", "Segment_Combined_K562", "TFBS_Conserved",
"TS_miRNA", "TCGA", "ReMap_Public_TFBS", "ReMap_Encode_TFBS",
"ReMap_PublicAndEncode_TFBS", "ReMap_Public_mergedTFBS",
"ReMap_Encode_mergedTFBS", "ReMap_PublicAndEncode_mergedTFBS",
"Blueprint_BoneMarrow_Histone", "Blueprint_CellLine_Histone",
"Blueprint_CordBlood_Histone", "Blueprint_Thymus_Histone",
"Blueprint_VenousBlood_Histone", "Blueprint_DNaseI",
"Blueprint_Methylation_hyper", "Blueprint_Methylation_hypo",
"EpigenomeAtlas_15Segments_E029", "EpigenomeAtlas_15Segments_E030",
"EpigenomeAtlas_15Segments_E031",
"EpigenomeAtlas_15Segments_E032", "EpigenomeAtlas_15Segments_E033",
"EpigenomeAtlas_15Segments_E034", "EpigenomeAtlas_15Segments_E035",
"EpigenomeAtlas_15Segments_E036", "EpigenomeAtlas_15Segments_E037",
"EpigenomeAtlas_15Segments_E038", "EpigenomeAtlas_15Segments_E039",
"EpigenomeAtlas_15Segments_E040", "EpigenomeAtlas_15Segments_E041",
"EpigenomeAtlas_15Segments_E042", "EpigenomeAtlas_15Segments_E043",
"EpigenomeAtlas_15Segments_E044", "EpigenomeAtlas_15Segments_E045",
"EpigenomeAtlas_15Segments_E046",
"EpigenomeAtlas_15Segments_E047", "EpigenomeAtlas_15Segments_E048",
"EpigenomeAtlas_15Segments_E050", "EpigenomeAtlas_15Segments_E051",
"EpigenomeAtlas_15Segments_E062", "CpG_anno", "Genic_anno",
"FANTOM5_CAT_Cell", "FANTOM5_CAT_Tissue", "FANTOM5_CAT_DO",
"FANTOM5_CAT_EFO", "FANTOM5_CAT_HPO", "FANTOM5_CAT_MESH",
"FANTOM5_CAT_PICS", "EpigenomeAtlas_DNaseNarrow",
"EpigenomeAtlas_DNaseBroad"), verbose = T, RData.location =
"http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{GR.annotation}{
%% ~~Describe \code{GR.annotation} here~~
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
function (GR.annotation = c(NA, "Uniform_TFBS",
"ENCODE_TFBS_ClusteredV3",
"ENCODE_TFBS_ClusteredV3_CellTypes", "Uniform_DNaseI_HS",
"ENCODE_DNaseI_ClusteredV3", "ENCODE_DNaseI_ClusteredV3_CellTypes",
"Broad_Histone", "SYDH_Histone", "UW_Histone", "FANTOM5_Enhancer_Cell",
"FANTOM5_Enhancer_Tissue", "FANTOM5_Enhancer_Extensive",
"FANTOM5_Enhancer", "Segment_Combined_Gm12878",
"Segment_Combined_H1hesc",
"Segment_Combined_Helas3", "Segment_Combined_Hepg2",
"Segment_Combined_Huvec",
"Segment_Combined_K562", "TFBS_Conserved", "TS_miRNA", "TCGA",
"ReMap_Public_TFBS", "ReMap_Encode_TFBS", "ReMap_PublicAndEncode_TFBS",
"ReMap_Public_mergedTFBS", "ReMap_Encode_mergedTFBS",
"ReMap_PublicAndEncode_mergedTFBS",
"Blueprint_BoneMarrow_Histone", "Blueprint_CellLine_Histone",
"Blueprint_CordBlood_Histone", "Blueprint_Thymus_Histone",
"Blueprint_VenousBlood_Histone", "Blueprint_DNaseI",
"Blueprint_Methylation_hyper",
"Blueprint_Methylation_hypo", "EpigenomeAtlas_15Segments_E029",
"EpigenomeAtlas_15Segments_E030", "EpigenomeAtlas_15Segments_E031",
"EpigenomeAtlas_15Segments_E032", "EpigenomeAtlas_15Segments_E033",
"EpigenomeAtlas_15Segments_E034", "EpigenomeAtlas_15Segments_E035",
"EpigenomeAtlas_15Segments_E036", "EpigenomeAtlas_15Segments_E037",
"EpigenomeAtlas_15Segments_E038", "EpigenomeAtlas_15Segments_E039",
"EpigenomeAtlas_15Segments_E040", "EpigenomeAtlas_15Segments_E041",
"EpigenomeAtlas_15Segments_E042", "EpigenomeAtlas_15Segments_E043",
"EpigenomeAtlas_15Segments_E044", "EpigenomeAtlas_15Segments_E045",
"EpigenomeAtlas_15Segments_E046", "EpigenomeAtlas_15Segments_E047",
"EpigenomeAtlas_15Segments_E048", "EpigenomeAtlas_15Segments_E050",
"EpigenomeAtlas_15Segments_E051", "EpigenomeAtlas_15Segments_E062",
"CpG_anno", "Genic_anno", "FANTOM5_CAT_Cell", "FANTOM5_CAT_Tissue",
"FANTOM5_CAT_DO", "FANTOM5_CAT_EFO", "FANTOM5_CAT_HPO",
"FANTOM5_CAT_MESH",
"FANTOM5_CAT_PICS", "EpigenomeAtlas_DNaseNarrow",
"EpigenomeAtlas_DNaseBroad"),
verbose = T, RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
grl <- NULL
if (class(GR.annotation) == "GenomicRangesList") {
grl <- GR.annotation
}
else {
lgr <- NULL
if (class(GR.annotation) == "list") {
lgr <- GR.annotation
if (!all(sapply(lgr, function(x) class(x) == "GRanges"))) {
lgr <- NULL
}
}
else if (class(GR.annotation) == "GRanges") {
lgr <- list(Customised = GR.annotation)
}
else {
if (length(GR.annotation) > 1) {
GR.annotation <- GR.annotation[1]
}
if (!is.na(GR.annotation)) {
lgr <- xRDataLoader(RData.customised = GR.annotation,
verbose = verbose, RData.location = RData.location,
guid = guid)
}
}
if (class(lgr) == "list") {
lgr <- base::Filter(base::Negate(is.null), lgr)
grl <- GenomicRanges::GRangesList(lgr)
}
}
invisible(grl)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
