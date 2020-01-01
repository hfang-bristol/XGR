\name{xDefineEQTL}
\alias{xDefineEQTL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% ~~function to do ... ~~
}
\description{
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xDefineEQTL(data = NULL, include.eQTL = c(NA, "JKscience_CD14",
"JKscience_LPS2", "JKscience_LPS24", "JKscience_IFN", "JKscience_TS2A",
"JKscience_TS2A_CD14", "JKscience_TS2A_LPS2", "JKscience_TS2A_LPS24",
"JKscience_TS2A_IFN", "JKscience_TS2B", "JKscience_TS2B_CD14",
"JKscience_TS2B_LPS2", "JKscience_TS2B_LPS24", "JKscience_TS2B_IFN",
"JKscience_TS3A", "JKng_bcell", "JKng_bcell_cis", "JKng_bcell_trans",
"JKng_mono", "JKng_mono_cis", "JKng_mono_trans", "JKpg_CD4",
"JKpg_CD4_cis", "JKpg_CD4_trans", "JKpg_CD8", "JKpg_CD8_cis",
"JKpg_CD8_trans",
"JKnc_neutro", "JKnc_neutro_cis", "JKnc_neutro_trans",
"WESTRAng_blood", "WESTRAng_blood_cis", "WESTRAng_blood_trans",
"JK_nk", "JK_nk_cis", "JK_nk_trans", "GTEx_V4_Adipose_Subcutaneous",
"GTEx_V4_Artery_Aorta", "GTEx_V4_Artery_Tibial",
"GTEx_V4_Esophagus_Mucosa", "GTEx_V4_Esophagus_Muscularis",
"GTEx_V4_Heart_Left_Ventricle", "GTEx_V4_Lung",
"GTEx_V4_Muscle_Skeletal", "GTEx_V4_Nerve_Tibial",
"GTEx_V4_Skin_Sun_Exposed_Lower_leg", "GTEx_V4_Stomach",
"GTEx_V4_Thyroid", "GTEx_V4_Whole_Blood",
"GTEx_V6p_Adipose_Subcutaneous",
"GTEx_V6p_Adipose_Visceral_Omentum", "GTEx_V6p_Adrenal_Gland",
"GTEx_V6p_Artery_Aorta", "GTEx_V6p_Artery_Coronary",
"GTEx_V6p_Artery_Tibial",
"GTEx_V6p_Brain_Anterior_cingulate_cortex_BA24",
"GTEx_V6p_Brain_Caudate_basal_ganglia",
"GTEx_V6p_Brain_Cerebellar_Hemisphere", "GTEx_V6p_Brain_Cerebellum",
"GTEx_V6p_Brain_Cortex", "GTEx_V6p_Brain_Frontal_Cortex_BA9",
"GTEx_V6p_Brain_Hippocampus", "GTEx_V6p_Brain_Hypothalamus",
"GTEx_V6p_Brain_Nucleus_accumbens_basal_ganglia",
"GTEx_V6p_Brain_Putamen_basal_ganglia",
"GTEx_V6p_Breast_Mammary_Tissue",
"GTEx_V6p_Cells_EBVtransformed_lymphocytes",
"GTEx_V6p_Cells_Transformed_fibroblasts", "GTEx_V6p_Colon_Sigmoid",
"GTEx_V6p_Colon_Transverse",
"GTEx_V6p_Esophagus_Gastroesophageal_Junction",
"GTEx_V6p_Esophagus_Mucosa", "GTEx_V6p_Esophagus_Muscularis",
"GTEx_V6p_Heart_Atrial_Appendage", "GTEx_V6p_Heart_Left_Ventricle",
"GTEx_V6p_Liver", "GTEx_V6p_Lung", "GTEx_V6p_Muscle_Skeletal",
"GTEx_V6p_Nerve_Tibial", "GTEx_V6p_Ovary", "GTEx_V6p_Pancreas",
"GTEx_V6p_Pituitary",
"GTEx_V6p_Prostate", "GTEx_V6p_Skin_Not_Sun_Exposed_Suprapubic",
"GTEx_V6p_Skin_Sun_Exposed_Lower_leg",
"GTEx_V6p_Small_Intestine_Terminal_Ileum", "GTEx_V6p_Spleen",
"GTEx_V6p_Stomach", "GTEx_V6p_Testis", "GTEx_V6p_Thyroid",
"GTEx_V6p_Uterus", "GTEx_V6p_Vagina", "GTEx_V6p_Whole_Blood",
"eQTLGen", "eQTLGen_cis", "eQTLGen_trans", "scRNAseq_eQTL_Bcell",
"scRNAseq_eQTL_CD4", "scRNAseq_eQTL_CD8", "scRNAseq_eQTL_cMono",
"scRNAseq_eQTL_DC", "scRNAseq_eQTL_Mono", "scRNAseq_eQTL_ncMono",
"scRNAseq_eQTL_NK",
"scRNAseq_eQTL_PBMC", "jpRNAseq_eQTL_Bcell", "jpRNAseq_eQTL_CD4",
"jpRNAseq_eQTL_CD8", "jpRNAseq_eQTL_Mono", "jpRNAseq_eQTL_NK",
"jpRNAseq_eQTL_PBMC", "Pi_eQTL_Bcell", "Pi_eQTL_Blood", "Pi_eQTL_CD14",
"Pi_eQTL_CD4", "Pi_eQTL_CD8", "Pi_eQTL_IFN", "Pi_eQTL_LPS2",
"Pi_eQTL_LPS24", "Pi_eQTL_Monocyte", "Pi_eQTL_Neutrophil",
"Pi_eQTL_NK", "Pi_eQTL_shared_CD14", "Pi_eQTL_shared_IFN",
"Pi_eQTL_shared_LPS2", "Pi_eQTL_shared_LPS24", "Osteoblast_eQTL"),
eQTL.customised = NULL, verbose = TRUE, RData.location =
"http://galahad.well.ox.ac.uk/bigdata", guid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
%% ~~Describe \code{data} here~~
}
\item{include.eQTL}{
%% ~~Describe \code{include.eQTL} here~~
}
\item{eQTL.customised}{
%% ~~Describe \code{eQTL.customised} here~~
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
function (data = NULL, include.eQTL = c(NA, "JKscience_CD14",
"JKscience_LPS2", "JKscience_LPS24", "JKscience_IFN", "JKscience_TS2A",
"JKscience_TS2A_CD14", "JKscience_TS2A_LPS2", "JKscience_TS2A_LPS24",
"JKscience_TS2A_IFN", "JKscience_TS2B", "JKscience_TS2B_CD14",
"JKscience_TS2B_LPS2", "JKscience_TS2B_LPS24", "JKscience_TS2B_IFN",
"JKscience_TS3A", "JKng_bcell", "JKng_bcell_cis", "JKng_bcell_trans",
"JKng_mono", "JKng_mono_cis", "JKng_mono_trans", "JKpg_CD4",
"JKpg_CD4_cis", "JKpg_CD4_trans", "JKpg_CD8", "JKpg_CD8_cis",
"JKpg_CD8_trans", "JKnc_neutro", "JKnc_neutro_cis",
"JKnc_neutro_trans",
"WESTRAng_blood", "WESTRAng_blood_cis", "WESTRAng_blood_trans",
"JK_nk", "JK_nk_cis", "JK_nk_trans", "GTEx_V4_Adipose_Subcutaneous",
"GTEx_V4_Artery_Aorta", "GTEx_V4_Artery_Tibial",
"GTEx_V4_Esophagus_Mucosa",
"GTEx_V4_Esophagus_Muscularis", "GTEx_V4_Heart_Left_Ventricle",
"GTEx_V4_Lung", "GTEx_V4_Muscle_Skeletal", "GTEx_V4_Nerve_Tibial",
"GTEx_V4_Skin_Sun_Exposed_Lower_leg", "GTEx_V4_Stomach",
"GTEx_V4_Thyroid", "GTEx_V4_Whole_Blood",
"GTEx_V6p_Adipose_Subcutaneous",
"GTEx_V6p_Adipose_Visceral_Omentum", "GTEx_V6p_Adrenal_Gland",
"GTEx_V6p_Artery_Aorta", "GTEx_V6p_Artery_Coronary",
"GTEx_V6p_Artery_Tibial",
"GTEx_V6p_Brain_Anterior_cingulate_cortex_BA24",
"GTEx_V6p_Brain_Caudate_basal_ganglia",
"GTEx_V6p_Brain_Cerebellar_Hemisphere", "GTEx_V6p_Brain_Cerebellum",
"GTEx_V6p_Brain_Cortex", "GTEx_V6p_Brain_Frontal_Cortex_BA9",
"GTEx_V6p_Brain_Hippocampus", "GTEx_V6p_Brain_Hypothalamus",
"GTEx_V6p_Brain_Nucleus_accumbens_basal_ganglia",
"GTEx_V6p_Brain_Putamen_basal_ganglia",
"GTEx_V6p_Breast_Mammary_Tissue",
"GTEx_V6p_Cells_EBVtransformed_lymphocytes",
"GTEx_V6p_Cells_Transformed_fibroblasts", "GTEx_V6p_Colon_Sigmoid",
"GTEx_V6p_Colon_Transverse",
"GTEx_V6p_Esophagus_Gastroesophageal_Junction",
"GTEx_V6p_Esophagus_Mucosa", "GTEx_V6p_Esophagus_Muscularis",
"GTEx_V6p_Heart_Atrial_Appendage", "GTEx_V6p_Heart_Left_Ventricle",
"GTEx_V6p_Liver", "GTEx_V6p_Lung", "GTEx_V6p_Muscle_Skeletal",
"GTEx_V6p_Nerve_Tibial", "GTEx_V6p_Ovary", "GTEx_V6p_Pancreas",
"GTEx_V6p_Pituitary", "GTEx_V6p_Prostate",
"GTEx_V6p_Skin_Not_Sun_Exposed_Suprapubic",
"GTEx_V6p_Skin_Sun_Exposed_Lower_leg",
"GTEx_V6p_Small_Intestine_Terminal_Ileum",
"GTEx_V6p_Spleen", "GTEx_V6p_Stomach", "GTEx_V6p_Testis",
"GTEx_V6p_Thyroid", "GTEx_V6p_Uterus", "GTEx_V6p_Vagina",
"GTEx_V6p_Whole_Blood", "eQTLGen", "eQTLGen_cis", "eQTLGen_trans",
"scRNAseq_eQTL_Bcell", "scRNAseq_eQTL_CD4", "scRNAseq_eQTL_CD8",
"scRNAseq_eQTL_cMono", "scRNAseq_eQTL_DC", "scRNAseq_eQTL_Mono",
"scRNAseq_eQTL_ncMono", "scRNAseq_eQTL_NK", "scRNAseq_eQTL_PBMC",
"jpRNAseq_eQTL_Bcell", "jpRNAseq_eQTL_CD4", "jpRNAseq_eQTL_CD8",
"jpRNAseq_eQTL_Mono", "jpRNAseq_eQTL_NK", "jpRNAseq_eQTL_PBMC",
"Pi_eQTL_Bcell", "Pi_eQTL_Blood", "Pi_eQTL_CD14", "Pi_eQTL_CD4",
"Pi_eQTL_CD8", "Pi_eQTL_IFN", "Pi_eQTL_LPS2", "Pi_eQTL_LPS24",
"Pi_eQTL_Monocyte", "Pi_eQTL_Neutrophil", "Pi_eQTL_NK",
"Pi_eQTL_shared_CD14",
"Pi_eQTL_shared_IFN", "Pi_eQTL_shared_LPS2", "Pi_eQTL_shared_LPS24",
"Osteoblast_eQTL"), eQTL.customised = NULL, verbose = TRUE,
RData.location = "http://galahad.well.ox.ac.uk/bigdata",
guid = NULL)
{
default.include.eQTL <- c("JKscience_CD14", "JKscience_LPS2",
"JKscience_LPS24", "JKscience_IFN", "JKscience_TS2A",
"JKscience_TS2A_CD14", "JKscience_TS2A_LPS2", "JKscience_TS2A_LPS24",
"JKscience_TS2A_IFN", "JKscience_TS2B", "JKscience_TS2B_CD14",
"JKscience_TS2B_LPS2", "JKscience_TS2B_LPS24", "JKscience_TS2B_IFN",
"JKscience_TS3A", "JKng_bcell", "JKng_bcell_cis", "JKng_bcell_trans",
"JKng_mono", "JKng_mono_cis", "JKng_mono_trans", "JKpg_CD4",
"JKpg_CD4_cis", "JKpg_CD4_trans", "JKpg_CD8", "JKpg_CD8_cis",
"JKpg_CD8_trans", "JKnc_neutro", "JKnc_neutro_cis",
"JKnc_neutro_trans",
"WESTRAng_blood", "WESTRAng_blood_cis", "WESTRAng_blood_trans",
"JK_nk", "JK_nk_cis", "JK_nk_trans", "GTEx_V4_Adipose_Subcutaneous",
"GTEx_V4_Artery_Aorta", "GTEx_V4_Artery_Tibial",
"GTEx_V4_Esophagus_Mucosa",
"GTEx_V4_Esophagus_Muscularis", "GTEx_V4_Heart_Left_Ventricle",
"GTEx_V4_Lung", "GTEx_V4_Muscle_Skeletal", "GTEx_V4_Nerve_Tibial",
"GTEx_V4_Skin_Sun_Exposed_Lower_leg", "GTEx_V4_Stomach",
"GTEx_V4_Thyroid", "GTEx_V4_Whole_Blood", "eQTLdb_NK",
"eQTLdb_CD14", "eQTLdb_LPS2", "eQTLdb_LPS24", "eQTLdb_IFN",
"GTEx_V6p_Adipose_Subcutaneous", "GTEx_V6p_Adipose_Visceral_Omentum",
"GTEx_V6p_Adrenal_Gland", "GTEx_V6p_Artery_Aorta",
"GTEx_V6p_Artery_Coronary",
"GTEx_V6p_Artery_Tibial",
"GTEx_V6p_Brain_Anterior_cingulate_cortex_BA24",
"GTEx_V6p_Brain_Caudate_basal_ganglia",
"GTEx_V6p_Brain_Cerebellar_Hemisphere",
"GTEx_V6p_Brain_Cerebellum", "GTEx_V6p_Brain_Cortex",
"GTEx_V6p_Brain_Frontal_Cortex_BA9", "GTEx_V6p_Brain_Hippocampus",
"GTEx_V6p_Brain_Hypothalamus",
"GTEx_V6p_Brain_Nucleus_accumbens_basal_ganglia",
"GTEx_V6p_Brain_Putamen_basal_ganglia",
"GTEx_V6p_Breast_Mammary_Tissue",
"GTEx_V6p_Cells_EBVtransformed_lymphocytes",
"GTEx_V6p_Cells_Transformed_fibroblasts",
"GTEx_V6p_Colon_Sigmoid", "GTEx_V6p_Colon_Transverse",
"GTEx_V6p_Esophagus_Gastroesophageal_Junction",
"GTEx_V6p_Esophagus_Mucosa",
"GTEx_V6p_Esophagus_Muscularis", "GTEx_V6p_Heart_Atrial_Appendage",
"GTEx_V6p_Heart_Left_Ventricle", "GTEx_V6p_Liver", "GTEx_V6p_Lung",
"GTEx_V6p_Muscle_Skeletal", "GTEx_V6p_Nerve_Tibial",
"GTEx_V6p_Ovary", "GTEx_V6p_Pancreas", "GTEx_V6p_Pituitary",
"GTEx_V6p_Prostate", "GTEx_V6p_Skin_Not_Sun_Exposed_Suprapubic",
"GTEx_V6p_Skin_Sun_Exposed_Lower_leg",
"GTEx_V6p_Small_Intestine_Terminal_Ileum",
"GTEx_V6p_Spleen", "GTEx_V6p_Stomach", "GTEx_V6p_Testis",
"GTEx_V6p_Thyroid", "GTEx_V6p_Uterus", "GTEx_V6p_Vagina",
"GTEx_V6p_Whole_Blood", "eQTLGen", "eQTLGen_cis", "eQTLGen_trans",
"scRNAseq_eQTL_Bcell", "scRNAseq_eQTL_CD4", "scRNAseq_eQTL_CD8",
"scRNAseq_eQTL_cMono", "scRNAseq_eQTL_DC", "scRNAseq_eQTL_Mono",
"scRNAseq_eQTL_ncMono", "scRNAseq_eQTL_NK", "scRNAseq_eQTL_PBMC",
"jpRNAseq_eQTL_Bcell", "jpRNAseq_eQTL_CD4", "jpRNAseq_eQTL_CD8",
"jpRNAseq_eQTL_Mono", "jpRNAseq_eQTL_NK", "jpRNAseq_eQTL_PBMC",
"Pi_eQTL_Bcell", "Pi_eQTL_Blood", "Pi_eQTL_CD14", "Pi_eQTL_CD4",
"Pi_eQTL_CD8", "Pi_eQTL_IFN", "Pi_eQTL_LPS2", "Pi_eQTL_LPS24",
"Pi_eQTL_Monocyte", "Pi_eQTL_Neutrophil", "Pi_eQTL_NK",
"Pi_eQTL_shared_CD14", "Pi_eQTL_shared_IFN", "Pi_eQTL_shared_LPS2",
"Pi_eQTL_shared_LPS24", "Osteoblast_eQTL")
ind <- match(default.include.eQTL, include.eQTL)
include.eQTL <- default.include.eQTL[!is.na(ind)]
df_SGS <- NULL
if (length(include.eQTL) > 0) {
if (sum(grepl("GTEx_V4_", include.eQTL, perl = TRUE)) >
0) {
GTEx_V4 <- xRDataLoader(RData.customised = "GTEx_V4",
RData.location = RData.location, guid = guid,
verbose = verbose)
}
if (sum(grepl("GTEx_V6p_", include.eQTL, perl = TRUE)) >
0) {
GTEx_V6p <- xRDataLoader(RData.customised = "GTEx_V6p",
RData.location = RData.location, guid = guid,
verbose = verbose)
}
if (sum(grepl("scRNAseq_eQTL_", include.eQTL, perl = TRUE)) >
0) {
scRNAseq_eQTL <- xRDataLoader(RData.customised = "scRNAseq_eQTL",
RData.location = RData.location, guid = guid,
verbose = verbose)
}
if (sum(grepl("jpRNAseq_eQTL_", include.eQTL, perl = TRUE)) >
0) {
jpRNAseq_eQTL <- xRDataLoader(RData.customised = "jpRNAseq_eQTL",
RData.location = RData.location, guid = guid,
verbose = verbose)
}
res_list <- lapply(include.eQTL, function(x) {
if (verbose) {
now <- Sys.time()
message(sprintf("Processing \%s ...", x), appendLF = TRUE)
}
if (sum(grep("JKscience_TS2A", x, perl = TRUE)) >
0) {
cis <- xRDataLoader(RData.customised = "JKscience_TS2A",
RData.location = RData.location, guid = guid,
verbose = verbose)
if (x == "JKscience_TS2A") {
minFDR <- apply(cis[, c(9:12)], 1, base::min,
na.rm = TRUE)
df <- data.frame(SNP = cis[, 1], Gene = cis[,
4], Sig = minFDR, stringsAsFactors = FALSE)
}
else {
if (x == "JKscience_TS2A_CD14") {
j <- 9
}
else if (x == "JKscience_TS2A_LPS2") {
j <- 10
}
else if (x == "JKscience_TS2A_LPS24") {
j <- 11
}
else if (x == "JKscience_TS2A_IFN") {
j <- 12
}
ind <- which(!is.na(cis[, j]) & cis[, j] <
0.05)
df <- data.frame(SNP = cis[ind, 1], Gene = cis[ind,
4], Sig = cis[ind, j], stringsAsFactors = FALSE)
}
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("JKscience_TS2B", x, perl = TRUE)) >
0) {
cis <- xRDataLoader(RData.customised = "JKscience_TS2B",
RData.location = RData.location, guid = guid,
verbose = verbose)
if (x == "JKscience_TS2B") {
minFDR <- apply(cis[, c(9:12)], 1, base::min,
na.rm = TRUE)
df <- data.frame(SNP = cis[, 1], Gene = cis[,
4], Sig = minFDR, stringsAsFactors = FALSE)
}
else {
if (x == "JKscience_TS2B_CD14") {
j <- 9
}
else if (x == "JKscience_TS2B_LPS2") {
j <- 10
}
else if (x == "JKscience_TS2B_LPS24") {
j <- 11
}
else if (x == "JKscience_TS2B_IFN") {
j <- 12
}
ind <- which(!is.na(cis[, j]) & cis[, j] <
0.05)
df <- data.frame(SNP = cis[ind, 1], Gene = cis[ind,
4], Sig = cis[ind, j], stringsAsFactors = FALSE)
}
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (x == "JKscience_TS3A") {
trans <- xRDataLoader(RData.customised = "JKscience_TS3A",
RData.location = RData.location, guid = guid,
verbose = verbose)
minFDR <- apply(trans[, c(9:12)], 1, base::min,
na.rm = TRUE)
df <- data.frame(SNP = trans[, 1], Gene = trans[,
4], Sig = minFDR, stringsAsFactors = FALSE)
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (x == "JKscience_CD14" | x == "JKscience_LPS2" |
x == "JKscience_LPS24" | x == "JKscience_IFN") {
cis <- xRDataLoader(RData.customised = "JKscience_TS2A",
RData.location = RData.location, guid = guid,
verbose = verbose)
trans <- xRDataLoader(RData.customised = "JKscience_TS3A",
RData.location = RData.location, guid = guid,
verbose = verbose)
df <- rbind(cis, trans)
if (x == "JKscience_CD14") {
j <- 9
}
else if (x == "JKscience_LPS2") {
j <- 10
}
else if (x == "JKscience_LPS24") {
j <- 11
}
else if (x == "JKscience_IFN") {
j <- 12
}
ind <- which(!is.na(df[, j]) & df[, j] < 0.05)
df <- data.frame(SNP = df[ind, 1], Gene = df[ind,
4], Sig = df[ind, j], stringsAsFactors = FALSE)
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("JKng_bcell", x, perl = TRUE)) >
0) {
res_ls <- xRDataLoader(RData.customised = "JKng_bcell",
RData.location = RData.location, guid = guid,
verbose = verbose)
df_cis <- data.frame(SNP = res_ls$cis[, 1], Gene = res_ls$cis[,
2], Sig = res_ls$cis[, 5], stringsAsFactors = FALSE)
df_trans <- data.frame(SNP = res_ls$trans[, 1],
Gene = res_ls$trans[, 2], Sig = res_ls$trans[,
5], stringsAsFactors = FALSE)
if (x == "JKng_bcell") {
df <- rbind(df_cis, df_trans)
}
else if (x == "JKng_bcell_cis") {
df <- df_cis
}
else if (x == "JKng_bcell_trans") {
df <- df_trans
}
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("JKng_mono", x, perl = TRUE)) >
0) {
res_ls <- xRDataLoader(RData.customised = "JKng_mono",
RData.location = RData.location, guid = guid,
verbose = verbose)
df_cis <- data.frame(SNP = res_ls$cis[, 1], Gene = res_ls$cis[,
2], Sig = res_ls$cis[, 5], stringsAsFactors = FALSE)
df_trans <- data.frame(SNP = res_ls$trans[, 1],
Gene = res_ls$trans[, 2], Sig = res_ls$trans[,
5], stringsAsFactors = FALSE)
if (x == "JKng_mono") {
df <- rbind(df_cis, df_trans)
}
else if (x == "JKng_mono_cis") {
df <- df_cis
}
else if (x == "JKng_mono_trans") {
df <- df_trans
}
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("JKpg_CD4", x, perl = TRUE)) >
0) {
res_ls <- xRDataLoader(RData.customised = "JKpg_CD4",
RData.location = RData.location, guid = guid,
verbose = verbose)
df_cis <- data.frame(SNP = res_ls$cis[, 1], Gene = res_ls$cis[,
2], Sig = res_ls$cis[, 6], stringsAsFactors = FALSE)
df_trans <- data.frame(SNP = res_ls$trans[, 1],
Gene = res_ls$trans[, 2], Sig = res_ls$trans[,
6], stringsAsFactors = FALSE)
if (x == "JKpg_CD4") {
df <- rbind(df_cis, df_trans)
}
else if (x == "JKpg_CD4_cis") {
df <- df_cis
}
else if (x == "JKpg_CD4_trans") {
df <- df_trans
}
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("JKpg_CD8", x, perl = TRUE)) >
0) {
res_ls <- xRDataLoader(RData.customised = "JKpg_CD8",
RData.location = RData.location, guid = guid,
verbose = verbose)
df_cis <- data.frame(SNP = res_ls$cis[, 1], Gene = res_ls$cis[,
2], Sig = res_ls$cis[, 6], stringsAsFactors = FALSE)
df_trans <- data.frame(SNP = res_ls$trans[, 1],
Gene = res_ls$trans[, 2], Sig = res_ls$trans[,
6], stringsAsFactors = FALSE)
if (x == "JKpg_CD8") {
df <- rbind(df_cis, df_trans)
}
else if (x == "JKpg_CD8_cis") {
df <- df_cis
}
else if (x == "JKpg_CD8_trans") {
df <- df_trans
}
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("JKnc_neutro", x, perl = TRUE)) >
0) {
res_ls <- xRDataLoader(RData.customised = "JKnc_neutro",
RData.location = RData.location, guid = guid,
verbose = verbose)
df_cis <- data.frame(SNP = res_ls$cis[, 1], Gene = res_ls$cis[,
2], Sig = res_ls$cis[, 6], stringsAsFactors = FALSE)
df_trans <- data.frame(SNP = res_ls$trans[, 1],
Gene = res_ls$trans[, 2], Sig = res_ls$trans[,
6], stringsAsFactors = FALSE)
if (x == "JKnc_neutro") {
df <- rbind(df_cis, df_trans)
}
else if (x == "JKnc_neutro_cis") {
df <- df_cis
}
else if (x == "JKnc_neutro_trans") {
df <- df_trans
}
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("WESTRAng_blood", x, perl = TRUE)) >
0) {
res_ls <- xRDataLoader(RData.customised = "WESTRAng_blood",
RData.location = RData.location, guid = guid,
verbose = verbose)
df_cis <- data.frame(SNP = res_ls$cis[, 1], Gene = res_ls$cis[,
2], Sig = res_ls$cis[, 6], stringsAsFactors = FALSE)
df_trans <- data.frame(SNP = res_ls$trans[, 1],
Gene = res_ls$trans[, 2], Sig = res_ls$trans[,
6], stringsAsFactors = FALSE)
if (x == "WESTRAng_blood") {
df <- rbind(df_cis, df_trans)
}
else if (x == "WESTRAng_blood_cis") {
df <- df_cis
}
else if (x == "WESTRAng_blood_trans") {
df <- df_trans
}
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("eQTLGen", x, perl = TRUE)) > 0) {
res_ls <- xRDataLoader(RData.customised = "eQTLGen",
RData.location = RData.location, guid = guid,
verbose = verbose)
df_cis <- data.frame(SNP = res_ls$cis[, 1], Gene = res_ls$cis[,
2], Sig = res_ls$cis[, 3], stringsAsFactors = FALSE)
df_trans <- data.frame(SNP = res_ls$trans[, 1],
Gene = res_ls$trans[, 2], Sig = res_ls$trans[,
3], stringsAsFactors = FALSE)
if (x == "eQTLGen") {
df <- rbind(df_cis, df_trans)
}
else if (x == "eQTLGen_cis") {
df <- df_cis
}
else if (x == "eQTLGen_trans") {
df <- df_trans
}
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("JK_nk", x, perl = TRUE)) > 0) {
res_ls <- xRDataLoader(RData.customised = "JK_nk",
RData.location = RData.location, guid = guid,
verbose = verbose)
df_cis <- data.frame(SNP = res_ls$cis[, 1], Gene = res_ls$cis[,
2], Sig = res_ls$cis[, 6], stringsAsFactors = FALSE)
df_trans <- data.frame(SNP = res_ls$trans[, 1],
Gene = res_ls$trans[, 2], Sig = res_ls$trans[,
6], stringsAsFactors = FALSE)
if (x == "JK_nk") {
df <- rbind(df_cis, df_trans)
}
else if (x == "JK_nk_cis") {
df <- df_cis
}
else if (x == "JK_nk_trans") {
df <- df_trans
}
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("GTEx_V4_", x, perl = TRUE)) >
0) {
y <- gsub("GTEx_V4_", "", x)
cis <- ""
eval(parse(text = paste("cis <- GTEx_V4$", y,
sep = "")))
df <- data.frame(SNP = cis[, 1], Gene = cis[,
2], Sig = cis[, 5], stringsAsFactors = FALSE)
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("GTEx_V6p_", x, perl = TRUE)) >
0) {
y <- gsub("GTEx_V6p_", "", x)
cis <- ""
eval(parse(text = paste("cis <- GTEx_V6p$", y,
sep = "")))
df <- data.frame(SNP = cis[, 1], Gene = cis[,
2], Sig = cis[, 5], stringsAsFactors = FALSE)
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("scRNAseq_eQTL_", x, perl = TRUE)) >
0) {
y <- gsub("scRNAseq_eQTL_", "", x)
cis <- ""
eval(parse(text = paste("cis <- scRNAseq_eQTL$",
y, sep = "")))
df <- data.frame(SNP = cis[, 1], Gene = cis[,
2], Sig = cis[, 5], stringsAsFactors = FALSE)
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("jpRNAseq_eQTL_", x, perl = TRUE)) >
0) {
y <- gsub("jpRNAseq_eQTL_", "", x)
cis <- ""
eval(parse(text = paste("cis <- jpRNAseq_eQTL$",
y, sep = "")))
df <- data.frame(SNP = cis[, 1], Gene = cis[,
2], Sig = cis[, 3], stringsAsFactors = FALSE)
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (x == "Osteoblast_eQTL") {
cis <- xRDataLoader(RData.customised = "Osteoblast_eQTL",
RData.location = RData.location, guid = guid,
verbose = verbose)
df <- data.frame(SNP = cis[, 2], Gene = cis[,
3], Sig = cis[, 7], stringsAsFactors = FALSE)
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("eQTLdb_", x, perl = TRUE)) > 0) {
cis <- xRDataLoader(RData.customised = x, RData.location =
RData.location,
guid = guid, verbose = verbose)
df <- data.frame(SNP = cis[, 1], Gene = cis[,
2], Sig = cis[, 5], stringsAsFactors = FALSE)
df <- cbind(df, Context = rep(x, nrow(df)), stringsAsFactors = FALSE)
}
else if (sum(grep("Pi_eQTL_", x, perl = TRUE)) >
0) {
res_ls <- xRDataLoader(RData.customised = x,
RData.location = RData.location, guid = guid,
verbose = verbose)
df_cis <- data.frame(SNP = res_ls$cis[, 1], Gene = res_ls$cis[,
2], Sig = res_ls$cis[, 6], stringsAsFactors = FALSE)
df_trans <- data.frame(SNP = res_ls$trans[, 1],
Gene = res_ls$trans[, 2], Sig = res_ls$trans[,
6], stringsAsFactors = FALSE)
df <- rbind(df_cis, df_trans)
df <- data.frame(df, Context = x, stringsAsFactors = FALSE)
}
else {
df <- NULL
}
return(df)
})
SGS <- do.call(rbind, res_list)
df_SGS <- SGS[!is.na(SGS[, 1]) & !is.na(SGS[, 2]), ]
if (verbose) {
now <- Sys.time()
message(sprintf("\%d eGenes are built-in", length(unique(df_SGS[,
2]))), appendLF = TRUE)
}
}
df_SGS_customised <- NULL
if (!is.null(eQTL.customised)) {
if (is.vector(eQTL.customised)) {
df <- utils::read.delim(file = eQTL.customised, header = TRUE,
row.names = NULL, stringsAsFactors = FALSE)
}
else if (is.matrix(eQTL.customised) | is.data.frame(eQTL.customised)) {
df <- eQTL.customised
}
if (!is.null(df)) {
colnames(df) <- c("SNP", "Gene", "Sig", "Context")
SGS_customised <- df
df_SGS_customised <- SGS_customised[!is.na(SGS_customised[,
1]) & !is.na(SGS_customised[, 2]), ]
if (verbose) {
now <- Sys.time()
message(sprintf("\%d eGenes are customised for \%d contexts",
length(unique(df_SGS_customised[, 2])),
length(unique(df_SGS_customised[,
4]))), appendLF = TRUE)
}
}
}
df_SGS <- do.call(rbind, list(df_SGS, df_SGS_customised))
if (!is.null(df_SGS)) {
df_SGS <- df_SGS[df_SGS[, 1] != "" & df_SGS[, 2] != "",
]
}
if (!is.null(data)) {
data <- gsub("_", ":", data, perl = TRUE)
data <- gsub("imm:", "chr", data, perl = TRUE)
data <- unique(data)
ind <- match(df_SGS[, 1], data)
df_SGS <- data.frame(df_SGS[!is.na(ind), ])
if (verbose) {
now <- Sys.time()
message(sprintf("A total of \%d input SNPs with \%d eGenes",
length(data), length(unique(df_SGS[, 2]))), appendLF = TRUE)
}
}
invisible(df_SGS)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
