#' Function to conduct region-based enrichment analysis using genomic annotations via sampling
#'
#' \code{xLDenricher} is supposed to conduct region-based enrichment analysis for the input genomic region data (genome build h19), using genomic annotations (eg active chromatin, transcription factor binding sites/motifs, conserved sites). Enrichment analysis is achieved by comparing the observed overlaps against the expected overlaps which are estimated from the null distribution. The null distribution is generated via sampling, that is, randomly generating samples for data genomic regions from background genomic regions. Background genomic regions can be provided by the user; by default, the annotatable genomic regions will be used. 
#'
#' @param bLD a bLD object, containing a set of blocks based on which to generate a null distribution
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 150) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 150) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised GR object directly
#' @param num.samples the number of samples randomly generated
#' @param respect how to respect the properties of to-be-sampled LD blocks. It can be one of 'maf' (respecting the maf of the best SNP), 'distance' (respecting the distance of the best SNP to the nearest gene), and 'both' (respecting the maf and distance)
#' @param restrict.chr logical to restrict to the same chromosome. By default, it sets to false
#' @param preserve how to preserve the resulting null LD block. It can be one of 'boundary' (preserving the boundary of the LD block), and 'exact' (exactly preserving the relative SNP locations within the LD block)
#' @param seed an integer specifying the seed
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param GR.annotation the genomic regions of annotation data. By default, it is 'NA' to disable this option. Pre-built genomic annotation data are detailed in the section 'Note'. Beyond pre-built annotation data, the user can specify the customised input. To do so, first save your RData file (a list of GR objects, each is an GR object correponding to an annotation) into your local computer. Then, tell "GR.annotation" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' a data frame with 13 columns:
#' \itemize{
#'  \item{\code{name}: the annotation name}
#'  \item{\code{nAnno}: the number of regions from annotation data}
#'  \item{\code{nOverlap}: the observed number of LD blocks overlapped with annotation data}
#'  \item{\code{fc}: fold change}
#'  \item{\code{zscore}: z-score}
#'  \item{\code{pvalue}: p-value}
#'  \item{\code{adjp}: adjusted p-value. It is the p value but after being adjusted for multiple comparisons}
#'  \item{\code{or}: a vector containing odds ratio}
#'  \item{\code{CIl}: a vector containing lower bound confidence interval for the odds ratio}
#'  \item{\code{CIu}: a vector containing upper bound confidence interval for the odds ratio}
#'  \item{\code{nData}: the number of input LD blocks}
#'  \item{\code{nExpect}: the expected number of LD blocks overlapped with annotation data}
#'  \item{\code{std}: the standard deviation of expected number of LD blocks overlapped with annotation data}
#' }
#' @note The genomic annotation data are described below according to the data sources and data types.\cr
#' 1. ENCODE Transcription Factor ChIP-seq data
#' \itemize{
#'  \item{\code{Uniform_TFBS}: a list (690 combinations of cell types and transcription factors) of GenomicRanges objects; each is an GR object containing uniformly identified peaks per cell type per transcription factor.}
#'  \item{\code{ENCODE_TFBS_ClusteredV3}: a list (161 transcription factors) of GenomicRanges objects; each is an GR object containing clustered peaks per transcription factor, along with a meta-column 'cells' telling cell types associtated with a clustered peak.}
#'  \item{\code{ENCODE_TFBS_ClusteredV3_CellTypes}: a list (91 cell types) of a list (transcription factors) of GenomicRanges objects. Each cell type is a list (transcription factor) of GenomicRanges objects; each is an GR object containing clustered peaks per transcription factor.}
#' }
#' 2. ENCODE DNaseI Hypersensitivity site data
#' \itemize{
#'  \item{\code{Uniform_DNaseI_HS}: a list (125 cell types) of GenomicRanges objects; each is an GR object containing uniformly identified peaks per cell type.}
#'  \item{\code{ENCODE_DNaseI_ClusteredV3}: an GR object containing clustered peaks, along with a meta-column 'num_cells' telling how many cell types associtated with a clustered peak.}
#'  \item{\code{ENCODE_DNaseI_ClusteredV3_CellTypes}: a list (125 cell types) of GenomicRanges objects; each is an GR object containing clustered peaks per cell type.}
#' }
#' 3. ENCODE Histone Modification ChIP-seq data from different sources
#' \itemize{
#'  \item{\code{Broad_Histone}: a list (156 combinations of cell types and histone modifications) of GenomicRanges objects; each is an GR object containing identified peaks per cell type and per histone modification. This dataset was generated from ENCODE/Broad Institute.}
#'  \item{\code{SYDH_Histone}: a list (29 combinations of cell types and histone modifications) of GenomicRanges objects; each is an GR object containing identified peaks per cell type and per histone modification. This dataset was generated from ENCODE/Stanford/Yale/Davis/Harvard.}
#'  \item{\code{UW_Histone}: a list (172 combinations of cell types and histone modifications) of GenomicRanges objects; each is an GR object containing identified peaks per cell type and per histone modification. This dataset was generated from ENCODE/University of Washington.}
#' }
#' 4. FANTOM5 expressed enhancer atlas
#' \itemize{
#'  \item{\code{FANTOM5_Enhancer_Cell}: a list (71 cell types) of GenomicRanges objects; each is an GR object containing enhancers specifically expressed in a cell type.}
#'  \item{\code{FANTOM5_Enhancer_Tissue}: a list (41 tissues) of GenomicRanges objects; each is an GR object containing enhancers specifically expressed in a tissue.}
#'  \item{\code{FANTOM5_Enhancer_Extensive}: a list (5 categories of extensitive enhancers) of GenomicRanges objects; each is an GR object containing extensitive enhancers. They are: "Extensive_ubiquitous_enhancers_cells" for ubiquitous enhancers expressed over the entire set of cell types; "Extensive_ubiquitous_enhancers_organs" for ubiquitous enhancers expressed over the entire set of tissues; "Extensive_enhancers_tss_associations" for TSS-enhancer associations(RefSeq promoters only); "Extensive_permissive_enhancers" and "Extensive_robust_enhancers" for permissive and robust enhancer sets.}
#'  \item{\code{FANTOM5_Enhancer}: a list (117 cell types/tissues/categories) of GenomicRanges objects; each is an GR object.}
#' }
#' 5. ENCODE combined (ChromHMM and Segway) Genome Segmentation data
#' \itemize{
#'  \item{\code{Segment_Combined_Gm12878}: a list (7 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the cell line GM12878 (a lymphoblastoid cell line).}
#'  \item{\code{Segment_Combined_H1hesc}: a list (7 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the cell line H1-hESC (H1 human embryonic stem cells).}
#'  \item{\code{Segment_Combined_Helas3}: a list (7 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the cell line HeLa S3.}
#'  \item{\code{Segment_Combined_Hepg2}: a list (7 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the cell line HepG2 (liver hepatocellular carcinoma).}
#'  \item{\code{Segment_Combined_Huvec}: a list (7 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the cell line HUVEC (Human Umbilical Vein Endothelial Cells).}
#'  \item{\code{Segment_Combined_K562}: a list (7 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the cell line K562 (human erythromyeloblastoid leukemia cell line).}
#' }
#' 6. Conserved TFBS
#' \itemize{
#'  \item{\code{TFBS_Conserved}: a list (245 PWM) of GenomicRanges objects; each is an GR object containing human/mouse/rat conserved TFBS for each PWM.}
#' }
#' 7. TargetScan miRNA regulatory sites
#' \itemize{
#'  \item{\code{TS_miRNA}: a list (153 miRNA) of GenomicRanges objects; each is an GR object containing miRNA regulatory sites for each miRNA.}
#' }
#' 8. TCGA exome mutation data
#' \itemize{
#'  \item{\code{TCGA}: a list (11 tumor types) of GenomicRanges objects; each is an GR object containing exome mutation across tumor patients of the same tumor type.}
#' }
#' 9. ReMap integration of transcription factor ChIP-seq data (publicly available and ENCODE)
#' \itemize{
#'  \item{\code{ReMap_Public_TFBS}: a list (395 combinations of GSE studies and transcription factors and cell types) of GenomicRanges objects; each is an GR object containing identified peaks per GSE study per transcripton factor per cell type.}
#'  \item{\code{ReMap_Public_mergedTFBS}: a list (131 transcription factors under GSE studies) of GenomicRanges objects; each is an GR object containing merged peaks per transcripton factor.}
#'  \item{\code{ReMap_PublicAndEncode_mergedTFBS}: a list (237 transcription factors under GSE studies and ENCODE) of GenomicRanges objects; each is an GR object containing merged peaks per transcripton factor.}
#'  \item{\code{ReMap_Encode_TFBS}: a list (155 transcription factors under ENCODE) of GenomicRanges objects; each is an GR object containing identified peaks per transcripton factor.}
#' }
#' 10. Blueprint Histone Modification ChIP-seq data
#' \itemize{
#'  \item{\code{Blueprint_BoneMarrow_Histone}: a list (132 combinations of histone modifications and samples) of GenomicRanges objects; each is an GR object containing identified peaks per histone per sample (from bone marrow).}
#'  \item{\code{Blueprint_CellLine_Histone}: a list (38 combinations of histone modifications and cell lines) of GenomicRanges objects; each is an GR object containing identified peaks per histone per cell line.}
#'  \item{\code{Blueprint_CordBlood_Histone}: a list (126 combinations of histone modifications and samples) of GenomicRanges objects; each is an GR object containing identified peaks per histone per sample (from cord blood).}
#'  \item{\code{Blueprint_Thymus_Histone}: a list (5 combinations of histone modifications and samples) of GenomicRanges objects; each is an GR object containing identified peaks per histone per sample (from thymus).}
#'  \item{\code{Blueprint_VenousBlood_Histone}: a list (296 combinations of histone modifications and samples) of GenomicRanges objects; each is an GR object containing identified peaks per histone per sample (from venous blood).}
#' }
#' 11. BLUEPRINT DNaseI Hypersensitivity site data
#' \itemize{
#'  \item{\code{Blueprint_DNaseI}: a list (36 samples) of GenomicRanges objects; each is an GR object containing identified peaks per sample.}
#' }
#' 12. BLUEPRINT DNA Methylation data
#' \itemize{
#'  \item{\code{Blueprint_Methylation_hyper}: a list (206 samples) of GenomicRanges objects; each is an GR object containing hyper-methylated CpG regions per sample.}
#'  \item{\code{Blueprint_Methylation_hypo}: a list (206 samples) of GenomicRanges objects; each is an GR object containing hypo-methylated CpG regions per sample.}
#' }
#' 13. Roadmap Epigenomics Core 15-state Genome Segmentation data for primary cells (blood and T cells)
#' \itemize{
#' \item{\code{EpigenomeAtlas_15Segments_E033}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E033 (Primary T cells from cord blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E034}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E034 (Primary T cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E037}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E037 (Primary T helper memory cells from peripheral blood 2).}
#' \item{\code{EpigenomeAtlas_15Segments_E038}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E038 (Primary T helper naive cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E039}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E039 (Primary T helper naive cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E040}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E040 (Primary T helper memory cells from peripheral blood 1).}
#' \item{\code{EpigenomeAtlas_15Segments_E041}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E041 (Primary T helper cells PMA-I stimulated).}
#' \item{\code{EpigenomeAtlas_15Segments_E042}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E042 (Primary T helper 17 cells PMA-I stimulated).}
#' \item{\code{EpigenomeAtlas_15Segments_E043}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E043 (Primary T helper cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E044}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E044 (Primary T regulatory cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E045}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E045 (Primary T cells effector/memory enriched from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E047}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E047 (Primary T killer naive cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E048}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E048 (Primary T killer memory cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E062}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E062 (Primary mononuclear cells from peripheral blood).}
#' }
#' 14. Roadmap Epigenomics Core 15-state Genome Segmentation data for primary cells (HSC and B cells)
#' \itemize{
#' \item{\code{EpigenomeAtlas_15Segments_E029}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E029 (Primary monocytes from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E030}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E030 (Primary neutrophils from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E031}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E031 (Primary B cells from cord blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E032}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E032 (Primary B cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E035}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E035 (Primary hematopoietic stem cells).}
#' \item{\code{EpigenomeAtlas_15Segments_E036}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E036 (Primary hematopoietic stem cells short term culture).}
#' \item{\code{EpigenomeAtlas_15Segments_E046}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E046 (Primary Natural Killer cells from peripheral blood).}
#' \item{\code{EpigenomeAtlas_15Segments_E050}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E050 (Primary hematopoietic stem cells G-CSF-mobilized Female).}
#' \item{\code{EpigenomeAtlas_15Segments_E051}: a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome E051 (Primary hematopoietic stem cells G-CSF-mobilized Male).}
#' }
#' 15. CpG annotation
#' \itemize{
#'  \item{\code{CpG_anno}: a list (4 categories) of GenomicRanges objects; each is an GR object. They are exclusive, including (in order) "CpG_islands", "CpG_shores" (2Kb upstream/downstream from the ends of the CpG islands), "CpG_shelves" (2Kb upstream/downstream of the farthest upstream/downstream limits of the CpG shores), and "CpG_inter" (the remaining inter-CGI genomic regions 'open sea'). }
#' }
#' 16. Genic annotation
#' \itemize{
#'  \item{\code{Genic_anno}: a list (12 categories) of GenomicRanges objects; each is an GR object. They are not exclusively, including "Genic_1to5kb" (1-5Kb upstream of TSS), "Genic_promoters" (1Kb upstream of TSS), "Genic_5UTRs", "Genic_firstexons" (first exons), "Genic_exons", "Genic_exonintronboundaries", "Genic_introns", "Genic_intronexonboundaries", "Genic_cds", "Genic_3UTRs", "Genic_intergenic" (the intergenic regions exclude the previous list of annotations), and "Genic_lncRNA" (GENCODE long non-coding RNA (lncRNA) transcripts). }
#' }
#' 17. FANTOM5 sample-ontology-enriched CAT genes
#' \itemize{
#'  \item{\code{FANTOM5_CAT_Cell}: a list (173 cell types) of GenomicRanges objects; each is an GR object containing CAT genes specifically expressed in a cell type.}
#'  \item{\code{FANTOM5_CAT_Tissue}: a list (174 tissues) of GenomicRanges objects; each is an GR object containing CAT genes specifically expressed in a tissue.}
#' }
#' 18. FANTOM5 trait-associated CAT genes
#' \itemize{
#'  \item{\code{FANTOM5_CAT_DO}: a list (299 traits grouped by disease ontology) of GenomicRanges objects; each is an GR object containing CAT genes harboring at least one trait-associated SNP.}
#'  \item{\code{FANTOM5_CAT_EFO}: a list (93 traits grouped by experiment factor ontology) of GenomicRanges objects; each is an GR object containing CAT genes harboring at least one trait-associated SNP.}
#'  \item{\code{FANTOM5_CAT_HPO}: a list (176 traits grouped by human phenotype ontology) of GenomicRanges objects; each is an GR object containing CAT genes harboring at least one trait-associated SNP.}
#'  \item{\code{FANTOM5_CAT_MESH}: a list (210 traits grouped by Medical Subject Headings) of GenomicRanges objects; each is an GR object containing CAT genes harboring at least one trait-associated SNP.}
#'  \item{\code{FANTOM5_CAT_PICS}: a list (39 traits grouped by PICS dieases) of GenomicRanges objects; each is an GR object containing CAT genes harboring at least one trait-associated SNP.}
#' }
#' @export
#' @seealso \code{\link{xEnrichViewer}}
#' @include xLDenricher.r
#' @examples
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#'
#' \dontrun{
#' # a) provide the seed SNPs with the significance info
#' ## load ImmunoBase
#' data(ImmunoBase)
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' data <- GenomicRanges::mcols(gr)[,c(1,3)]
#'
#' # b) get LD block (EUR population)
#' bLD <- xLDblock(data, include.LD="EUR", LD.r2=0.8, RData.location=RData.location)
#' 
#' ## c) perform enrichment analysis using FANTOM expressed enhancers
#' eTerm <- xLDenricher(bLD, GR.annotation="FANTOM5_Enhancer_Cell", num.samples=20000, preserve=c("boundary","exact"), RData.location=RData.location)
#'
#' ## d) view enrichment results for the top significant terms
#' xEnrichViewer(eTerm)
#'
#' ## e) barplot of enriched terms
#' bp <- xEnrichBarplot(eTerm, top_num='auto', displayBy="fdr")
#' bp
#'
#' ## f) save enrichment results to the file called 'LD_enrichments.txt'
#' output <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
#' utils::write.table(output, file="LD_enrichments.txt", sep="\t", row.names=FALSE)
#' 
#' ## g) compare boundary and exact
#' GR.SNP <- xRDataLoader("dbSNP_GWAS", RData.location=RData.location)
#' GR.annotation <- xRDataLoader("FANTOM5_Enhancer_Cell", RData.location=RData.location)
#' eTerm_boundary <- xLDenricher(bLD, GR.SNP=GR.SNP, GR.annotation=GR.annotation, num.samples=20000, preserve="boundary", RData.location=RData.location)
#' eTerm_exact <- xLDenricher(bLD, GR.SNP=GR.SNP, GR.annotation=GR.annotation, num.samples=20000, preserve="exact", RData.location=RData.location)
#' list_eTerm <- list(boundary=eTerm_boundary, exact=eTerm_exact)
#' bp <- xEnrichCompare(list_eTerm, displayBy="zscore")
#' }

xLDenricher <- function(bLD, GR.SNP=c("dbSNP_GWAS","dbSNP_Common","dbSNP_Single"), num.samples=2000, respect=c("maf","distance","both"), restrict.chr=F, preserve=c("boundary","exact"), seed=825, p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), GR.annotation=c(NA,"Uniform_TFBS","ENCODE_TFBS_ClusteredV3","ENCODE_TFBS_ClusteredV3_CellTypes", "Uniform_DNaseI_HS","ENCODE_DNaseI_ClusteredV3","ENCODE_DNaseI_ClusteredV3_CellTypes", "Broad_Histone","SYDH_Histone","UW_Histone","FANTOM5_Enhancer_Cell","FANTOM5_Enhancer_Tissue","FANTOM5_Enhancer_Extensive","FANTOM5_Enhancer","Segment_Combined_Gm12878","Segment_Combined_H1hesc","Segment_Combined_Helas3","Segment_Combined_Hepg2","Segment_Combined_Huvec","Segment_Combined_K562","TFBS_Conserved","TS_miRNA","TCGA", "ReMap_Public_TFBS","ReMap_Public_mergedTFBS","ReMap_PublicAndEncode_mergedTFBS","ReMap_Encode_TFBS", "Blueprint_BoneMarrow_Histone","Blueprint_CellLine_Histone","Blueprint_CordBlood_Histone","Blueprint_Thymus_Histone","Blueprint_VenousBlood_Histone","Blueprint_DNaseI","Blueprint_Methylation_hyper","Blueprint_Methylation_hypo","EpigenomeAtlas_15Segments_E029", "EpigenomeAtlas_15Segments_E030", "EpigenomeAtlas_15Segments_E031", "EpigenomeAtlas_15Segments_E032", "EpigenomeAtlas_15Segments_E033", "EpigenomeAtlas_15Segments_E034", "EpigenomeAtlas_15Segments_E035", "EpigenomeAtlas_15Segments_E036", "EpigenomeAtlas_15Segments_E037", "EpigenomeAtlas_15Segments_E038", "EpigenomeAtlas_15Segments_E039", "EpigenomeAtlas_15Segments_E040", "EpigenomeAtlas_15Segments_E041", "EpigenomeAtlas_15Segments_E042", "EpigenomeAtlas_15Segments_E043", "EpigenomeAtlas_15Segments_E044", "EpigenomeAtlas_15Segments_E045", "EpigenomeAtlas_15Segments_E046", "EpigenomeAtlas_15Segments_E047", "EpigenomeAtlas_15Segments_E048", "EpigenomeAtlas_15Segments_E050", "EpigenomeAtlas_15Segments_E051", "EpigenomeAtlas_15Segments_E062", "CpG_anno","Genic_anno", "FANTOM5_CAT_Cell","FANTOM5_CAT_Tissue","FANTOM5_CAT_DO","FANTOM5_CAT_EFO","FANTOM5_CAT_HPO","FANTOM5_CAT_MESH","FANTOM5_CAT_PICS"), verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    respect <- match.arg(respect)
    preserve <- match.arg(preserve)
    p.adjust.method <- match.arg(p.adjust.method)
    
	#######################################################
	if(verbose){
		message(sprintf("First, generate null distribution via doing %d sampling (%s) ...", num.samples, as.character(Sys.time())), appendLF=T)
	}
	
    if(verbose){
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("'xLDsampling' is being called (%s):", as.character(Sys.time())), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
	grl <- xLDsampling(bLD=bLD, GR.SNP=GR.SNP, num.samples=num.samples, respect=respect, restrict.chr=restrict.chr, preserve=preserve, seed=seed, verbose=verbose, RData.location=RData.location)
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xLDsampling' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
	
	#######################################################
	
	if(verbose){
		message(sprintf("Second, load GR annotations (%s) ...", as.character(now <- Sys.time())), appendLF=T)
	}
	
	if(class(GR.annotation) == "list"){
		aGR <- GR.annotation
	}else{
		if(is.na(GR.annotation)){
			stop("Please specify annotation RData!\n")
		}else{
			if(length(GR.annotation)>1){
				message("\tONLY the first specified annotation RData will be used!\n")
				GR.annotation <- GR.annotation[1]
			}
			aGR <- xRDataLoader(RData.customised=GR.annotation, verbose=verbose, RData.location=RData.location)
			if(is.null(aGR)){
				stop("Your specified annotation RData does not exist!\n")
			}
		}
	}
	aGRL <- GenomicRanges::GRangesList(aGR)	
	
	#######################################################
	
	if(verbose){
		message(sprintf("Third, perform enrichment analysis (%s) ...", as.character(now <- Sys.time())), appendLF=T)
	}
	
	queryHits <- subjectHits <- B <- n <- best <- NULL
	
	## observed
	### gr_observed
	if(preserve=="boundary"){
		gr_best <- bLD$best
		df_best <- GenomicRanges::as.data.frame(gr_best)
		df_best$start <- df_best$start + df_best$upstream
		df_best$end <- df_best$end + df_best$downstream
		gr_observed <- xGR(df_best, format='data.frame', add.name=F)
		gr_observed$best <- rownames(df_best)
		
	}else if(preserve=="exact"){
		grl_block <- bLD$block
		gr_block <- BiocGenerics::unlist(grl_block,use.names=F)
		gr_observed <- gr_block[,'best']
		
	}
	
	### obs
	system.time({
	q2r <- as.data.frame(GenomicRanges::findOverlaps(gr_observed, aGRL))
	q2r$best <- gr_observed$best[q2r$queryHits]
	df_n_observed <- as.data.frame(q2r %>% dplyr::group_by(subjectHits) %>% dplyr::summarise(n=length(unique(best))))
	obs <- rep(0, length(aGRL))
	obs[df_n_observed$subjectHits] <- df_n_observed$n
	names(obs) <- names(aGRL)
	})
	
	if(verbose){
		message(sprintf("\t%d observed (%s)", length(obs), as.character(Sys.time())), appendLF=T)
	}
	
	## expected
	### gr_expected
	gr_expected <- BiocGenerics::unlist(grl)
	### a2B
	system.time({
	q2r <- as.data.frame(GenomicRanges::findOverlaps(gr_expected, aGRL))
	q2r$best <- gr_expected$best[q2r$queryHits]
	q2r$B <- gr_expected$B[q2r$queryHits]
	})
	system.time({
	df_n_expected <- as.data.frame(q2r %>% dplyr::group_by(subjectHits,B) %>% dplyr::summarise(n=length(unique(best))))
	a2B <- as.matrix(xSparseMatrix(df_n_expected, rows=1:length(aGRL), columns=1:length(grl), verbose=T))
	rownames(a2B) <- names(aGRL)
	})

	if(verbose){
		message(sprintf("\t%d x %d of the expected matrix (%s)", nrow(a2B), ncol(a2B), as.character(now <- Sys.time())), appendLF=T)
	}

	### exp
	exp_mean <- apply(a2B, 1, mean)
	exp_std <- apply(a2B, 1, stats::sd)
	
	## ratio
	ratio <- obs/exp_mean
	
    ## zscore
    zscore <- (obs-exp_mean)/exp_std

    ## pvalue
    obs_matrix <- matrix(rep(obs,ncol(a2B)), ncol=ncol(a2B))
    pvalue <- apply((obs_matrix - a2B)<=0, 1, sum) / ncol(a2B)
	####################
	
    zscore[is.na(zscore)] <- 0
    zscore[is.infinite(zscore)] <- max(zscore[!is.infinite(zscore)])
    pvalue[is.na(ratio)] <- 1
    ratio[is.na(ratio)] <- 1
 	
 	####################
 	or <- CIl <- CIu <- NA
 	####################
 
	enrichment_df <- data.frame(names(aGRL), sapply(aGR,length), length(gr_observed), obs, exp_mean, exp_std, ratio, zscore, pvalue, or, CIl, CIu, row.names=NULL, stringsAsFactors=F)
	colnames(enrichment_df) <- c("name", "nAnno", "nData", "nOverlap", "nExpect", "std", "fc", "zscore", "pvalue", "or", "CIl", "CIu")

	## Adjust P-values for multiple comparisons
	pvals <- enrichment_df$pvalue
	adjpvals <- stats::p.adjust(pvals, method=p.adjust.method)
	enrichment_df$adjp <- adjpvals

	####################################################################################
	
	enrichment_df$zscore <- signif(enrichment_df$zscore, digits=3)
	
	pvals <- enrichment_df$pvalue
	adjpvals <- enrichment_df$adjp
	pvals <- signif(pvals, digits=2)
	adjpvals <- signif(adjpvals, digits=2)
	
	# scientific notations
	pvals  <- base::sapply(pvals, function(x){
		if(x < 0.1 & x!=0){
			as.numeric(format(x,scientific=T))
		}else{
			x
		}
	})
	
	adjpvals <- base::sapply(adjpvals, function(x){
		if(x < 0.1 & x!=0){
			as.numeric(format(x,scientific=T))
		}else{
			x
		}
	})
	
	enrichment_df$pvalue <- pvals
	enrichment_df$adjp <- adjpvals
	
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
	res_df <- enrichment_df[, c("name", "nAnno", "nOverlap", "fc", "zscore", "pvalue", "adjp","or", "CIl", "CIu", "nData", "nExpect","std")]
	
	invisible(res_df)
}
