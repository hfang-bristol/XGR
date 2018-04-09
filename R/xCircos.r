#' Function to visualise a network as a circos plot
#'
#' \code{xCircos} is used to visualise a network as a circos plot. The network must be a 'igraph' object. The degree of similarity between SNPs (or genes) is visualised by the colour of links. This function can be used either to visualise the most similar links or to plot links involving an input SNP (or gene).
#'
#' @param g an object of class "igraph". For example, it stores semantic similarity results with nodes for genes/SNPs and edges for pair-wise semantic similarity between them 
#' @param entity the entity of similarity analysis for which results are being plotted. It can be either "SNP" or "Gene"
#' @param top_num the top number of similarity edges to be plotted
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param rescale logical to indicate whether the edge values are rescaled to the range [0,1]. By default, it sets to true
#' @param nodes.query nodes in query for which edges attached to them will be displayed. By default, it sets to NULL meaning no such restriction
#' @param ideogram logical to indicate whether chromosome banding is plotted
#' @param chr.exclude a character vector of chromosomes to exclude from the plot, e.g. c("chrX", "chrY"). By defautl, it is 'auto' meaning those chromosomes without data will be excluded. If NULL, no chromosome is excluded
#' @param entity.label.cex the font size of genes/SNPs labels. Default is 0.8
#' @param entity.label.side the position of genes/SNPs labels relative to chromosome ideogram. It can be "out" (by default) or "in"
#' @param entity.label.track an integer specifying the plot track for genes/SNPs labels. Default is 1
#' @param entity.label.query which genes/SNPs in query will be labelled. By default, it sets to NULL meaning all will be displayed. If labels in query can not be found, then all will be displayed
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 146) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 146) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' a circos plot with edge weights between input snps/genes represented by the colour of the links
#' @note none
#' @export
#' @import RCircos
#' @seealso \code{\link{xSocialiserGenes}}, \code{\link{xSocialiserSNPs}}
#' @include xCircos.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
#' 
#' library(RCircos)
#'
#' # provide genes and SNPs reported in GWAS studies
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' 
#' # 1) SNP-based similarity analysis using GWAS Catalog traits (mapped to EF)
#' ## Get lead SNPs reported in AS GWAS
#' example.snps <- names(ImmunoBase$AS$variants)
#' SNP.g <- xSocialiserSNPs(example.snps, include.LD=NA, RData.location=RData.location)
#' # Circos plot of the EF-based SNP similarity network
#' #out.file <- "SNP_Circos.pdf"
#' #pdf(file=out.file, height=12, width=12, compress=TRUE)
#' xCircos(g=SNP.g, entity="SNP", RData.location=RData.location)
#' #dev.off()
#' # Circos plot involving nodes 'rs6871626'
#' xCircos(g=SNP.g, entity="SNP", nodes.query="rs6871626", RData.location=RData.location)
#'
#' # 2) Gene-based similarity analysis using Disease Ontology (DO)
#' ## Get genes within 10kb away from AS GWAS lead SNPs
#' example.genes <- names(which(ImmunoBase$AS$genes_variants<=10000))
#' gene.g <- xSocialiserGenes(example.genes, ontology="DO", RData.location=RData.location)
#' # Circos plot of the DO-based gene similarity network
#' #out.file <- "Gene_Circos.pdf"
#' #pdf(file=out.file, height=12, width=12, compress=TRUE)
#' xCircos(g=gene.g, entity="Gene", chr.exclude="chrY", RData.location=RData.location)
#' #dev.off()
#'
#' # 3) Advanced usages: Gene-SNP pairs from trans-eQTL mapping
#' JKscience_TS3A <- xRDataLoader(RData.customised='JKscience_TS3A', RData.location=RData.location)
#' ## extract the significant trans-eQTL in IFN
#' ind <- -1*log10(JKscience_TS3A$IFN_fdr)
#' ind <- which(!is.na(ind) & ind>2)
#' relations <- JKscience_TS3A[ind, c("Symbol","variant","IFN_fdr")]
#' relations <- data.frame(from=relations$Symbol, to=relations$variant, weight=-log10(relations$IFN_fdr))
#' ig_Gene2SNP <- igraph::graph.data.frame(d=relations, directed=TRUE)
#' # Circos plot of the eQTL (Gene-SNP) network
#' #out.file <- "eQTL_Circos.pdf"
#' #pdf(file=out.file, height=12, width=12, compress=TRUE)
#' xCircos(g=ig_Gene2SNP, entity="Both", top_num=NULL, nodes.query=c("GAD1","TNFRSF1B"), chr.exclude=NULL, entity.label.side="out", RData.location=RData.location)
#' #dev.off()
#' }

xCircos <- function(g, entity=c("SNP","Gene","Both"), top_num=50, colormap=c("yr","bwr","jet","gbr","wyr","br","rainbow","wb","lightyellow-orange"), rescale=T, nodes.query=NULL, ideogram=T, chr.exclude="auto", entity.label.cex=0.7, entity.label.side=c("in","out"), entity.label.track=1, entity.label.query=NULL, GR.SNP=c("dbSNP_GWAS","dbSNP_Common"), GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
	entity <- match.arg(entity)
	entity.label.side <- match.arg(entity.label.side)
	
	flag_package <- F
    pkgs <- c("RCircos")
    if(all(pkgs %in% rownames(utils::installed.packages()))){
        tmp <- sapply(pkgs, function(pkg) {
            #suppressPackageStartupMessages(require(pkg, character.only=T))
            requireNamespace(pkg, quietly=T)
        })
        if(all(tmp)){
        	flag_package <- T
        }
    }
	if(!flag_package){
		stop("The package 'RCircos' is not available.\n")
	}
	
  	## Check input g
  	if (class(g) != "igraph") {
    	stop("The function must apply to a 'igraph' object.\n")
  	}

	## Convert from igraph into data.frame
  	df <- igraph::get.data.frame(g, what="edges")
	
	## check the weight and sort the weight
	if(is.null(df$weight)){
		df$weight <- rep(1, nrow(df))
		## force NOT to rescale weight
		rescale <- F
	}else{
		df$weight <- as.numeric(df$weight)
	}
	df <- df[with(df,order(-weight)), ]
	
	## restrict to nodes in query
	if(!is.null(nodes.query)){
		ind_from <- which(!is.na(match(df[,1], nodes.query)))
		ind_to <- which(!is.na(match(df[,2], nodes.query)))
		ind <- union(ind_from, ind_to)
		if(length(ind)>0){
			df <- df[ind,]
			
			if(verbose){
				ind <- match(nodes.query, union(df[,1], df[,2]))
				nodes.query <- nodes.query[!is.na(ind)]
				now <- Sys.time()
				message(sprintf("Circos plot restricted to nodes '%s' (%s) ...", paste(nodes.query,collapse=','), as.character(now)), appendLF=T)
			}
		}
	}
	
	## keep the top edges
  	if(is.null(top_num)){
    	top_num <- nrow(df)
  	}
  	if(top_num > nrow(df)){
    	top_num <- nrow(df)
  	}
  	top_num <- as.integer(top_num)
  	df <- df[1:top_num, ]
  
  	## load positional information
	if(verbose){
		now <- Sys.time()
		message(sprintf("Loading positional information for %s (%s) ...", entity, as.character(now)), appendLF=T)
	}

  	if(entity=="SNP" | entity=="Both"){
		if(class(GR.SNP) == "GRanges"){
			pos_snp <- GR.SNP
		}else{
			pos_snp <- xRDataLoader(RData.customised=GR.SNP[1], verbose=verbose, RData.location=RData.location)
			if(is.null(pos_snp)){
				GR.SNP <- "dbSNP_GWAS"
				if(verbose){
					message(sprintf("Instead, %s will be used", GR.SNP), appendLF=T)
				}
				pos_snp <- xRDataLoader(RData.customised=GR.SNP, verbose=verbose, RData.location=RData.location)
			}		
		}
  	}
  	
  	if(entity == "Gene" | entity=="Both"){
		if(class(GR.Gene) == "GRanges"){
			pos_gene <- GR.Gene
		}else{
			pos_gene <- xRDataLoader(RData.customised=GR.Gene[1], verbose=verbose, RData.location=RData.location)
			if(is.null(pos_gene)){
				GR.Gene <- "UCSC_knownGene"
				if(verbose){
					message(sprintf("Instead, %s will be used", GR.Gene), appendLF=T)
				}
				pos_gene <- xRDataLoader(RData.customised=GR.Gene, verbose=verbose, RData.location=RData.location)
			}
		}
  	}
  	
  	if(entity=="SNP"){
  		pos <- pos_snp
  	}else if(entity=="Gene"){
  		pos <- pos_gene
  	}else if(entity=="Both"){
    	## Combined both
    	GenomicRanges::mcols(pos_gene) <- NULL
    	GenomicRanges::mcols(pos_snp) <- NULL
    	pos <- c(pos_gene, pos_snp)
  	}

  	## Convert into format required for Circos plot
  	allnames <- names(pos)

    A <- match(df$from, allnames)
	B <- match(df$to, allnames)
	#flag <- complete.cases(cbind(A, B))
	flag <- !is.na(A) & !is.na(B)
	AA <- A[flag]
	BB <- B[flag]
	input.data.A <- GenomicRanges::as.data.frame(pos[AA], row.names=NULL)
	input.data.B <- GenomicRanges::as.data.frame(pos[BB], row.names=NULL)
	input.data <- cbind.data.frame(input.data.A[, 1:3], input.data.B[, 1:3])
	if(is.null(df$weight)){
		input.data$similarity <- rep(1, sum(flag))
	}else{
		input.data$similarity <- as.numeric(as.character(df$weight[flag]))
	}
	label.data <- rbind(input.data.A[, 1:3], input.data.B[, 1:3])
	label.data$Name <- c(df$from[flag], df$to[flag])
  	
  	## decide on which chromosomes will be excluded
  	if(!is.null(chr.exclude)){
  		chr.exclude <- chr.exclude[!is.na(chr.exclude)]
		if(length(chr.exclude)==0){
			chr.exclude <- NULL
		}else if(sum(chr.exclude=='auto')>0){
			flag <- levels(label.data$seqnames) %in% as.character(unique(label.data$seqnames))
			chr.exclude <- levels(label.data$seqnames)[!flag]
		}
  	}
  	
  	## Load human chromosome ideogram
	if(verbose){
		now <- Sys.time()
		message(sprintf("Loading human chromosome banding information (hg19) (%s) ...", as.character(now)), appendLF=T)
	}
	
	#data(UCSC.HG19.Human.CytoBandIdeogram, package="RCircos")
	eval(parse(text="data(UCSC.HG19.Human.CytoBandIdeogram)"))
	
  	cyto.info <- ""
  	eval(parse(text=paste("cyto.info <- UCSC.HG19.Human.CytoBandIdeogram", sep="")))
  	if(ideogram==F) {
    	cyto.info$Stain <- rep("gpos100", nrow(cyto.info))
  	}
  	
  	## Set RCircos core components
	if(verbose){
		now <- Sys.time()
		message(sprintf("Initialising RCircos Core Components (%s) ...", as.character(now)), appendLF=T)
	}
  	num.inside <- 1
  	num.outside <- 1
  	suppressMessages(RCircos.Set.Core.Components(cyto.info, chr.exclude, num.inside, num.outside))

  	## Reset parameters
  	params <- RCircos.Get.Plot.Parameters()  
  	if(0){
  	params$track.padding <- 0 # 0.02
  	params$track.height <- 0.05 # 0.1
  	
  	params$chr.ideog.pos <- 1
  	params$highlight.pos <- 1.1
  	params$chr.name.pos <- 1.1
  	#params$plot.radius <- 0.9
  	params$track.out.start <- 1.2
  	params$highlight.width <- 0 
	}
  	params$text.size <- entity.label.cex
  	RCircos.Reset.Plot.Parameters(params)

  	## Initialise graphic device, plot chromosome ideogram
	if(verbose){
		now <- Sys.time()
		message(sprintf("Plotting chromosome ideogram (%s) ...", as.character(now)), appendLF=T)
	}
  	RCircos.Set.Plot.Area()
  	RCircos.Chromosome.Ideogram.Plot()

  	## Plot link data coloured according to the similarity output
	if(verbose){
		now <- Sys.time()
		message(sprintf("Plotting link data (%s) ...", as.character(now)), appendLF=T)
	}
	
	## Also rescale similarity into the [0,1] range
	if(rescale){
		sim <- input.data$similarity
		if(verbose){
			now <- Sys.time()
			message(sprintf("Also rescale similarity into the [0,1] range (%s)", as.character(now)), appendLF=T)
		}
		# rescale to [0 1]
		input.data$similarity <- (sim - min(sim))/(max(sim) - min(sim))
	}
	
	palette.name <- xColormap(colormap=colormap)
	#cut_index <- as.numeric(cut(input.data$similarity, breaks=seq(0, 1, 0.05)))
	#cut_index[is.na(cut_index)] <- 1
	cut_index <- as.numeric(cut(input.data$similarity, breaks=min(input.data$similarity)+max(input.data$similarity)*seq(0.05, 1, 0.05)))
	#cut_index[is.na(cut_index)] <- max(input.data$similarity)
  	input.data$PlotColor <- palette.name(20)[cut_index]
  	input.data <- input.data[order(input.data$similarity, decreasing=F), ]
  	RCircos.Link.Plot(input.data, track.num=1, FALSE)

  	## Label SNPs/genes in outside track
	if(verbose){
		now <- Sys.time()
		message(sprintf("Adding SNP and/or gene names (%s) ...", as.character(now)), appendLF=T)
	}
  	name.col <- "Name"
  	side <- entity.label.side
  	track.num <- entity.label.track
  	label.data <- label.data[!duplicated(label.data$Name), ]
  	if(!is.null(entity.label.query)){
  		ind <- match(label.data$Name, entity.label.query)
  		if(sum(!is.na(ind)) >= 1){
  			label.data <- label.data[!is.na(ind), ]
  		}
  	}
  	if(verbose){
  		RCircos.Gene.Name.Plot(label.data, name.col, track.num, side)
  	}else{
  		suppressMessages(RCircos.Gene.Name.Plot(label.data, name.col, track.num, side))
  	}
  	
  	invisible()
}
