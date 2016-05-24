#' Function to visualise a network as a circos plot
#'
#' \code{xCircos} is used to visualise a network as a circos plot. The network must be a 'igraph' object.
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
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 146) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 146) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "RData.location"
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_genes', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Even the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location"
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
#' # Load the library
#' library(XGR)
#' library(RCircos)
#' RData.location="~/Sites/SVN/github/RDataCentre/Portal"
#' 
#' # provide genes and SNPs reported in AS GWAS studies
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase')
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
#' } 

xCircos <- function(g, entity=c("SNP","Gene"), top_num=50, colormap=c("yr","bwr","jet","gbr","wyr","br","rainbow","wb","lightyellow-orange"), rescale=T, nodes.query=NULL, ideogram=T, chr.exclude="auto", entity.label.cex=0.7, GR.SNP="dbSNP_GWAS", GR.Gene="UCSC_genes", verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/Portal")
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
	entity <- match.arg(entity)
	
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

  	if(entity=="SNP") {
    	pos <- xRDataLoader(RData.customised=GR.SNP, verbose=verbose, RData.location=RData.location)
    	if(is.null(pos)){
    		GR.SNP <- "dbSNP_GWAS"
			if(verbose){
				message(sprintf("Instead, %s will be used", GR.SNP), appendLF=T)
			}
    		pos <- xRDataLoader(RData.customised=GR.SNP, verbose=verbose, RData.location=RData.location)
    	}
  	}else if(entity == "Gene") {
    	pos <- xRDataLoader(RData.customised=GR.Gene, verbose=verbose, RData.location=RData.location)
    	if(is.null(pos)){
    		GR.Gene <- "UCSC_genes"
			if(verbose){
				message(sprintf("Instead, %s will be used", GR.Gene), appendLF=T)
			}
    		pos <- xRDataLoader(RData.customised=GR.Gene, verbose=verbose, RData.location=RData.location)
    	}
  	}else{
    	stop("Please indicate whether your analysis entity was SNPs or genes.\n")
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
  	RCircos.Set.Core.Components(cyto.info, chr.exclude, num.inside, num.outside)
  
  	## Reset parameters
  	params <- RCircos.Get.Plot.Parameters()
  	params$track.padding <- 0 # 0.02
  	params$track.height <- 0.05 # 0.1
  	
  	params$chr.ideog.pos <- 1
  	params$highlight.pos <- 1.1
  	params$chr.name.pos <- 1.1
  	params$plot.radius <- 0.9
  	params$track.out.start <- 1.2
  	params$highlight.width <- 0
  	
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
	
	palette.name <- supraHex::visColormap(colormap=colormap)
  	input.data$PlotColor <- palette.name(20)[as.numeric(cut(input.data$similarity, breaks=seq(0, 1, 0.05)))]
  	input.data <- input.data[order(input.data$similarity, decreasing=F), ]
  	RCircos.Link.Plot(input.data, track.num=1, FALSE)

  	## Label SNPs/genes in outside track
	if(verbose){
		now <- Sys.time()
		message(sprintf("Adding SNP or gene names (%s) ...", as.character(now)), appendLF=T)
	}
  	name.col <- "Name"
  	side <- "out"
  	track.num <- 1
  	label.data <- label.data[!duplicated(label.data$Name), ]
  	if(verbose){
  		RCircos.Gene.Name.Plot(label.data, name.col, track.num, side)
  	}else{
  		suppressMessages(RCircos.Gene.Name.Plot(label.data, name.col, track.num, side))
  	}
  	
  	invisible()
}
