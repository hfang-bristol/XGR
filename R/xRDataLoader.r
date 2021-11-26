#' Function to load the package built-in RData
#'
#' \code{xRDataLoader} is supposed to load the package built-in RData.
#'
#' @param RData which built-in RData to load. It can be one of "GWAS2EF", "GWAS_LD", "IlluminaHumanHT", "IlluminaOmniExpress", "ig.DO", "ig.EF", "ig.GOBP", "ig.GOCC", "ig.GOMF", "ig.HPCM", "ig.HPMA", "ig.HPMI", "ig.HPPA", "ig.MP", "org.Hs.eg", "org.Hs.egDGIdb", "org.Hs.egDO", "org.Hs.egGOBP", "org.Hs.egGOCC", "org.Hs.egGOMF", "org.Hs.egHPCM", "org.Hs.egHPMA", "org.Hs.egHPMI", "org.Hs.egHPPA", "org.Hs.egMP", "org.Hs.egMsigdbC1", "org.Hs.egMsigdbC2BIOCARTA", "org.Hs.egMsigdbC2CGP", "org.Hs.egMsigdbC2CPall", "org.Hs.egMsigdbC2CP", "org.Hs.egMsigdbC2KEGG", "org.Hs.egMsigdbC2REACTOME", "org.Hs.egMsigdbC3MIR", "org.Hs.egMsigdbC3TFT", "org.Hs.egMsigdbC4CGN", "org.Hs.egMsigdbC4CM", "org.Hs.egMsigdbC5BP", "org.Hs.egMsigdbC5CC", "org.Hs.egMsigdbC5MF", "org.Hs.egMsigdbC6", "org.Hs.egMsigdbC7", "org.Hs.egMsigdbH", "org.Hs.egPS", "org.Hs.egSF", "org.Hs.egPfam", "org.Hs.string", "org.Hs.PCommons_DN", "org.Hs.PCommons_UN", "org.Hs.egGTExV4", "org.Hs.egGTExV6"
#' @param RData.customised a file name for RData-formatted file. By default, it is NULL. It is designed when the user wants to import customised RData that are not listed in the above argument 'RData'. However, this argument can be always used even for those RData that are listed in the argument 'RData' 
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param RData.location the characters to tell the location of built-in RData files. By default, it remotely locates at \url{http://galahad.well.ox.ac.uk/bigdata}; the development version locates at \url{http://galahad.well.ox.ac.uk/bigdata}. For the user equipped with fast internet connection, this option can be just left as default. But it is always advisable to download these files locally. Especially when the user needs to run this function many times, there is no need to ask the function to remotely download every time (also it will unnecessarily increase the runtime). For examples, these files (as a whole or part of them) can be first downloaded into your current working directory, and then set this option as: \eqn{RData.location="."}. Surely, the location can be anywhere as long as the user provides the correct path pointing to (otherwise, the script will have to remotely download each time)
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. For example, 'gskpn' (see 'https://osf.io/gskpn'). If a valid provided and the query matched, it has priority over the one specified via RData.location
#' @return 
#' any use-specified variable that is given on the right side of the assigement sign '<-', which contains the loaded RData. If the data cannot be loaded, it returns NULL.
#' @note If there are no use-specified variable that is given on the right side of the assigement sign '<-', then no RData will be loaded onto the working environment. To enable 'guid', please also install a package "osfr" via \code{BiocManager::install("osfr",dependencies=TRUE)}.
#' @export
#' @import dnet
#' @import igraph
#' @import ggplot2
#' @import RCircos
#' @importFrom GenomicRanges findOverlaps distance mcols seqnames as.data.frame GRangesList GRanges split start end
#' @importFrom IRanges IRanges width pintersect reduce
#' @importFrom grDevices colorRampPalette dev.cur rgb dev.new rainbow hcl extendrange dev.off pdf col2rgb jpeg
#' @importFrom graphics plot lines legend contour text par hist curve abline
#' @importFrom supraHex visColormap visTreeBootstrap visHeatmapAdv
#' @importFrom stats sd median mad ecdf na.omit predict prcomp lm quantile as.dist hclust cor as.dendrogram order.dendrogram wilcox.test coef p.adjust dist
#' @importFrom BiocGenerics unlist start end
#' @importFrom tibble tibble enframe as_tibble
#' @importFrom dplyr select filter arrange mutate group_by summarise desc n arrange_all slice left_join pull bind_rows semi_join transmute distinct n_distinct
#' @importFrom purrr map_chr
#' @importFrom readr write_delim read_delim
#' @importFrom ggnetwork ggnetwork geom_nodes geom_edges
#' @importFrom ggrepel geom_text_repel geom_label_repel GeomTextRepel
#' @importFrom Matrix Diagonal colSums Matrix t
#' @importFrom MASS fitdistr
#' @importFrom osfr osf_retrieve_node osf_ls_files osf_download
#' @importFrom methods is
#' @seealso \code{\link{xRDataLoader}}
#' @include xRDataLoader.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase')
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg')
#' ig.HPPA <- xRDataLoader(RData='ig.HPPA')
#' org.Hs.egHPPA <- xRDataLoader(RData='org.Hs.egHPPA')
#' org.Hs.egHPPA <- xRDataLoader(RData.customised='org.Hs.egHPPA')
#' org.Hs.egHPPA <- xRDataLoader(RData.customised='org.Hs.egHPPA')
#' 
#' # from OSF
#' org.Mm.egKEGG <- xRDataLoader(RData='org.Mm.egKEGG', guid='gskpn')
#' org.Mm.string_high <- xRDataLoader(RData='org.Mm.string_high', guid='gskpn')
#' }

xRDataLoader <- function(RData=c(NA,"GWAS2EF", "GWAS_LD", "IlluminaHumanHT", "IlluminaOmniExpress", "ig.DO", "ig.EF", "ig.GOBP", "ig.GOCC", "ig.GOMF", "ig.HPCM", "ig.HPMA", "ig.HPMI", "ig.HPPA", "ig.MP", "org.Hs.eg", "org.Hs.egDGIdb", "org.Hs.egDO", "org.Hs.egGOBP", "org.Hs.egGOCC", "org.Hs.egGOMF", "org.Hs.egHPCM", "org.Hs.egHPMA", "org.Hs.egHPMI", "org.Hs.egHPPA", "org.Hs.egMP", "org.Hs.egMsigdbC1", "org.Hs.egMsigdbC2BIOCARTA", "org.Hs.egMsigdbC2CGP", "org.Hs.egMsigdbC2CPall", "org.Hs.egMsigdbC2CP", "org.Hs.egMsigdbC2KEGG", "org.Hs.egMsigdbC2REACTOME", "org.Hs.egMsigdbC3MIR", "org.Hs.egMsigdbC3TFT", "org.Hs.egMsigdbC4CGN", "org.Hs.egMsigdbC4CM", "org.Hs.egMsigdbC5BP", "org.Hs.egMsigdbC5CC", "org.Hs.egMsigdbC5MF", "org.Hs.egMsigdbC6", "org.Hs.egMsigdbC7", "org.Hs.egMsigdbH", "org.Hs.egPS", "org.Hs.egSF", "org.Hs.egPfam", "org.Hs.string", "org.Hs.PCommons_DN", "org.Hs.PCommons_UN"), RData.customised=NULL, verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{
	
    startT <- Sys.time()
    if(verbose){
    	message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    	message("", appendLF=TRUE)
    }
	
	RData <- RData[1]
    if(is.na(RData) & !is.null(RData.customised)){
		RData <- RData.customised
	}else if(is.na(RData) & is.null(RData.customised)){
		stop("There is no input! Please input one of two parameters ('RData' or 'RData.customised').\n")
	}
	RData <- gsub('.RData$', "", RData, ignore.case=TRUE, perl=TRUE)
	RData <- gsub(".RDa$", "", RData, ignore.case=TRUE, perl=TRUE)
	
	######################################################################################
	# obtain from Open Science Frame (OSF)
	######################################################################################
	flag_osf <- FALSE
	# check in order: 
	# 1) whether 5-digit guid (global unique identifier, eg 'gskpn') is provided
	# 2) wether the package 'osfr' is installed
	# 3) whether provided guid (for a project on OSF) can be retrieved (via osfr::osf_retrieve_node)
	# 4) whether to-be-queried RData file is there (via osfr::osf_ls_files)
	if(!is.null(guid) && nchar(guid)==5){
		pkgs <- c("osfr")
    	if(all(pkgs %in% rownames(utils::installed.packages()))){
        	tmp <- sapply(pkgs, function(pkg) {
            	requireNamespace(pkg, quietly=TRUE)
        	})
        	if(all(tmp)){
        		
        		######################################
				## temporarily mask the package "osfr"
				prj <- fls <- res <- NULL
				if(!is(suppressWarnings(try(prj<-osfr::osf_retrieve_node(guid), TRUE)), "try-error")){
					target <- paste0(RData,".RData")
					fls <- osfr::osf_ls_files(prj, type="file", pattern=target, n_max=Inf)
					if(nrow(fls)>0){
						ind <- match(fls$name, target)
						ind <- ind[!is.na(ind)]
						if(length(ind)==1){
							fl <- fls[ind,]
							
							res <- fl %>% osfr::osf_download(path=tempdir(),conflicts='overwrite')
							#res %>% osf_open()
							# verify the file downloaded locally
							if(file.exists(res$local_path)){
								out <- get(load(res$local_path))
								load_RData <- sprintf("'%s' at %s", prj$name, paste0('https://osf.io/',prj$id))
								RData <- target
								flag_osf <- TRUE
							}
						
						}
					}
				}
				######################################

			}
        }
	}
	
	######################################################################################
	# define function: my_https_downloader
	######################################################################################
	my_https_downloader <- function (url, method=c("auto","internal","wininet","libcurl","wget","curl"), quiet=TRUE, mode=c("w","wb","a","ab"), cacheOK=TRUE, extra=getOption("download.file.extra")){
	
		## https://stat.ethz.ch/R-manual/R-devel/library/utils/html/download.file.html
		method <- match.arg(method)
		mode <- match.arg(mode)
	
		## specify the temporary image files
		tdir <- tempdir()
		destfile <- file.path(tdir, "temp.RData")
		## remove the existing temporary RData file
		unlink(destfile, recursive=TRUE, force=TRUE)
	
		if(base::grepl("^https?://", url)){
			isR32 <- base::getRversion() >= "3.2"
			if(.Platform$OS.type == "windows"){
				if(isR32){
					method <- "wininet"
				}else{
					seti2 <- utils::"setInternet2"
					internet2_start <- seti2(NA)
					if(!internet2_start){
						on.exit(suppressWarnings(seti2(internet2_start)))
						suppressWarnings(seti2(TRUE))
					}
					method <- "internal"
				}
				#suppressWarnings(utils::download.file(url, destfile=destfile, method=method, quiet=quiet, mode=mode, cacheOK=cacheOK, extra=extra))
			}else{
				if(isR32 && capabilities("libcurl")){
					method <- "libcurl"
				}else if(nzchar(Sys.which("wget")[1])){
					method <- "wget"
				}else if(nzchar(Sys.which("curl")[1])){
					method <- "curl"
					orig_extra_options <- getOption("download.file.extra")
					on.exit(options(download.file.extra = orig_extra_options))
					options(download.file.extra = paste("-L", orig_extra_options))
				}else if(nzchar(Sys.which("lynx")[1])) {
					method <- "lynx"
				}else{
					stop("no download method found")
				}
				#suppressWarnings(utils::download.file(url, destfile=destfile, method=method, quiet=quiet, mode=mode, cacheOK=cacheOK, extra=extra))
			}
		}else{
			#suppressWarnings(utils::download.file(url, destfile=destfile, method=method, quiet=quiet, mode=mode, cacheOK=cacheOK, extra=extra))
		}
		
		if(is(suppressWarnings(try(utils::download.file(url, destfile=destfile, method=method, quiet=quiet, mode=mode, cacheOK=cacheOK, extra=extra), TRUE)),"try-error")){
			res_RData <- NULL
			res_flag <- FALSE
		}
		
		if(file.exists(destfile) & file.info(destfile)$size!=0){
		
			if(is(suppressWarnings(try(load(destfile), TRUE)),"try-error")){
				res_RData <- NULL
				res_flag <- FALSE				
			}else{
				res_RData <- get(load(destfile))
				res_flag <- TRUE
			}
			
		}else{
			res_RData <- NULL
			res_flag <- FALSE
		}
		
		res <- list(RData = res_RData,
					flag = res_flag)
		
		invisible(res)
	}
	######################################################################################
	
	######################################################################################	
	## obtain locally or remotely (other than OSF)
	######################################################################################
	if(!flag_osf){
	
		###############################
		## make sure there is no "/" at the end
		path_host <- gsub("/$", "", RData.location)
		if(path_host=="" || length(path_host)==0 || is.na(path_host)){
			path_host <- "https://github.com/hfang-bristol/RDataCentre/blob/master/Portal"
		}
	
		## load 
		load_remote <- paste(path_host, "/", RData, ".RData", sep="")
		load_local1 <- file.path(path_host, paste("data/", RData, ".RData", sep=""))
		load_local2 <- file.path(path_host, paste(RData, ".RData", sep=""))
		load_package <- RData
	
		#####################################################################
		## first, load data from the package itself (NOW DISABLE THIS OPTION)
		#####################################################################
		#if(length(suppressWarnings(tryCatch(eval(parse(text=paste("data(",load_package,", package='XGR')",sep=""))), error=function(e) e, warning=function(w) w)))==2){
	
		if(1){
			## second, load local R files
			RData_local <- c(load_local1, load_local2)
			load_flag <- sapply(RData_local, function(x){
				if(.Platform$OS.type=="windows") x <- gsub("/", "\\\\", x)
				ifelse(file.exists(x), TRUE, FALSE)
			})
			## otherwise, load remote R files
			if(sum(load_flag)==0){
			
				flag_failed <- FALSE
				if(length(grep('^https',load_remote,perl=TRUE))){
					if(length(grep('github',load_remote,perl=TRUE))){
						load_remote <- paste(load_remote, "?raw=true", sep="")
					}
					res <- my_https_downloader(load_remote, mode="wb")
					if(res$flag==FALSE){
						flag_failed <- TRUE
					}else{
						eval(parse(text=paste(RData, " <- res$RData", sep="")))
					}
				}else{
					res <- my_https_downloader(load_remote, mode="wb")
					if(res$flag==FALSE){
						flag_failed <- TRUE
					}else{
						eval(parse(text=paste(RData, " <- res$RData", sep="")))
					}
				}
			
				if(flag_failed){
			
					load_remotes <- c(
					paste("https://github.com/hfang-bristol/RDataCentre/blob/master/Portal/", RData, ".RData?raw=true", sep=""),
					paste("http://galahad.well.ox.ac.uk/bigdata/", RData, ".RData", sep=""),
					paste("http://galahad.well.ox.ac.uk/bigdata/", RData, ".RData", sep="")
					)
				
					for(i in 1:length(load_remotes)){
						load_remote <- load_remotes[i]
						if(verbose){
							now <- Sys.time()
							message(sprintf("Attempt to download from %s (at %s)", load_remote, as.character(now)), appendLF=TRUE)
						}
						res <- my_https_downloader(load_remote, mode="wb")
						if(res$flag==TRUE){
							break
						}
					}
				
					if(res$flag==FALSE){
						warnings("Built-in Rdata files cannot be loaded. Please check your internet connection or their location in your local machine.\n")
						eval(parse(text=paste(RData, " <- res$RData", sep="")))
					}else{
						eval(parse(text=paste(RData, " <- res$RData", sep="")))
					}
				}
		
				load_RData <- load_remote
				out <- base::get(RData)
			
			}else{
				load_RData <- RData_local[load_flag]
				out <- base::get(load(load_RData))
			}
		}else{
			load_RData <- sprintf("package 'XGR' version %s", utils::packageVersion("XGR"))
			out <- base::get(RData)
		}
	
	}
	
	####################################################################################
	
    if(verbose){
        now <- Sys.time()
        if(!is.null(out)){
			message(sprintf("'%s' (from %s) has been loaded into the working environment (at %s)", RData, load_RData, as.character(now)), appendLF=TRUE)
		}else{
			message(sprintf("'%s' CANNOT be loaded (at %s)", RData, as.character(now)), appendLF=TRUE)
		}
    }
    
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(verbose){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    invisible(out)
}
