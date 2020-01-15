#' Function to read RDS files
#'
#' \code{xRDS} is supposed to read RDS files.
#'
#' @param RDS which RDS to load. To support the remote reading of a compressed RDS file, it must be compressed via the gzip method
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param placeholder the characters to tell the placeholder of RDS files
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. For example, 'gskpn' (see 'https://osf.io/gskpn'). If a valid provided and the query matched, it has priority over the one specified via placeholder
#' @return 
#' the loaded RDS. If the data cannot be loaded, it returns NULL.
#' @note To enable 'guid', please also install a package "osfr" via \code{BiocManager::install(c("remotes","centerforopenscience/osfr"),dependencies=T)}.
#' @export
#' @seealso \code{\link{xRDS}}
#' @include xRDS.r
#' @examples
#' \dontrun{
#' org.Hs.eg <- xRDS('org.Hs.eg')
#' ig.HPPA <- xRDS('ig.HPPA')
#' org.Hs.egHPPA <- xRDS('org.Hs.egHPPA')
#' 
#' # from OSF
#' org.Mm.egKEGG <- xRDS('org.Mm.egKEGG', guid='gskpn')
#' org.Mm.string_high <- xRDS('org.Mm.string_high', guid='gskpn')
#' }

xRDS <- function(RDS=NULL, verbose=T, placeholder=NULL, guid=NULL)
{
	
    startT <- Sys.time()
    if(verbose){
    	message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    	message("", appendLF=TRUE)
    }
	
    if(is.null(RDS)){
		stop("Please provide the RDS file.\n")
	}
	RDS <- gsub('.RDS$', "", RDS, ignore.case=T, perl=T)
	RDS <- gsub(".rds$", "", RDS, ignore.case=T, perl=T)
	
	######################################################################################
	# obtain from Open Science Frame (OSF)
	######################################################################################
	flag_osf <- F
	# check in order: 
	# 1) whether 5-digit guid (global unique identifier, eg 'gskpn') is provided
	# 2) wether the package 'osfr' is installed
	# 3) whether provided guid (for a project on OSF) can be retrieved (via osfr::osf_retrieve_node)
	# 4) whether to-be-queried RDS file is there (via osfr::osf_ls_files)
	if(!is.null(guid) && nchar(guid)==5){
		pkgs <- c("osfr")
    	if(all(pkgs %in% rownames(utils::installed.packages()))){
        	tmp <- sapply(pkgs, function(pkg) {
            	requireNamespace(pkg, quietly=T)
        	})
        	if(all(tmp)){
        		
        		######################################
				## temporarily mask the package "osfr"
				prj <- fls <- res <- NULL
				
				if(all(class(suppressWarnings(try(eval(parse(text=paste0('prj<-osfr::osf_retrieve_node(guid)'))), T))) != "try-error")){
					target <- paste0(RDS,".RDS")
					eval(parse(text=noquote(paste0('fls <- osfr::osf_ls_files(prj, type="file", pattern=target, n_max=Inf)'))))
					if(nrow(fls)>0){
						ind <- match(fls$name, target)
						ind <- ind[!is.na(ind)]
						if(length(ind)==1){
							fl <- fls[ind,]
						
							## specify the temporary file
							destfile <- file.path(tempdir(), fl$name)
							eval(parse(text=paste0('res <- fl %>% osfr::osf_download(overwrite=T, path=destfile)')))
							#res %>% osf_open()
							# verify the file downloaded locally
							if(file.exists(res$local_path)){
								out <- get(load(res$local_path))
								load_RDS <- sprintf("'%s' at %s", prj$name, paste0('https://osf.io/',prj$id))
								RDS <- target
								flag_osf <- T
							}
						
						}
					}
				}
				######################################

			}
        }
	}
	
	######################################################################################	
	## obtain locally or remotely (other than OSF)
	######################################################################################
	if(!flag_osf & !is.null(placeholder)){
		
		out <- NULL
		###############################
		## make sure there is no "/" at the end
		placeholder <- gsub("/$", "", placeholder)
		
		if(grepl("^https?://", placeholder)){
			load_remote <- paste0(placeholder, "/", RDS, ".RDS")
			if(class(suppressWarnings(try(out <- readRDS(gzcon(url(load_remote))), T)))=="try-error"){
				out <- NULL
			}else{
				load_RDS <- load_remote
			}

		}else{
			load_local <- file.path(placeholder, paste0(RDS, ".RDS"))
			if(.Platform$OS.type=="windows") load_local <- gsub("/", "\\\\", load_local)
			if(file.exists(load_local)){
				out <- readRDS(file.path(load_local))
				load_RDS <- load_local
			}
		}
	}
	
    if(verbose){
        if(!is.null(out)){
			message(sprintf("'%s' (from %s) has been loaded into the working environment (at %s)", RDS, load_RDS, as.character(Sys.time())), appendLF=T)
		}else{
			message(sprintf("'%s' CANNOT be loaded (at %s)", RDS, as.character(Sys.time())), appendLF=T)
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
