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
#' @note To enable 'guid', please also install a package "osfr" via \code{BiocManager::install("osfr",dependencies=TRUE)}.
#' @export
#' @importFrom osfr osf_retrieve_node osf_retrieve_node osf_download
#' @seealso \code{\link{xRDS}}
#' @include xRDS.r
#' @examples
#' guid <- 'gskpn'
#' \dontrun{
#' # from OSF
#' ig.HPPA <- xRDS('ig.HPPA', guid='gskpn')
#' }

xRDS <- function(RDS=NULL, verbose=TRUE, placeholder=NULL, guid=NULL)
{
	
    startT <- Sys.time()
    if(verbose){
    	message(sprintf("Starting ... (at %s)\n", as.character(startT)), appendLF=TRUE)
    }
	
    if(is.null(RDS)){
		stop("Please provide the RDS file.\n")
	}
	RDS <- gsub('.RDS$', "", RDS, ignore.case=TRUE, perl=TRUE)
	RDS <- gsub(".rds$", "", RDS, ignore.case=TRUE, perl=TRUE)
	
	out <- NULL
	
	######################################################################################
	# obtain from Open Science Frame (OSF)
	######################################################################################
	flag_osf <- FALSE
	# check in order: 
	# 1) whether 5-digit guid (global unique identifier, eg 'gskpn') is provided
	# 2) whether provided guid (for a project on OSF) can be retrieved (via osfr::osf_retrieve_node)
	# 3) whether to-be-queried RDS file is there (via osfr::osf_ls_files)
	if(!is.null(guid) && nchar(guid)==5){
		prj <- fls <- res <- NULL
		if(!is(suppressWarnings(try(prj<-osfr::osf_retrieve_node(guid), TRUE)),"try-error")){
			target <- paste0(RDS,".RDS")
			fls <- osfr::osf_ls_files(prj, type="file", pattern=target, n_max=Inf)
			if(nrow(fls)>0){
				ind <- match(fls$name, target)
				ind <- ind[!is.na(ind)]
				if(length(ind)==1){
					fl <- fls[ind,]
					res <- fl %>% osfr::osf_download(path=tempdir(), conflicts="overwrite")
					#res %>% osf_open()
					# verify the file downloaded locally
					if(file.exists(res$local_path)){
						out <- readRDS(res$local_path)
						load_RDS <- sprintf("'%s' at %s", prj$name, paste0('https://osf.io/',prj$id))
						RDS <- target
						flag_osf <- TRUE
					}
				}
			}
		}
	}
	
	######################################################################################	
	## obtain locally or remotely (other than OSF)
	######################################################################################
	if(!flag_osf & !is.null(placeholder)){
		
		###############################
		## make sure there is no "/" at the end
		placeholder <- gsub("/$", "", placeholder)
		
		if(grepl("^https?://", placeholder)){
			load_remote <- paste0(placeholder, "/", RDS, ".RDS")
			if(is(suppressWarnings(try(out <- readRDS(gzcon(url(load_remote))), TRUE)),"try-error")){
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
			message(sprintf("'%s' (from %s) successfully loaded (at %s)", RDS, load_RDS, as.character(Sys.time())), appendLF=TRUE)
		}else{
			message(sprintf("'%s' CANNOT be loaded (at %s)", RDS, as.character(Sys.time())), appendLF=TRUE)
		}
    }
    
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(verbose){
    	message(sprintf("\nEnded (at %s)", as.character(endT)), appendLF=TRUE)
    	message(sprintf("Runtime in total: %d secs\n", runTime), appendLF=TRUE)
    }
    
    invisible(out)
}
