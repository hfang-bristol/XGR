#' Function to visualise tuning parameters
#'
#' \code{xClassifyPara} is supposed to visualise tuning parameters. 
#'
#' @param sClass an object of class "sClass". Alternatively, it can be a data frame
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param contour logical to indicate whether 2D coutour plot should be also included
#' @param plot.3D logical to indicate whether to plot in 3D
#' @param theta.3D the azimuthal direction. By default, it is 40
#' @param phi.3D the colatitude direction. By default, it is 20
#' @param text.3D logical whether to text the point with the maximum value
#' @return 
#' a ggplot object (1D), a trellis object (2D), and NULl (3D)
#' @note none
#' @export
#' @seealso \code{\link{xClassifyPara}}
#' @include xClassifyPara.r
#' @examples
#' \dontrun{
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' gp <- xClassifyPara(sClass)
#' xClassifyPara(sClass, plot.3D=T, theta.3D=40, phi.3D=45)
#' }

xClassifyPara <-function(sClass, colormap="spectral", ncolors=50, contour=TRUE, plot.3D=FALSE, theta.3D=40, phi.3D=25, text.3D=TRUE)
{

    if(class(sClass) == "sClass"){
    	if(is.null(sClass$parameter)){
			warnings("The function must apply to a 'sClass' object with data frame 'parameter'.\n")
			return(NULL)
    	}else{
    		if(ncol(sClass$parameter)==3){
    			warnings("The function must apply to a 'sClass' object with data frame 'parameter'.\n")
				return(NULL)
    		}else{
    			df_fit <- sClass$parameter[,c(-2,-3)]
    		}
    	}
        
    }else if(class(sClass)=="data.frame" & nD!="auto"){
		df_fit <- sClass
		
    }else{
    	warnings("The function must apply to a 'sClass' object or data.frame.\n")
		return(NULL)
    }
	
	## dimension
	if(ncol(df_fit)==2){
		nD <- "1D"
		df <- data.frame(Parameter=df_fit[,2], Val=df_fit[,1], stringsAsFactors=FALSE)
		## get parameter name
		xlab <- colnames(df_fit)[1]
		ylab <- colnames(df_fit)[2]
		
	}else if(ncol(df_fit)==3){
		nD <- "2D"
		if(plot.3D){
			nD <- "3D"
		}
		## get parameter name
		xlab <- colnames(df_fit)[2]
		ylab <- colnames(df_fit)[3]
		zlab <- colnames(df_fit)[1]
		
	}
	
	## plot
	if(nD=="1D"){
		## define levels
		df$Parameter <- factor(df$Parameter, levels=rev(unique(df$Parameter)))
		
		Parameter <- Val <- NULL
		
		p <- ggplot(df, aes(x=Val, y=Parameter))
		p <- p + geom_point(fill='transparent',shape=21)
		
		# two lines
		ind <- which(df$Val==max(df$Val))
		p <- p + geom_vline(xintercept=df$Val[ind], color='black', linetype='dashed')
		p <- p + geom_hline(yintercept=as.numeric(df$Parameter[ind]), color='black', linetype='dashed')
		
		p <- p  + theme_bw() + theme(legend.position="none")
		p <- p + xlab(xlab) + ylab(ylab)
		p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
		## put arrows on both axes
		p <- p + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
		
		## x-axis position
    	res <- p + scale_x_continuous(position="top")
		
	}else if(nD=="2D"){
		## prepare data frame (x, y, z)
		x <- as.factor(df_fit[,2])
		y <- as.factor(df_fit[,3])
		z <- df_fit[,1]
		df_xyz <- data.frame(x=x, y=y, z=z)
		
		if(contour){
			res <- lattice::contourplot(z~x*y, df_xyz, xlab=xlab, ylab=ylab, scales=list(x=list(rot=90,cex=0.8)), region=TRUE, pretty=TRUE, col.regions=xColormap(colormap)(ncolors), colorkey=list(tck=0.8))
			
		}else{
			res <- lattice::levelplot(z~x*y, df_xyz, xlab=xlab, ylab=ylab, scales=list(x=list(rot=90,cex=0.8)), region=TRUE, pretty=TRUE, col.regions=xColormap(colormap)(ncolors), colorkey=list(tck=0.8))
			
		}
		
	}else if(nD=="3D"){
		
		## prepare data colvar matrix
		data_colvar <- as.matrix(xSparseMatrix(df_fit[,c(2,3,1)], verbose=FALSE))
	
		## prepare data matrix
		data <- as.matrix(xSparseMatrix(df_fit[,c(2,3,1)], verbose=FALSE))
		rows <- as.character(unique(df_fit[,2]))
		columns <- as.character(unique(df_fit[,3]))
		ind_rows <- match(rows, rownames(data))
		ind_columns <- match(columns, colnames(data))
		data <- data[ind_rows, ind_columns]
		data_colvar <- data_colvar[ind_rows, ind_columns]
		
		## x- and y-axis
		x_at <- 1:nrow(data)
		y_at <- 1:ncol(data)
		
		# highlight
		hightlight_ind <- which(data_colvar == max(data_colvar), arr.ind=TRUE)
		hightlight_x <- x_at[hightlight_ind[, 1]]
		hightlight_y <- y_at[hightlight_ind[, 2]]
		hightlight_zval <- round(data[hightlight_ind], 3)
		hightlight_xval <- rownames(data)[hightlight_x]
		hightlight_yval <- colnames(data)[hightlight_y]
		
		if(text.3D){
			plot <- FALSE
		}else{
			plot <- TRUE
		}
		
		## 3D plot
		zlim <- c(floor(min(data)*100)/100, ceiling(max(data)*100)/100)
		zlim <- c(1.5*zlim[1]-0.5*zlim[2], zlim[2])
		#zlim <- c(0.5, zlim[2])
		plot3D::persp3D(z=data, x=x_at, y=y_at, colvar=data_colvar, axes=TRUE, zlim=zlim, cex.axis=0.8, cex.lab=1.2, ticktype="simple", col=xColormap(colormap)(ncolors), colkey=list(side=4,length=0.3,width=0.6,shift=0.2,cex.axis=0.6,cex.clab=0.8,side.clab=2,tck=-0.3), xlab=paste0("\n",xlab), ylab=paste0("\n",ylab), zlab=paste0("\n",zlab), clab=c("","","",zlab), bty="b", facets=TRUE, curtain=FALSE, plot=plot, lighting=FALSE, lphi=90, theta=theta.3D, phi=phi.3D, contour=list(col="grey",lwd=0.8))
		if(text.3D){
			plot3D::scatter3D(x=hightlight_x, y=hightlight_y, z=zlim[1], type="h", colkey=FALSE, pch=19, cex=0.6, alpha=0.5, col="black", add=TRUE, plot=FALSE)
			label <- paste0(zlab,"=",hightlight_zval,"\n(",xlab,"=",hightlight_xval,")\n(",ylab,"=",hightlight_yval,")")
			plot3D::text3D(x=hightlight_x+0.5, y=hightlight_y, z=zlim[1], label=label, colkey=FALSE, cex=0.8, col="red", srt=0, add=TRUE, plot=TRUE)
		}
		
		res <- NULL
	}

    invisible(res)
}

    
