#' Function to visualise interpolated irregular data
#'
#' \code{xVisInterp} is supposed to visualise irregular data after bilinear or bicubic spline interpolation onto a grid. 
#'
#' @param ls_xyz a list with three components (x, y and z)
#' @param interpolation the method for the interpolation. It can be "linear" or "spline" interpolation
#' @param nx the dimension of output grid in x direction
#' @param ny the dimension of output grid in y direction
#' @param zlim the minimum and maximum z values, defaulting to the range of the finite values of z
#' @param nD an integer specifying the dimension of the visualisation. It can be one of '2D', '3D', and 'auto' (to display the input raw data as well in 2D)
#' @param colkey a logical (TRUE by default) or a 'list' with parameters for the color key (legend). List parameters should be one of 'side, plot, length, width, dist, shift, addlines, col.clab, cex.clab, side.clab, line.clab, adj.clab, font.clab'. The defaults for the parameters are 'side=4, plot=TRUE, length=1, width=1, dist=0, shift=0, addlines=FALSE, col.clab=NULL, cex.clab=par("cex.lab"), side.clab=NULL, line.clab=NULL, adj.clab=NULL, font.clab=NULL'. colkey=list(side=4,length=0.15,width=0.5,shift=0.35,dist=-0.15,cex.axis=0.6,cex.clab=0.8,side.clab=3)
#' @param contour a logical (FALSE by default) or a 'list' with parameters for the contour function. An optional parameter to this 'list' is the 'side' where the image should be plotted. Allowed values for 'side' are a z-value, or 'side = "zmin", "zmax"', for positioning at bottom or top respectively. The default is to put the image at the bottom
#' @param image a logical (FALSE by default) or a 'list' with parameters for the image2D function. An optional parameter to this 'list' is the 'side' where the image should be plotted. Allowed values for 'side' are a z-value, or 'side = "zmin", "zmax"', for positioning at bottom or top respectively. The default is to put the image at the bottom
#' @param clab a title for the colorbar. the label to be written on top of the color key; to lower it, 'clab' can be made a vector, with the first values empty strings.
#' @param nlevels the number of levels to partition the input matrix values. The same level has the same color mapped to
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param theta.3D the azimuthal direction. By default, it is 40
#' @param phi.3D the colatitude direction. By default, it is 20
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return 
#' invisible
#' @note none
#' @export
#' @seealso \code{\link{xVisInterp}}
#' @include xVisInterp.r
#' @examples
#' \dontrun{
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' g <- erdos.renyi.game(100, 1/10)
#' glayout <- layout_with_kk(g)
#' ls_xyz <- data.frame(x=glayout[,1], y=glayout[,2], z=degree(g))
#' 
#' # auto
#' ls_xyz.smooth <- xVisInterp(ls_xyz, nD="auto")
#' # 2D
#' ls_xyz.smooth <- xVisInterp(ls_xyz, nD="2D")
#' # 3D
#' xVisInterp(ls_xyz, nD="3D", theta.3D=40, phi.3D=20, clab="Value\n")
#' # 3D views of different angles
#' pdf("xVisInterp.pdf")
#' for(theta.3D in seq(0,360,10)){ xVisInterp(ls_xyz, nD="3D", contour=TRUE, image=TRUE, clab=paste0("theta:",theta.3D,"\n"), theta.3D=theta.3D, phi.3D=20)}
#' dev.off()
#' }

xVisInterp <-function(ls_xyz, interpolation=c("spline","linear"), nx=100, ny=100, zlim=NULL, nD=c("auto","2D","3D"), colkey=TRUE, contour=FALSE, image=FALSE, clab=c("Value",""), nlevels=20, colormap="terrain",
theta.3D=40, phi.3D=20, verbose=TRUE)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    interpolation <- match.arg(interpolation)
    nD <- match.arg(nD)
	
	if(class(ls_xyz)=='list'){
		ind <- names(ls_xyz) %in% c("x","y","z")
		ls_xyz <- ls_xyz[ind]
		if(length(ls_xyz)!=3){
			return(NULL)
		}
		
    }else if(class(ls_xyz)=="data.frame" | class(ls_xyz)=="matrix"){
    	if(ncol(ls_xyz)==3){
			ls_xyz <- ls_xyz[,1:3]
			colnames(ls_xyz) <- c("x","y","z")
		}else if(ncol(ls_xyz)==4){
			ls_xyz <- ls_xyz[,1:4]
			colnames(ls_xyz) <- c("x","y","z","label")
		}
		ls_xyz <- as.list(ls_xyz)
		
    }else{
    	stop("The function must apply to a 'list' object, or a 'data.frame'/'matrix' object.\n")
    }
	
	if(interpolation=="linear"){
		linear <- TRUE
	}else{
		linear <- FALSE
	}
	
    ## increase smoothness (using finer grid)
    x <- y <- z <- label <- NULL
    ls_xyz.smooth <- with(ls_xyz, akima::interp(x, y, z, nx=nx, ny=ny, linear=linear))
	
	col <- xColormap(colormap)(nlevels)
	
	if(class(colkey)!='list'){
		if(colkey){
			colkey <- list(side=4,length=0.15,width=0.5,shift=0.35,dist=-0.15,cex.axis=0.6,cex.clab=0.8,side.clab=3)
		}
	}
	
	if(nD=='auto'){
		si.zmin <- min(ls_xyz.smooth$z, na.rm=TRUE)
		si.zmax <- max(ls_xyz.smooth$z, na.rm=TRUE)
		breaks <- pretty(c(si.zmin,si.zmax),nlevels+1)
		colors <- xColormap(colormap)(length(breaks)-1)
		
		if(image){
			#graphics::image(ls_xyz.smooth, axes=FALSE, breaks=breaks, col=colors)
			graphics::image(ls_xyz.smooth, axes=FALSE, breaks=breaks, col=colors, xlim=grDevices::extendrange(ls_xyz.smooth$x,f=0.1), ylim=grDevices::extendrange(ls_xyz.smooth$y,f=0.1))
		}
		
		if(contour){
			if(image){
				graphics::contour(ls_xyz.smooth, add=TRUE, labcex=1, levels=breaks, col="thistle")
			}else{
				graphics::contour(ls_xyz.smooth, labcex=1, levels=breaks, col=colors)
			}
			graphics::points(ls_xyz, pch=20, cex=1, col="black")
			
		}else{
			if(!image){
				plot(y ~ x, data=ls_xyz, pch=20, col="blue", axes=FALSE, ann=FALSE)
			}else{
				graphics::points(ls_xyz, pch=20, cex=1)
			}
			
			if(is.null(ls_xyz$label)){
				with(ls_xyz, graphics::text(x, y, formatC(z,dig=2), adj=-0.2, cex=0.6))
			}else{
				with(ls_xyz, graphics::text(x, y, label, adj=-0.2, cex=0.6))
			}
			
		}
		
		
	}else if(nD=='2D'){
		if(class(contour)!='list'){
			if(contour){
				contour <- list(col="black",labcex=1,lwd=2,alpha=0.5,addbox=FALSE,nlevels=nlevels)
			}
		}
		
		graphics::par(mar=c(0,0,0,0))
		plot3D::image2D(z=ls_xyz.smooth$z, axes=FALSE, cex.axis=0.8, cex.lab=1.2, col=col, colkey=colkey, clab=clab, plot=TRUE, lighting=FALSE, lphi=90, contour=contour)
    
    }else{
    	
		if(is.null(zlim)){
			zlim <- c(floor(min(ls_xyz.smooth$z,na.rm=TRUE)*100)/100, ceiling(max(ls_xyz.smooth$z,na.rm=TRUE)*100)/100)
			
			if(verbose){
				now <- Sys.time()
				message(sprintf("The range of interpolated values: [%.3f, %.3f]", zlim[1], zlim[2]), appendLF=TRUE)
			}
			
			if(class(contour)!='list' & class(image)!='list'){
				if(contour | image){
					zlim <- c(1.5*zlim[1]-0.5*zlim[2], zlim[2])
					#zlim <- c(2*zlim[1]-zlim[2], zlim[2])
				}
			}
		}
    	
    	if(class(contour)=='list' | class(image)=='list'){
    		#image=list(col=xColormap("terrain")(20),axes=FALSE)
   			#contour=list(side=c("zmin","z")[1],labcex=1,col=c("thistle","black")[2],nlevels=20,lwd=0.8,box=FALSE)
   			
			# 3D (easiest)
			graphics::par(mar=c(0,0,0,0))
			plot3D::persp3D(z=ls_xyz.smooth$z, axes=FALSE, box=FALSE, zlim=zlim, cex.axis=0.8, cex.lab=1.2, ticktype=c("simple","detailed")[1], col=col, colkey=colkey, clab=clab, bty="b", facets=TRUE, curtain=FALSE, plot=TRUE, lighting=TRUE, lphi=90, theta=theta.3D, phi=phi.3D, d=1, image=image, contour=contour)
			
    	}else{
    		
			if(contour | image){
				plot <- FALSE
			}else{
				plot <- TRUE
			}
		
			graphics::par(mar=c(0,0,0,0))
			plot3D::persp3D(z=ls_xyz.smooth$z, axes=FALSE, box=FALSE, zlim=zlim, cex.axis=0.8, cex.lab=1.2, ticktype=c("simple","detailed")[1], col=col, colkey=colkey, clab=clab, bty="b", facets=TRUE, curtain=FALSE, plot=plot, lighting=TRUE, lphi=90, theta=theta.3D, phi=phi.3D, d=1)
		
			if(image){
				if(contour){
					plot2 <- FALSE
				}else{
					plot2 <- TRUE
				}
				plot3D::image3D(z=zlim[1], colvar=ls_xyz.smooth$z, col=col, box=FALSE, colkey=FALSE, add=TRUE, plot=plot2)
		
			}
		
			if(contour){
				if(image){
					col <- c("thistle","grey")[2]
				}
				plot3D::contour3D(z=zlim[1], colvar=ls_xyz.smooth$z, labcex=1, col=col, nlevels=nlevels, lwd=0.8, box=FALSE, addbox=FALSE, colkey=FALSE, add=TRUE, plot=TRUE)
			}
		}
    }

    invisible(ls_xyz.smooth)
}

    
