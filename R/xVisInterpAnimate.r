#' Function to animate the visualisation of interpolated irregular data
#'
#'
#' \code{xVisInterpAnimate} is supposed to animate the visualisation of interpolated irregular data. The output can be a pdf file containing a list of frames/images, a mp4 video file or a gif file. To support video output file, the software 'ffmpeg' must be first installed (also put its path into the system PATH variable; see Note). To support gif output file, the software 'ImageMagick' must be first installed (also put its path into the system PATH variable; see Note).
#'
#' @param ls_xyz a list with 3 required components (x, y and z) and an optional component (label)
#' @param interpolation the method for the interpolation. It can be "linear" or "spline" interpolation
#' @param nx the dimension of output grid in x direction
#' @param ny the dimension of output grid in y direction
#' @param zlim the minimum and maximum z values, defaulting to the range of the finite values of z
#' @param colkey a logical (TRUE by default) or a 'list' with parameters for the color key (legend). List parameters should be one of 'side, plot, length, width, dist, shift, addlines, col.clab, cex.clab, side.clab, line.clab, adj.clab, font.clab'. The defaults for the parameters are 'side=4, plot=TRUE, length=1, width=1, dist=0, shift=0, addlines=FALSE, col.clab=NULL, cex.clab=par("cex.lab"), side.clab=NULL, line.clab=NULL, adj.clab=NULL, font.clab=NULL'. colkey=list(side=4,length=0.15,width=0.5,shift=0.35,dist=-0.15,cex.axis=0.6,cex.clab=0.8,side.clab=3)
#' @param contour a logical (FALSE by default) or a 'list' with parameters for the contour function. An optional parameter to this 'list' is the 'side' where the image should be plotted. Allowed values for 'side' are a z-value, or 'side = "zmin", "zmax"', for positioning at bottom or top respectively. The default is to put the image at the bottom
#' @param image a logical (FALSE by default) or a 'list' with parameters for the image2D function. An optional parameter to this 'list' is the 'side' where the image should be plotted. Allowed values for 'side' are a z-value, or 'side = "zmin", "zmax"', for positioning at bottom or top respectively. The default is to put the image at the bottom
#' @param clab a title for the colorbar. the label to be written on top of the color key; to lower it, 'clab' can be made a vector, with the first values empty strings.
#' @param nlevels the number of levels to partition the input matrix values. The same level has the same color mapped to
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param label.pch a numeric value specifying the graphiics symbol (by default, 17 for upward triangle). This argument only works when the labelling is enabled
#' @param label.text.cex a numeric value specifying the text size. This argument only works when the labelling is enabled
#' @param xy.swap logical to indicate whether to wrap x and y. By default, it sets to false
#' @param theta.3D the starting azimuthal direction. By default, it is 0
#' @param phi.3D the colatitude direction. By default, it is 20
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param filename the without-extension part of the name of the output file. By default, it is 'xVisInterpAnimate'
#' @param filetype the type of the output file, i.e. the extension of the output file name. It can be one of either 'pdf' for the pdf file, 'mp4' for the mp4 video file, 'gif' for the gif file
#' @param image.type the type of the image files temporarily generated. It can be one of either 'jpg' or 'png'. These temporary image files are used for producing mp4/gif output file. The reason doing so is to accommodate that sometimes only one of image types is supported so that you can choose the right one
#' @param image.bg the background color for each frame/image. This argument only works when producing mp4 video or gif file
#' @param height.device a numeric value specifying the height (or width) of device/frame/image.
#' @param num.frame a numeric value specifying the number of frames/images. By default, it sets to the number of columns in the input data matrix
#' @param sec_per_frame a numeric value specifying how long (seconds) it takes to stream a frame/image. This argument only works when producing mp4 video or gif file
#' @param res the resolution for each frame/image. This argument only works when producing mp4 video or gif file
#' @return 
#' If specifying the output file name (see argument 'filename' above), the output file is either 'filename.pdf' or 'filename.mp4' or 'filename.gif' in the current working directory. If no output file name specified, by default the output file is either 'xVisInterpAnimate.pdf' or 'xVisInterpAnimate.mp4' or 'xVisInterpAnimate.gif'
#' @note When producing mp4 video, this function requires the installation of the software 'ffmpeg' at \url{https://www.ffmpeg.org}. Shell command lines for ffmpeg installation in Terminal (for both Linux and Mac) are:
#' \itemize{
#' \item{1) \code{wget -O ffmpeg.tar.gz http://www.ffmpeg.org/releases/ffmpeg-2.7.1.tar.gz}}
#' \item{2) \code{mkdir ~/ffmpeg | tar xvfz ffmpeg.tar.gz -C ~/ffmpeg --strip-components=1}}
#' \item{3) \code{cd ffmpeg}}
#' \item{4a) # Assuming you want installation with a ROOT (sudo) privilege: \cr\code{./configure --disable-yasm}}
#' \item{4b) # Assuming you want local installation without ROOT (sudo) privilege: \cr\code{./configure --disable-yasm --prefix=$HOME/ffmpeg}}
#' \item{5) \code{make}}
#' \item{6) \code{make install}}
#' \item{7) # add the system PATH variable to your ~/.bash_profile file if you follow 4b) route: \cr\code{export PATH=$HOME/ffmpeg:$PATH}}
#' \item{8) # make sure ffmpeg has been installed successfully: \cr\code{ffmpeg -h}}
#' }
#' When producing gif file, this function requires the installation of the software 'ImageMagick' at \url{http://www.imagemagick.org}. Shell command lines for ImageMagick installation in Terminal are:
#' \itemize{
#' \item{1) \code{wget http://www.imagemagick.org/download/ImageMagick.tar.gz}}
#' \item{2) \code{mkdir ~/ImageMagick | tar xvzf ImageMagick.tar.gz -C ~/ImageMagick --strip-components=1}}
#' \item{3) \code{cd ImageMagick}}
#' \item{4) \code{./configure --prefix=$HOME/ImageMagick}}
#' \item{5) \code{make}}
#' \item{6) \code{make install}}
#' \item{7) # add the system PATH variable to your ~/.bash_profile file. \cr For Linux: \cr\code{export MAGICK_HOME=$HOME/ImageMagick} \cr\code{export PATH=$MAGICK_HOME/bin:$PATH} \cr\code{export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:+$LD_LIBRARY_PATH:}$MAGICK_HOME/lib} \cr For Mac: \cr\code{export MAGICK_HOME=$HOME/ImageMagick} \cr\code{export PATH=$MAGICK_HOME/bin:$PATH} \cr\code{export DYLD_LIBRARY_PATH=$MAGICK_HOME/lib/}}
#' \item{8a) # check configuration: \cr\code{convert -list configure}}
#' \item{8b) # check image format supported: \cr\code{identify -list format}}
#' \item{Tips: \cr Prior to 4), please make sure \code{libjpeg} and \code{libpng} are installed. If NOT, for Mac try this: \cr\code{brew install libjpeg libpng} \cr To check whether ImageMagick does work, please get additional information from: \cr\code{identify -list format} \cr\code{convert -list configure} \cr On details, please refer to \url{http://www.imagemagick.org/script/advanced-unix-installation.php}}
#' }
#' @export
#' @seealso \code{\link{visNetMul}}
#' @include xVisInterpAnimate.r
#' @examples
#' \dontrun{
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' g <- erdos.renyi.game(20, 1/10)
#' glayout <- layout_with_kk(g)
#' ls_xyz <- data.frame(x=glayout[,1], y=glayout[,2], z=degree(g), label=degree(g))
#' 
#' # 3D views of different angles
#' # output as a pdf file
#' xVisInterpAnimate(ls_xyz, image=TRUE, filetype="pdf")
#' # output as a mp4 file
#' xVisInterpAnimate(ls_xyz, filetype="mp4")
#' # output as a gif file
#' xVisInterpAnimate(ls_xyz, filetype="gif", num.frame=72, sec_per_frame=0.5)
#' }

xVisInterpAnimate <- function (ls_xyz, interpolation=c("spline","linear"), nx=100, ny=100, zlim=NULL, colkey=TRUE, contour=FALSE, image=FALSE, clab=c("Value"), nlevels=20, colormap="terrain", label.pch=17, label.text.cex=0.8, xy.swap=FALSE, theta.3D=0, phi.3D=20, verbose=TRUE, filename="xVisInterpAnimate", filetype=c("pdf", "mp4", "gif"), image.type=c("jpg","png"), image.bg="transparent", height.device=NULL, num.frame=36, sec_per_frame=1, res=72)
{
    
    interpolation <- match.arg(interpolation)
    filetype <- match.arg(filetype)
    image.type <- match.arg(image.type)
    
    if(is.null(filename)){
        outputfile <- paste("xVisInterpAnimate", filetype, sep=".")
    }else{
        outputfile <- paste(filename, filetype, sep=".")
    }

    if(filetype=="pdf"){
        
        if(is.null(height.device)){
            height.device <- 7
        }
        
        grDevices::pdf(outputfile, width=height.device, height=height.device, bg=image.bg)
        for(theta.3D.move in seq(theta.3D, 360+theta.3D, 360/num.frame)){ 
        	theta.3D.move <- theta.3D.move %% 360
        	
			if(verbose){
				now <- Sys.time()
				message(sprintf("Viewed at %d polar angle and %d azimuthal angle (%s)", phi.3D, theta.3D.move, as.character(now)), appendLF=TRUE)
			}
        	
        	clab.move <- paste0("Viewer's angle (polar:",phi.3D,", azimuthal:",theta.3D.move,")\n",clab,"\n\n\n")
        	xVisInterp(ls_xyz=ls_xyz, interpolation=interpolation, nx=nx, ny=ny, zlim=zlim, nD="3D", colkey=colkey, contour=contour, image=image, clab=clab.move, nlevels=nlevels, colormap=colormap, label.pch=label.pch, label.text.cex=label.text.cex, xy.swap=xy.swap, theta.3D=theta.3D.move, phi.3D=phi.3D, verbose=FALSE)
        }
        grDevices::dev.off()
        
        if(file.exists(file.path(getwd(), outputfile))){
            message(sprintf("Congratulations! A file '%s' (in the directory %s) has been created!", outputfile, getwd()), appendLF=T)
        }
        
        invisible()
        
    }else if(filetype=="mp4" | filetype=="gif"){
        
        ## num.frame: how many frames in total
        ## sec_per_frame: seconds per frame
        ## frame_per_sec: frames per second
        frame_per_sec <- 1/sec_per_frame
        
        ## specify the temporary image files
        tdir <- tempdir()
        if(image.type=='png'){
        	image_files <- file.path(tdir, "Rplot%06d.png")
        	## remove the existing temporary png files
        	unlink(file.path(tdir, "Rplot*.png"), recursive=T, force=T)
        }else if(image.type=='jpg'){
        	image_files <- file.path(tdir, "Rplot%06d.jpg")
        	## remove the existing temporary jpg files
        	unlink(file.path(tdir, "Rplot*.jpg"), recursive=T, force=T)
        }
        tfile <- file.path(tdir, base::basename(outputfile))
        unlink(tfile, recursive=T, force=T)
        
        if(is.null(height.device)){
            height.device <- 480*res/72
        }
        
        if(image.type=='png'){
        	grDevices::png(image_files, width=height.device, height=height.device, res=res, bg=image.bg)
        }else if(image.type=='jpg'){
        	grDevices::jpeg(image_files, width=height.device, height=height.device, res=res, bg=image.bg)
        }
        for(theta.3D.move in seq(theta.3D, 360+theta.3D, 360/num.frame)){ 
        	theta.3D.move <- theta.3D.move %% 360
        	
			if(verbose){
				now <- Sys.time()
				message(sprintf("Viewed at %d polar angle and %d azimuthal angle (%s)", phi.3D, theta.3D.move, as.character(now)), appendLF=TRUE)
			}
        	
        	clab.move <- paste0("Viewer's angle (polar:",phi.3D,", azimuthal:",theta.3D.move,")\n",clab,"\n\n\n")
        	xVisInterp(ls_xyz=ls_xyz, interpolation=interpolation, nx=nx, ny=ny, zlim=zlim, nD="3D", colkey=colkey, contour=contour, image=image, label.pch=label.pch, label.text.cex=label.text.cex, xy.swap=xy.swap, clab=clab.move, nlevels=nlevels, colormap=colormap, theta.3D=theta.3D.move, phi.3D=phi.3D, verbose=FALSE)
        }
        grDevices::dev.off()
        
        if(filetype=="mp4"){
			ffmpeg1 <- paste("ffmpeg -y -v quiet -r", frame_per_sec, "-i", image_files, "-q:v 1", file.path(tdir, outputfile))
			ffmpeg2 <- paste("$HOME/ffmpeg -y -v quiet -r", frame_per_sec, "-i", image_files, "-q:v 1", file.path(tdir, outputfile))
			ffmpeg_local <- c(ffmpeg1, ffmpeg2)
			cmd_flag <- 1
			for(i in 1:length(ffmpeg_local)){
				cmd <- try(system(ffmpeg_local[i]), silent=TRUE)
				if(cmd==0){
					cmd_flag <- 0
					message(sprintf("Executing this command: '%s'\n", ffmpeg_local[i]), appendLF=T)
					break
				}
			}
			
		}else if(filetype=="gif"){
		
			## http://www.r-bloggers.com/animate-gif-images-in-r-imagemagick/
			## -delay ticks: '100 ticks' corresponds to 1 second
			## ticks/100: seconds per image/frame
			## 100/ticks: images/frames per second
		
			image_files <- paste('Rplot','*.', image.type, sep='')
			convert1 <- paste("convert -delay", 100*sec_per_frame, file.path(tdir, image_files), file.path(tdir, base::basename(outputfile)))
			convert2 <- paste("$HOME/ImageMagick/bin/convert -delay", 100*sec_per_frame, file.path(tdir, image_files), tfile)
			convert_local <- c(convert1, convert2)
			cmd_flag <- 1
			for(i in 1:length(convert_local)){
				cmd <- try(system(convert_local[i]), silent=TRUE)
				if(cmd==0){
					cmd_flag <- 0
					message(sprintf("Executing this command: '%s'\n", convert_local[i]), appendLF=T)
					break
				}
			}
			
		}
		
        if(cmd_flag==0){
            if(file.exists(tfile)){
                file.copy(from=tfile, to=outputfile, overwrite=T, recursive=F, copy.mode=T)
                message(sprintf("Congratulations! A file '%s' (in the directory %s) has been created!", outputfile, getwd()), appendLF=T)
            }
        }else{
            stop("Unfortunately, fail to produce the file. Please install ffmpeg or ImageMagick first. Also make sure its path being put into the system PATH variable (see Help). Alternatively, produce the pdf file instead\n")
        }

        invisible(cmd)
    }
    
}
