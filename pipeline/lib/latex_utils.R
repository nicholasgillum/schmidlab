#Sweave and latex
latex_escape <- function(filename) {
	pretty = paste("\\\\verb|", filename,"|")
	return(pretty);
}

sanatize<-function(dirty){
	return(gsub('_', "\\\\\\\\_", dirty))
}

table <- function(filenames) {
	x ="Input Files ";
	for(i in 1:length(filenames)) {
		x=paste(x, " &", latex_escape(files[i]), "\\\\\\\\", sep="")
	}
	return(x);
}

getPairedMAFiguresBg <- function(tmpDir, config)  {
	drawCommand=""
	n=config$count
	if(config$quiet$backgroundCorrection){
		return(" ")
	}
	for(i in 0:(n-1)) {
		FILENAME=paste(tmpDir,"/maAfterBg", i, ".png", sep="" )
		drawCommand = paste(drawCommand, 
				"\\\\begin{figure}[htb!]\n",
				"\\\\centering\n",
				"\\\\includegraphics[width=6in]{", FILENAME, "}\n",
				"\\\\caption{MA plot of ", sanatize(config$targets$Name[i]),
				" after background correction}\n",
				"\\\\end{figure}\\\\pagebreak\n",
				sep="")
	}
	return(drawCommand)
}

getMAFiguresWAN <- function(tmpDir, config)  {
	if(config$vsn | config$quiet$withinArrayNormalization) {
		return(" ")
	}
	drawCommand=""
	for(i in 1:config$count) {
		FILENAME=paste(config$files$tmpDir,"/maAfterWAN", i-1, ".png", sep="" )
		drawCommand = paste(drawCommand, 
				"\\\\begin{figure}[htb!]\n",
				"\\\\centering",
				"\\\\includegraphics[width=6in]{", FILENAME, "}\n",
				"\\\\caption{", sanatize(config$targets$Name[i]), 
				" after within array normalization by ", 
				config$algorithms$withinArray,
				"}\n \\\\end{figure}\\\\pagebreak\n",
				sep="")
	}
	return(drawCommand)
}

getPairedMAFiguresBAN <- function(config)  {
	tmpDir = config$files$tmpDir
	n      = config$count
	method = config$algorithms$betweenArray
	
	drawCommand=""
	for(i in 0:(n-1)) {
		FILENAME=paste(tmpDir,"/maAfterBAN", i, ".png", sep="" )
		drawCommand = paste(drawCommand, 
				"\\\\begin{figure}[htb!]\n",
				"\\\\centering\n",
				"\\\\includegraphics[width=6in]{", FILENAME, "}\n",
				"\\\\caption{MA plot of ",sanatize(config$targets$Name[i]), 
				" after between array normalization by ", method, "}\n",
				"\\\\end{figure}\\\\pagebreak\n",
				sep="")
	}
	return(drawCommand)
}

#For the document
getBackgroundText <- function(config) {
	return(config$bgText$main)
}

getBgMaPlotText <- function(config) {
	if(config$quiet$backgroundCorrection) {
		return(config$bgText$quiet)
	} else {
		return(config$bgText$verbose)
	}
}


getBoxplotCommandOne <-function(config) {
	if(config$vsn) {return(" ")}
	return(paste(
					"\\\\begin{figure}[htb!]\n",
					"\\\\centering\n",
					"\\\\includegraphics{",	config$files$tmpDir, "/boxPlot1.png}\n",
					"\\\\caption{Boxplots of the signal intensities of each signal ",
					"channel of the microarrays. Raw data after background correction.",
					"}\\\\end{figure}\n\\\\pagebreak\n",sep=""));
}

getWithinArrayNormText <-function(config){
	if(config$vsn) {
		return(config$wanText$vsn)
	}
	return(gsub("xxx", config$algorithms$withinArray, config$wanText$other))
}

getBoxplotCommandTwo <-function(config) {
	if(config$vsn) {return(" ")}
	return(paste("\\\\begin{figure}[htb!]\n\\\\centering\n",
					"\\\\includegraphics{", config$files$tmpDir,
					"/boxPlot2.png}\n \\\\caption{Boxplots of the signal intensities ",
					"of each signal channel of the microarrays after within array ",
					"normalization.}\\\\end{figure}\\\\pagebreak",sep=""));
}

getWanText <- function(config) {
	if(config$quiet$withinArrayNormalization) {
		return(config$wanText$quiet)
	} else {
		return(config$wanText$verbose)
	}
}


#From carma web
if( !isGeneric("drawBoxplot") )
	setGeneric("drawBoxplot", function(x,new.plot=TRUE,xlim=NULL,at=0,exclude.flagged=FALSE,...)
				standardGeneric("drawBoxplot"))

setMethod("drawBoxplot","MadbSet",
		function(x,new.plot=TRUE,xlim=NULL,at=0,exclude.flagged=FALSE,...){
			if(exclude.flagged){
				W <- getWeights(x)
				for(i in 1:ncol(exprs(x))){
					exprs(x)[W[,i]==0,i] <- NA
				}
			}
			boxplots(x=exprs(x),new.plot=new.plot,xlim=xlim,at=at,...)
		}
)

boxplots <- function(x,new.plot=TRUE,xlim=NULL,at=0,col=NULL,log2.transform=FALSE,ylim=NULL,...){
	data <- x
#	if(class(x)=="exprSet" | class(x)=="EexprSet"){
#		data <- exprs(x)
#	}
	if(log2.transform){
		data <- log2(data)
	}
	if(is.null(ncol(data))){
		data <- as.matrix(data)
	}
	## remove those nasty infinite values
	for(i in 1:ncol(data)){
		data[is.infinite(data[,i]),i] <- NA
	}
	
	CN <- colnames(data)
	if(is.null(CN)){
		CN <- 1:ncol(data)
	}
	if(is.null(col)){
		col=rep(0,ncol(data))
	}
	if(length(col)!=ncol(data)){
		col=rep(col[1],ncol(data))
	}
	if(new.plot){
		if(is.null(xlim)){
			xlim=c(0.5,(ncol(data)+0.5))
		}
		par(omi=0.7*c(min(2,max(strwidth(CN,units="i"))),0,0,0),cex.axis=0.7)
		if(is.null(ylim)){
			ylim=c(min(data,na.rm=TRUE),max(data,na.rm=TRUE))
		}
		plot(1,1,pch=NA,ylab="expression",xaxt="n",ylim=ylim,xlim=xlim,xlab=NA)
	}
	axis(side=1,labels=CN,at=(1+at):(ncol(data)+at),las=3)
	for(i in 1:ncol(data)){
		boxplot(data[,i],at=(i+at),add=TRUE,col=col[i],range=0,...)
	}
}