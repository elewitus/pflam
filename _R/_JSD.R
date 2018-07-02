#read libraries
#(pracma)
###library(phytools)
##library(Matrix)
##library(PerformanceAnalytics)
##library(diversitree)
##library(RPANDA)
##library(ggplot2)
##library(gridExtra)
##library(TreePar)
##library(DDD)
##library(scatterplot3d)
##library(mclust)
##library(gplots)
##library(fpc)
##library(pvclust)
##library(geiger)
###library(classInt)
###library(TreeSim)
##library(igraph)#!
###library(Heatplus)
##library(RColorBrewer)
###library(vegan)
###library(cluster)
###library(cladoRcpp)
###library(phylobase)
###library(rexpokit)
###library(optimx)
###library(FD)
###library(xtable)
###library(BioGeoBEARS)#!
###library(clusterSim)
###library(ade4)
##library(TESS)#!
###library(iteRates)
##library(moments)
##library(MCMCpack)
##library(gtools)
###library(FactoMineR)
###library(corrgram)
###library(rfishbase)
#data(fishbase)

##gaussian kernel
sigma = 0.1
gKernel <- function(x) 1/(sigma*sqrt(2*pi)) * exp(-(x^2)/2*sigma^2)
kernelG <- function(x, mean=0, sd=1) dnorm(x, mean = mean, sd = sd)

##kernel density estimate
dens <- function(x, bw = bw.nrd0, kernel = kernelG, n = 4096,
                from = min(x) - 3*sd, to = max(x) + 3*sd, adjust = 1,
                ...) {
  if(has.na <- any(is.na(x))) {
    x <- na.omit(x)
    if(length(x) == 0)
        stop("no finite or non-missing data!")
  }
  sd <- (if(is.numeric(bw)) bw[1] else bw(x)) * adjust
  X <- seq(from, to, len = n)
  M <- outer(X, x, kernel, sd = sd, ...)
  structure(list(x = X, y = rowMeans(M), bw = sd,
                 call = match.call(), n = length(x),
                 data.name = deparse(substitute(x)),
                 has.na = has.na), class =  "density")
}

integr <- function(x, f)
{
       
       # var is numeric
       if (!is.numeric(x))
       {
              stop('The variable of integration "x" is not numeric.')
       }

       # integrand is numeric
       if (!is.numeric(f))
       {
              stop('The integrand "f" is not numeric.')
       }

       # length(var)=length(int)
       if (length(x) != length(f))
       {
              stop('The lengths of the variable of integration and the integrand do not match.')
       }

      # get lengths of var and integrand
       n = length(x)

       # trapezoidal integration
       integral = 0.5*sum((x[2:n] - x[1:(n-1)]) * (f[2:n] + f[1:(n-1)]))

       # print definite integral
       return(integral)
}


##compute JSD distance matrix
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
	KLD <- function(x,y) sum(x*log(x/y))
	JSD <- function(x,y) sqrt(0.5*KLD(x,(x+y)/2)+0.5*KLD(y,(x+y)/2))
	matrixColSize <- length(colnames(inMatrix))
	matrixRowSize <- length(rownames(inMatrix))
	colnames <- colnames(inMatrix)
	resultsMatrix <- matrix(0, matrixColSize, matrixColSize) 
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
	for(i in 1:matrixColSize) {
		for(j in 1:matrixColSize) { 
			resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
			as.vector(inMatrix[,j]))
		}
	}
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix) 
 }

JSDist <- function(x,y) sqrt(dist.JSD(x,y))


##find number of clusters, cluster data
pam.clustering=function(x,k) { # x=dist matrix, k= No. clusters
require(cluster)
      cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
            return(cluster)
                        }


#simulate random matrices for non-tree JS distances
makeMat = function(x) {
    set.seed(2718)
    mat = matrix(ncol=x, nrow=x)
    for (i in 1:x) {
        mat[, i] = rnorm(x)
    }
    return(mat)
}

##kurtosis
kurtosis.sub <-
    function (x, na.rm = FALSE, method = c("excess","moment", "fisher"), ...)
{
    
    method = match.arg(method)

    stopifnot(NCOL(x) == 1)

    # Warnings:
    if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
        warning("argument is not numeric or logical: returning NA")
        return(as.numeric(NA))}

    # Remove NAs:
    if (na.rm) x = x[!is.na(x)]

    # Kurtosis:
    n = length(x)
    if (is.integer(x)) x = as.numeric(x)
    if (method == "moment") {
        kurtosis = sum((x-mean(x))^4/as.numeric(var(x))^2)/length(x)
    }
     if (method == "excess") {
        kurtosis = sum((x-mean(x))^4/var(x)^2)/length(x) - 3
    }

    if (method == "fisher") {
        kurtosis = ((n+1)*(n-1)*((sum(x^4)/n)/(sum(x^2)/n)^2 -
            (3*(n-1))/(n+1)))/((n-2)*(n-3))
    }

    # Return Value:
    kurtosis
}


#plot gaussian
#gaussX   <- seq(5,15,length=1000)
#gaussY   <- dnorm(gaussX,mean=10, sd=3)
#gau <- plot(gaussX,gaussY, type="l", lwd=1)

`%ni%` <- Negate(`%in%`)
drop.tip.ni <- function(phy, tip, trim.internal = TRUE, subtree = FALSE,
         root.edge = 0, rooted = is.rooted(phy)) {
  Ntip <- length(phy$tip.label)
  if (is.character(tip)) 
    tip <- which(phy$tip.label %ni% tip)

  phy <- reorder(phy)
  NEWROOT <- ROOT <- Ntip + 1
  Nnode <- phy$Nnode
  Nedge <- nrow(phy$edge)

  wbl <- !is.null(phy$edge.length)
  edge1 <- phy$edge[, 1]
  edge2 <- phy$edge[, 2]
  keep <- !(edge2 %in% tip)  

  ints <- edge2 > Ntip
  repeat {
    sel <- !(edge2 %in% edge1[keep]) & ints & keep
    if (!sum(sel)) 
      break
    keep[sel] <- FALSE
  }

  phy2 <- phy
  phy2$edge <- phy2$edge[keep, ]
  if (wbl) 
    phy2$edge.length <- phy2$edge.length[keep]
  TERMS <- !(phy2$edge[, 2] %in% phy2$edge[, 1])
  oldNo.ofNewTips <- phy2$edge[TERMS, 2]
  n <- length(oldNo.ofNewTips)
  idx.old <- phy2$edge[TERMS, 2]
  phy2$edge[TERMS, 2] <- rank(phy2$edge[TERMS, 2])
  phy2$tip.label <- phy2$tip.label[-tip]
  if (!is.null(phy2$node.label))
    phy2$node.label <-
      phy2$node.label[sort(unique(phy2$edge[, 1])) - Ntip]
  phy2$Nnode <- nrow(phy2$edge) - n + 1L
  i <- phy2$edge > n
  phy2$edge[i] <- match(phy2$edge[i], sort(unique(phy2$edge[i]))) + n
  storage.mode(phy2$edge) <- "integer"
  collapse.singles(phy2)
}

cbind.fill<-function(...){
    nm <- list(...) 
    nm <-lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

#colless
collessit <- function(x, metric = "colless") {
    if (metric == "colless") {
        xx <- as.treeshape(x)  # convert to apTreeshape format
        colless(xx, "yule")  # calculate colless' metric
    } else if (metric == "gamma") {
        gammaStat(x)
    } else stop("metric should be one of colless or gamma")
}

#kmeansBIC
kmeansBIC = function(fit){

m = ncol(fit$centers)
n = length(fit$cluster)
k = nrow(fit$centers)
D = fit$tot.withinss
return(data.frame(BIC = D + log(n)*m*k))
}


jb <- function(x,y){
	chol_x = chol(x)
	chol_y = chol(y)
	chol_xy = chol((x+y)/2)
		jb_div = log(prod(diag(chol_xy)^0.2)) - 
			log(prod(diag(chol_x)^0.2)*prod(diag(chol_y)^0.2))/2
return(jb_div)			
	}

tree.drop.tip<-function(tree,index){
    dat<-names(table(index))
    name.phylo<-tree$tip.label
    nb<-length(dat)
    nbo<-length(name.phylo)
    if(nbo<=nb){cat("Error! You want to remove more species than there is in the tree")
        break;
    }else{
        ind<-NULL
        for(i in 1:nbo){
            if(any(name.phylo[i]==dat)!=TRUE){ind[i]<-name.phylo[i]}
        }
        #cat(ind[!is.na(ind)])
        tree<-drop.tip(tree,ind)
        return(tree)
    }}



sactter.grid <- function (x, y = NULL, z = NULL, color = par("col"), pch = NULL, 
					main = NULL, sub = NULL, xlim = NULL, ylim = NULL, zlim = NULL, 
					xlab = NULL, ylab = NULL, zlab = NULL, scale.y = 1, angle = 40, 
					axis = TRUE, tick.marks = TRUE, label.tick.marks = TRUE, 
					x.ticklabs = NULL, y.ticklabs = NULL, z.ticklabs = NULL, 
					y.margin.add = 0, grid = TRUE, box = TRUE, lab = par("lab"), 
					lab.z = mean(lab[1:2]), type = "p", highlight.3d = FALSE, 
					mar = c(5, 3, 4, 3) + 0.1, bg = par("bg"), col.axis = par("col.axis"), 
					col.grid = "grey", col.lab = par("col.lab"), cex.symbols = par("cex"), 
					cex.axis = 0.8 * par("cex.axis"), cex.lab = par("cex.lab"), 
					font.axis = par("font.axis"), font.lab = par("font.lab"), 
					lty.axis = par("lty"), lty.grid = par("lty"), lty.hide = NULL, 
					lty.hplot = par("lty"), log = "", ...) 
{
	mem.par <- par(mar = mar)
	x.scal <- y.scal <- z.scal <- 1
	xlabel <- if (!missing(x)) 
		deparse(substitute(x))
	ylabel <- if (!missing(y)) 
		deparse(substitute(y))
	zlabel <- if (!missing(z)) 
		deparse(substitute(z))
	if (highlight.3d && !missing(color)) 
		warning("color is ignored when highlight.3d = TRUE")
	if (!is.null(d <- dim(x)) && (length(d) == 2) && (d[2] >= 
																											4)) 
		color <- x[, 4]
	else if (is.list(x) && !is.null(x$color)) 
		color <- x$color
	xyz <- xyz.coords(x = x, y = y, z = z, xlab = xlabel, ylab = ylabel, 
										zlab = zlabel, log = log)
	if (is.null(xlab)) {
		xlab <- xyz$xlab
		if (is.null(xlab)) 
			xlab <- ""
	}
	if (is.null(ylab)) {
		ylab <- xyz$ylab
		if (is.null(ylab)) 
			ylab <- ""
	}
	if (is.null(zlab)) {
		zlab <- xyz$zlab
		if (is.null(zlab)) 
			zlab <- ""
	}
	if (length(color) == 1) 
		color <- rep(color, length(xyz$x))
	else if (length(color) != length(xyz$x)) 
		stop("length(color) ", "must be equal length(x) or 1")
	angle <- (angle%%360)/90
	yz.f <- scale.y * abs(if (angle < 1) angle else if (angle > 
																												3) angle - 4 else 2 - angle)
	yx.f <- scale.y * (if (angle < 2) 
		1 - angle
										 else angle - 3)
	if (angle > 2) {
		temp <- xyz$x
		xyz$x <- xyz$y
		xyz$y <- temp
		temp <- xlab
		xlab <- ylab
		ylab <- temp
		temp <- xlim
		xlim <- ylim
		ylim <- temp
	}
	angle.1 <- (1 < angle && angle < 2) || angle > 3
	angle.2 <- 1 <= angle && angle <= 3
	dat <- cbind(as.data.frame(xyz[c("x", "y", "z")]), col = color)
	if (!is.null(xlim)) {
		xlim <- range(xlim)
		dat <- dat[xlim[1] <= dat$x & dat$x <= xlim[2], , drop = FALSE]
	}
	if (!is.null(ylim)) {
		ylim <- range(ylim)
		dat <- dat[ylim[1] <= dat$y & dat$y <= ylim[2], , drop = FALSE]
	}
	if (!is.null(zlim)) {
		zlim <- range(zlim)
		dat <- dat[zlim[1] <= dat$z & dat$z <= zlim[2], , drop = FALSE]
	}
	n <- nrow(dat)
	if (n < 1) 
		stop("no data left within (x|y|z)lim")
	y.range <- range(dat$y[is.finite(dat$y)])
	if (type == "p" || type == "h") {
		y.ord <- rev(order(dat$y))
		dat <- dat[y.ord, ]
		if (length(pch) > 1) 
			if (length(pch) != length(y.ord)) 
				stop("length(pch) ", "must be equal length(x) or 1")
		else pch <- pch[y.ord]
		if (length(bg) > 1) 
			if (length(bg) != length(y.ord)) 
				stop("length(bg) ", "must be equal length(x) or 1")
		else bg <- bg[y.ord]
		if (length(cex.symbols) > 1) 
			if (length(cex.symbols) != length(y.ord)) 
				stop("length(cex.symbols) ", "must be equal length(x) or 1")
		else cex.symbols <- cex.symbols[y.ord]
		daty <- dat$y
		daty[!is.finite(daty)] <- mean(daty[is.finite(daty)])
		if (highlight.3d && !(all(diff(daty) == 0))) 
			dat$col <- rgb(red = seq(0, 1, length = n) * (y.range[2] - 
																											daty)/diff(y.range), green = 0, blue = 0)
	}
	p.lab <- par("lab")
	y.range <- range(dat$y[is.finite(dat$y)], ylim)
	y.prty <- pretty(y.range, n = lab[2], min.n = max(1, min(0.5 * 
																													 	lab[2], p.lab[2])))
	y.scal <- round(diff(y.prty[1:2]), digits = 12)
	y.add <- min(y.prty)
	dat$y <- (dat$y - y.add)/y.scal
	y.max <- (max(y.prty) - y.add)/y.scal
	if (!is.null(ylim)) 
		y.max <- max(y.max, ceiling((ylim[2] - y.add)/y.scal))
	x.range <- range(dat$x[is.finite(dat$x)], xlim)
	x.prty <- pretty(x.range, n = lab[1], min.n = max(1, min(0.5 * 
																													 	lab[1], p.lab[1])))
	x.scal <- round(diff(x.prty[1:2]), digits = 12)
	dat$x <- dat$x/x.scal
	x.range <- range(x.prty)/x.scal
	x.max <- ceiling(x.range[2])
	x.min <- floor(x.range[1])
	if (!is.null(xlim)) {
		x.max <- max(x.max, ceiling(xlim[2]/x.scal))
		x.min <- min(x.min, floor(xlim[1]/x.scal))
	}
	x.range <- range(x.min, x.max)
	z.range <- range(dat$z[is.finite(dat$z)], zlim)
	z.prty <- pretty(z.range, n = lab.z, min.n = max(1, min(0.5 * 
																														lab.z, p.lab[2])))
	z.scal <- round(diff(z.prty[1:2]), digits = 12)
	dat$z <- dat$z/z.scal
	z.range <- range(z.prty)/z.scal
	z.max <- ceiling(z.range[2])
	z.min <- floor(z.range[1])
	if (!is.null(zlim)) {
		z.max <- max(z.max, ceiling(zlim[2]/z.scal))
		z.min <- min(z.min, floor(zlim[1]/z.scal))
	}
	z.range <- range(z.min, z.max)
	plot.new()
	if (angle.2) {
		x1 <- x.min + yx.f * y.max
		x2 <- x.max
	}
	else {
		x1 <- x.min
		x2 <- x.max + yx.f * y.max
	}
	plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
	temp <- strwidth(format(rev(y.prty))[1], cex = cex.axis/par("cex"))
	if (angle.2) 
		x1 <- x1 - temp - y.margin.add
	else x2 <- x2 + temp + y.margin.add
	plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
	if (angle > 2) 
		par(usr = par("usr")[c(2, 1, 3:4)])
	usr <- par("usr")
	title(main, sub, ...)
	if ("xy" %in% grid || grid) {
		i <- x.min:x.max
		segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + 
						 	z.min, col = col.grid, lty = lty.grid)
		i <- 0:y.max
		segments(x.min + (i * yx.f), i * yz.f + z.min, x.max + 
						 	(i * yx.f), i * yz.f + z.min, col = col.grid, lty = lty.grid)
	}
	if ("xz" %in% grid) {
		i <- x.min:x.max
		segments(i + (yx.f * y.max), yz.f * y.max + z.min, 
						 i + (yx.f * y.max), yz.f * y.max + z.max, 
						 col = col.grid, lty = lty.grid)
		temp <- yx.f * y.max
		temp1 <- yz.f * y.max
		i <- z.min:z.max
		segments(x.min + temp,temp1 + i, 
						 x.max + temp,temp1 + i , col = col.grid, lty = lty.grid)

	}
	
	if ("yz" %in% grid) {
		i <- 0:y.max
		segments(x.min + (i * yx.f), i * yz.f + z.min,  
						 x.min + (i * yx.f) ,i * yz.f + z.max,  
						 col = col.grid, lty = lty.grid)
		temp <- yx.f * y.max
		temp1 <- yz.f * y.max
		i <- z.min:z.max
		segments(x.min + temp,temp1 + i, 
						 x.min, i , col = col.grid, lty = lty.grid)
		

	
	}
	
	if (axis) {
		xx <- if (angle.2) 
			c(x.min, x.max)
		else c(x.max, x.min)
		if (tick.marks) {
			xtl <- (z.max - z.min) * (tcl <- -par("tcl"))/50
			ztl <- (x.max - x.min) * tcl/50
			mysegs <- function(x0, y0, x1, y1) segments(x0, 
																									y0, x1, y1, col = col.axis, lty = lty.axis)
			i.y <- 0:y.max
			mysegs(yx.f * i.y - ztl + xx[1], yz.f * i.y + z.min, 
						 yx.f * i.y + ztl + xx[1], yz.f * i.y + z.min)
			i.x <- x.min:x.max
			mysegs(i.x, -xtl + z.min, i.x, xtl + z.min)
			i.z <- z.min:z.max
			mysegs(-ztl + xx[2], i.z, ztl + xx[2], i.z)
			if (label.tick.marks) {
				las <- par("las")
				mytext <- function(labels, side, at, ...) mtext(text = labels, 
																												side = side, at = at, line = -0.5, col = col.lab, 
																												cex = cex.axis, font = font.lab, ...)
				if (is.null(x.ticklabs)) 
					x.ticklabs <- format(i.x * x.scal)
				mytext(x.ticklabs, side = 1, at = i.x)
				if (is.null(z.ticklabs)) 
					z.ticklabs <- format(i.z * z.scal)
				mytext(z.ticklabs, side = if (angle.1) 
					4
							 else 2, at = i.z, adj = if (0 < las && las < 
							 															3) 
							 	1
							 else NA)
				temp <- if (angle > 2) 
					rev(i.y)
				else i.y
				if (is.null(y.ticklabs)) 
					y.ticklabs <- format(y.prty)
				else if (angle > 2) 
					y.ticklabs <- rev(y.ticklabs)
				text(i.y * yx.f + xx[1], i.y * yz.f + z.min, 
						 y.ticklabs, pos = if (angle.1) 
						 	2
						 else 4, offset = 1, col = col.lab, cex = cex.axis/par("cex"), 
						 font = font.lab)
			}
		}
		mytext2 <- function(lab, side, line, at) mtext(lab, 
																									 side = side, line = line, at = at, col = col.lab, 
																									 cex = cex.lab, font = font.axis, las = 0)
		lines(c(x.min, x.max), c(z.min, z.min), col = col.axis, 
					lty = lty.axis)
		mytext2(xlab, 1, line = 1.5, at = mean(x.range))
		lines(xx[1] + c(0, y.max * yx.f), c(z.min, y.max * yz.f + 
																					z.min), col = col.axis, lty = lty.axis)
		mytext2(ylab, if (angle.1) 
			2
						else 4, line = 0.5, at = z.min + y.max * yz.f)
		lines(xx[c(2, 2)], c(z.min, z.max), col = col.axis, 
					lty = lty.axis)
		mytext2(zlab, if (angle.1) 
			4
						else 2, line = 1.5, at = mean(z.range))
		if (box) {
			if (is.null(lty.hide)) 
				lty.hide <- lty.axis
			temp <- yx.f * y.max
			temp1 <- yz.f * y.max
			lines(c(x.min + temp, x.max + temp), c(z.min + temp1, 
																						 z.min + temp1), col = col.axis, lty = lty.hide)
			lines(c(x.min + temp, x.max + temp), c(temp1 + z.max, 
																						 temp1 + z.max), col = col.axis, lty = lty.axis)
			temp <- c(0, y.max * yx.f)
			temp1 <- c(0, y.max * yz.f)
			lines(temp + xx[2], temp1 + z.min, col = col.axis, 
						lty = lty.hide)
			lines(temp + x.min, temp1 + z.max, col = col.axis, 
						lty = lty.axis)
			temp <- yx.f * y.max
			temp1 <- yz.f * y.max
			lines(c(temp + x.min, temp + x.min), c(z.min + temp1, 
																						 z.max + temp1), col = col.axis, lty = if (!angle.2) 
																						 	lty.hide
						else lty.axis)
			lines(c(x.max + temp, x.max + temp), c(z.min + temp1, 
																						 z.max + temp1), col = col.axis, lty = if (angle.2) 
																						 	lty.hide
						else lty.axis)
		}
	}
	x <- dat$x + (dat$y * yx.f)
	z <- dat$z + (dat$y * yz.f)
	col <- as.character(dat$col)
	if (type == "h") {
		z2 <- dat$y * yz.f + z.min
		segments(x, z, x, z2, col = col, cex = cex.symbols, 
						 lty = lty.hplot, ...)
		points(x, z, type = "p", col = col, pch = pch, bg = bg, 
					 cex = cex.symbols, ...)
	}
	else points(x, z, type = type, col = col, pch = pch, bg = bg, 
							cex = cex.symbols, ...)
	if (axis && box) {
		lines(c(x.min, x.max), c(z.max, z.max), col = col.axis, 
					lty = lty.axis)
		lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + 
						z.max, col = col.axis, lty = lty.axis)
		lines(xx[c(1, 1)], c(z.min, z.max), col = col.axis, 
					lty = lty.axis)
	}
	ob <- ls()
	rm(list = ob[!ob %in% c("angle", "mar", "usr", "x.scal", 
													"y.scal", "z.scal", "yx.f", "yz.f", "y.add", "z.min", 
													"z.max", "x.min", "x.max", "y.max", "x.prty", "y.prty", 
													"z.prty")])
	rm(ob)
	invisible(list(xyz.convert = function(x, y = NULL, z = NULL) {
		xyz <- xyz.coords(x, y, z)
		if (angle > 2) {
			temp <- xyz$x
			xyz$x <- xyz$y
			xyz$y <- temp
		}
		y <- (xyz$y - y.add)/y.scal
		return(list(x = xyz$x/x.scal + yx.f * y, y = xyz$z/z.scal + 
									yz.f * y))
	}, points3d = function(x, y = NULL, z = NULL, type = "p", 
												 ...) {
		xyz <- xyz.coords(x, y, z)
		if (angle > 2) {
			temp <- xyz$x
			xyz$x <- xyz$y
			xyz$y <- temp
		}
		y2 <- (xyz$y - y.add)/y.scal
		x <- xyz$x/x.scal + yx.f * y2
		y <- xyz$z/z.scal + yz.f * y2
		mem.par <- par(mar = mar, usr = usr)
		on.exit(par(mem.par))
		if (type == "h") {
			y2 <- z.min + yz.f * y2
			segments(x, y, x, y2, ...)
			points(x, y, type = "p", ...)
		} else points(x, y, type = type, ...)
	}, plane3d = function(Intercept, x.coef = NULL, y.coef = NULL, 
												lty = "dashed", lty.box = NULL, ...) {
		if (!is.atomic(Intercept) && !is.null(coef(Intercept))) Intercept <- coef(Intercept)
		if (is.null(lty.box)) lty.box <- lty
		if (is.null(x.coef) && length(Intercept) == 3) {
			x.coef <- Intercept[if (angle > 2) 3 else 2]
			y.coef <- Intercept[if (angle > 2) 2 else 3]
			Intercept <- Intercept[1]
		}
		mem.par <- par(mar = mar, usr = usr)
		on.exit(par(mem.par))
		x <- x.min:x.max
		ltya <- c(lty.box, rep(lty, length(x) - 2), lty.box)
		x.coef <- x.coef * x.scal
		z1 <- (Intercept + x * x.coef + y.add * y.coef)/z.scal
		z2 <- (Intercept + x * x.coef + (y.max * y.scal + y.add) * 
					 	y.coef)/z.scal
		segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, 
						 lty = ltya, ...)
		y <- 0:y.max
		ltya <- c(lty.box, rep(lty, length(y) - 2), lty.box)
		y.coef <- (y * y.scal + y.add) * y.coef
		z1 <- (Intercept + x.min * x.coef + y.coef)/z.scal
		z2 <- (Intercept + x.max * x.coef + y.coef)/z.scal
		segments(x.min + y * yx.f, z1 + y * yz.f, x.max + y * 
						 	yx.f, z2 + y * yz.f, lty = ltya, ...)
	}, box3d = function(...) {
		mem.par <- par(mar = mar, usr = usr)
		on.exit(par(mem.par))
		lines(c(x.min, x.max), c(z.max, z.max), ...)
		lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + 
						z.max, ...)
		lines(c(0, y.max * yx.f) + x.min, c(0, y.max * yz.f) + 
						z.max, ...)
		lines(c(x.max, x.max), c(z.min, z.max), ...)
		lines(c(x.min, x.min), c(z.min, z.max), ...)
		lines(c(x.min, x.max), c(z.min, z.min), ...)
	}))
}

"clade.members" <-
function(x, phy, tip.labels=FALSE, include.nodes=FALSE){
    
    # NEW2OLD: CONVERTED...
	
	# returns a vector of the tips that descend from an identified node
	if(class(phy) != "phylo") stop("Phylogeny required")
	
	NallNodes <- max(phy$edge)
	Ntips <- max(phy$edge) - phy$Nnode
	
	if(!(x %in% 1:NallNodes)) stop("Node not in range for phylogeny")
	
	# find the children of the node, append them to the vector of nodes (x)
	# and remove the parent, until all the nodes in the vector are tips...
	# now updated to keep track of parents...
	
	intN <- x[x>Ntips]
	descNode <- numeric(length=0)
	
	while(length(intN) > 0){
		
		minIntN <- min(intN)
		childOfMinIntN <- with(phy, edge[,2][which(edge[,1] == minIntN)])
		
		descNode <- c(descNode, minIntN)
		x <- c(x[x != minIntN], childOfMinIntN)
		
		intN <- x[x>Ntips]
	}
	
	RET <- unique(x)
	
	if(tip.labels) {
		RET <- phy$tip.label[x]
	} 
	
	if(include.nodes) {
		RET <- list(tips=RET, nodes=descNode)
	}
	
	return(RET)
	
}

"clade.members.list" <-
function(phy, tips=FALSE, tip.labels=FALSE, include.nodes=FALSE){
    
    # OLD2NEW CONVERTED
    
	# returns a list of vectors showing the tips
	# subtending from each node in the tree
	if(class(phy) != "phylo") stop("Phylogeny required")

	nodes <- 1:max(phy$edge)
	
	if(!tips) nodes <- nodes[nodes > length(nodes) - phy$Nnode]
	
	clade.list <- mapply(clade.members, nodes, MoreArgs=list(phy=phy, tip.labels=tip.labels, include.nodes=include.nodes), SIMPLIFY=FALSE)
	names(clade.list) <- nodes
	
	return(clade.list)
}
kurtosis.sub <- function(x, na.rm = FALSE, method = c("moment"), 
        ...) {
        method = match.arg(method)
        stopifnot(NCOL(x) == 1)
        if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
            warning("argument is not numeric or logical: returning NA")
            return(as.numeric(NA))
        }
        if (na.rm) 
            x = x[!is.na(x)]
        n = length(x)
        if (is.integer(x)) 
            x = as.numeric(x)
        if (method == "moment") {
            kurtosis = sum((x - mean(x))^4/as.numeric(var(x))^2)/length(x)
        }
        if (method == "excess") {
            kurtosis = sum((x - mean(x))^4/var(x)^2)/length(x) - 
                3
        }
        if (method == "fisher") {
            kurtosis = ((n + 1) * (n - 1) * ((sum(x^4)/n)/(sum(x^2)/n)^2 - 
                (3 * (n - 1))/(n + 1)))/((n - 2) * (n - 3))
        }
        kurtosis
    }
    skewness <- function(x, na.rm = FALSE) {
        if (is.matrix(x)) 
            apply(x, 2, skewness, na.rm = na.rm)
        else if (is.vector(x)) {
            if (na.rm) 
                x <- x[!is.na(x)]
            n <- length(x)
            (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
        }
        else if (is.data.frame(x)) 
            sapply(x, skewness, na.rm = na.rm)
        else skewness(as.vector(x), na.rm = na.rm)
    }
    peak_height <- function(x) {
        kernelG <- function(x, mean = 0, sd = 1) dnorm(x, mean = mean, 
            sd = sd)
        dens <- function(x, bw = bw.nrd0, kernel = kernelG, n = 4096, 
            from = min(x) - 3 * sd, to = max(x) + 3 * sd, adjust = 1, 
            ...) {
            if (has.na <- any(is.na(x))) {
                x <- na.omit(x)
                if (length(x) == 0) 
                  stop("no finite or non-missing data!")
            }
            sd <- (if (is.numeric(bw)) 
                bw[1]
            else bw(x)) * adjust
            X <- seq(from, to, len = n)
            M <- outer(X, x, kernel, sd = sd, ...)
            structure(list(x = X, y = rowMeans(M), bw = sd, call = match.call(), 
                n = length(x), data.name = deparse(substitute(x)), 
                has.na = has.na), class = "density")
        }
        integr <- function(x, f) {
            if (!is.numeric(x)) {
                stop("The variable of integration \"x\" is not numeric.")
            }
            if (!is.numeric(f)) {
                stop("The integrand \"f\" is not numeric.")
            }
            if (length(x) != length(f)) {
                stop("The lengths of the variable of integration and the integrand do not match.")
            }
            n = length(x)
            integral = 0.5 * sum((x[2:n] - x[1:(n - 1)]) * (f[2:n] + 
                f[1:(n - 1)]))
            return(integral)
        }
        d <- dens(log(x))
        dsc <- max(d$y/integr(d$x, d$y))
        return(dsc)
    }
getE<-function(phylo){
e = eigen(graph.laplacian(graph.adjacency(data.matrix(dist.nodes(phylo)), 
            weighted = T), normalized = F), symmetric = T, only.values = T)
	x = subset(e$values, e$values >= 0)
	d = dens(log(x))
	dsc = d$y/(integr(d$x,d$y))
	principal_eigenvalue <- max(x)
	skewness <- skewness(x)
    peak_height <- peak_height(x)
	gaps<-abs(diff(e$values))
	extra_col<-c(1:length(gaps))
	gapMat<-as.matrix(gaps)
	gapMatCol<-cbind(extra_col,gapMat)
	gapMatMax<-subset(gapMatCol,gapMatCol[,2]==max(gapMatCol[,2]))
	res<-list(d$x,dsc,eigengap=gapMatMax,principal_eigenvalue=principal_eigenvalue, 
            asymmetry=skewness, peakedness=peak_height)   
	return(res)
}

abtree<-function (ttree, sliceTime, drop.extinct = TRUE, plot = FALSE) 
{
    if (!inherits(ttree, "phylo")) {
        stop("ttree is not of class phylo")
    }
    tslice <- max(nodeHeights(ttree)) - sliceTime
    dnode <- node.depth.edgelength(ttree)
    cedge <- which((dnode[ttree$edge[, 1]] < tslice) & (dnode[ttree$edge[, 
        2]] >= tslice))
    droppers <- numeric()
    propPartTree <- prop.part(ttree)
    for (i in 1:length(cedge)) {
        desc <- ttree$edge[cedge[i], 2]
        if (desc > Ntip(ttree)) {
            desctip <- propPartTree[[desc - Ntip(ttree)]]
            droppers <- c(droppers, desctip[-1])
        }
    }
    stree <- drop.tip(ttree, droppers)
    dnode <- node.depth.edgelength(stree)
    cedge <- (dnode[stree$edge[, 2]] >= tslice)
    cnode_depth <- dnode[stree$edge[cedge, 1]]
    stree$edge.length[cedge] <- tslice - cnode_depth
    stree$root.time <- max(nodeHeights(ttree))
   stree1 <- stree
    return(stree1)
}
rootedge.to.singleton<-function(tree){
    cw<-reorder(tree,"cladewise")
    root.edge<-if(!is.null(cw$root.edge)) cw$root.edge else 0
    cw$edge[which(cw$edge>Ntip(cw))]<-cw$edge[which(cw$edge>Ntip(cw))]+1
    cw$edge<-rbind(Ntip(cw)+c(1,2),cw$edge)
    cw$Nnode<-cw$Nnode+1
    cw$edge.length<-c(root.edge,cw$edge.length)
    cw
}

##set colors
red=rgb(255, 0, 0, 100, names = NULL, maxColorValue = 255)
green=rgb(0,255, 0, 100, names = NULL, maxColorValue = 255)
blue=rgb(0, 0, 255, 100, names = NULL, maxColorValue = 255)
yellow=rgb(255,255, 0, 100, names = NULL, maxColorValue = 255)
orange=rgb(255,150, 0, 100, names = NULL, maxColorValue = 255)
cyan=rgb(50, 255, 255, 100, names = NULL, maxColorValue = 255)
purple=rgb(255, 0, 255, 100, names = NULL, maxColorValue = 255) 
gray=rgb(100,100, 200, 100, names = NULL, maxColorValue = 255) 
