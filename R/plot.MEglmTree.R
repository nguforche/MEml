#' Plot method for MEglmTree object
#' Take a MEglmTree and plot the barplot of fitted classes or boxplot of 
#' fitted probability 
#
#' @name Plot.MEglmTree
#
#' @param x  fited MEglmTree or MECtree object.
#' @param tp_args,ip_args,xlim,ylim see partykit
#' @param main,ep_args see partykid 
#' @param terminal_panel see partykit 
#' @param inner_panel see partykit 
#' @param drop_terminal see partykit 
#' @param debug see partykit 
#' @param edge_panel see partykit
#' @param node see partykit
#' @param tnex see partykit
#' @param nx,ny see partykit 
#' @param newpage,pop,gp,tree,low.is.green see partykit  
#' @param \dots Further arguments passed to partykit
# 
#' @author  Che Ngufor Ngufor.Che@@mayo.edu
#' @import  partykit Matrix vcd 
# 
#' @importFrom grid viewport gpar grid.layout
#' @importFrom grid pushViewport grid.rect grid.lines
#' @importFrom grid unit grid.text popViewport upViewport
#' @importFrom grid grid.points grid.yaxis grid.xaxis 
#' @importFrom grid grid.newpage plotViewport
#
#
NULL 
#
#' @rdname  Plot.MEglmTree
#' @method plot MEglmTree
#' @export
plot.MEglmTree <- function(x, main = NULL,
                       terminal_panel = node_terminal, tp_args = list(),
		       inner_panel = node_inner, ip_args = list(),
                       edge_panel = edge_simple, ep_args = list(),
		       drop_terminal = FALSE, tnex = 1, 
		       newpage = TRUE, pop = TRUE, gp = gpar(), ...)
{


if (!inherits(x, "MEglmTree")) stop("Object must be a \"MEglmTree \"'")
    x <- x$tree.fit 
    ### extract tree
    node <- node_party(x)
    ### total number of terminal nodes
    nx <- width(node)
    ### maximal depth of the tree
    ny <- depth(node, root = TRUE)

    ## setup newpage
    if (newpage) grid.newpage()

    ## setup root viewport
    root_vp <- viewport(layout = grid.layout(3, 3, 
    			heights = unit(c(ifelse(is.null(main), 0, 3), 1, 1), 
                                      c("lines", "null", "lines")),
    			widths = unit(c(1, 1, 1), 
                                     c("lines", "null", "lines"))), 
    			name = "root",
			gp = gp)       
    pushViewport(root_vp)
  
    ## viewport for main title (if any)
    if (!is.null(main)) {
        main_vp <- viewport(layout.pos.col = 2, layout.pos.row = 1, 
                            name = "main")
        pushViewport(main_vp)
        grid.text(y = unit(1, "lines"), main, just = "center")
        upViewport()
    }

    ## setup viewport for tree
    tree_vp <- viewport(layout.pos.col = 2, layout.pos.row = 2, 
    			xscale = c(0, nx), yscale = c(0, ny + (tnex - 1)), 
                        name = "tree")
    pushViewport(tree_vp)

    ### setup panel functions (if necessary)
    if(inherits(terminal_panel, "grapcon_generator"))
      terminal_panel <- do.call("terminal_panel", c(list(x), as.list(tp_args)))
    if(inherits(inner_panel, "grapcon_generator"))
      inner_panel <- do.call("inner_panel", c(list(x), as.list(ip_args)))
    if(inherits(edge_panel, "grapcon_generator"))
      edge_panel <- do.call("edge_panel", c(list(x), as.list(ep_args)))


    if((nx <= 1 & ny <= 1)) {
      pushViewport(plotViewport(margins = rep(1.5, 4), name = paste("Node", id_node(node), sep = "")))
      terminal_panel(node)
    } else {
      ## call the workhorse
      .plot_node(node,
        xlim = c(0, nx), ylim = c(0, ny - 0.5 + (tnex - 1)),
        nx = nx, ny = ny, 
        terminal_panel = terminal_panel,
        inner_panel = inner_panel,
        edge_panel = edge_panel,
        tnex = tnex,
        drop_terminal = drop_terminal,
        debug = FALSE)
    }
    upViewport()
    if (pop) popViewport() else upViewport()
}
#' @rdname  Plot.MEglmTree
#' @method plot MECTree
#' @export
plot.MECTree <- function(x,...){
	if(x$con.tree) {
	   class(x) <- "MEglmTree"
	   plot(x, ...)
	} else {
	   heat.tree(x$tree.fit, type=1, varlen=0, faclen=0,  
			branch=1,trace=1,under=FALSE,split.cex= 1.2,
			fallen.leaves=FALSE, extra=101,boxes.include.gap=FALSE, 
			shadow.col = "gray90")
}			
}
#
#' @rdname  Plot.MEglmTree
#' @method print MECTree
#' @export
print.MECTree <- function(x,...){
if (!inherits(x, "MECTree")) stop("Object must be a \"MECTree \"'")
    print("*** Mixed Effect Classification Trees ***")
    print(x$tree.fit)
    print("Estimated covariance matrix of fided effects:")
    print(vcov(x$glmer.fit))
    print("Estimated parameters with confidence intervals")
    print(x$glmer.CI)
    print(paste("Log likelihood: ", x$logLik))
}
#' @rdname  Plot.MEglmTree
#' @export
heat.tree <- function(tree, low.is.green=FALSE, ...) { # dots args passed to prp
y <- tree$frame$yval
#if(low.is.green)
#y <- -y
#max <- max(y)
#min <- min(y)
#cols <- rainbow(length(y), end=.36)[ifelse(y > y[1], (y-y[1]) * (99-50) / (max-y[1]) + 50,
#(y-min) * (50-1) / (y[1]-min) + 1)]
cols <- rainbow(length(y), end=.55, alpha = 0.8)[rank(y)]
prp(tree, branch.col=cols, box.col=cols, ...)
}
#
#' @rdname  Plot.MEglmTree
.plot_node <- function(node, xlim, ylim, nx, ny, 
               terminal_panel, inner_panel, edge_panel,
	       tnex = 2, drop_terminal = TRUE, debug = FALSE) {

    ### the workhorse for plotting trees
 
    ### set up viewport for terminal node
    if (is.terminal(node)) {
        x <- xlim[1] + diff(xlim)/2
        y <- ylim[1] + 0.5
       
        tn_vp <- viewport(x = unit(x, "native"),
                          y = unit(y, "native") - unit(0.5, "lines"),
                          width = unit(1, "native"), 
                          height = unit(tnex, "native") - unit(1, "lines"),
			  just = c("center", "top"),
                          name = paste("Node", id_node(node), sep = ""))
        pushViewport(tn_vp)
        if (debug)
            grid.rect(gp = gpar(lty = "dotted", col = 4))
        terminal_panel(node) 
        upViewport()
        return(NULL)
    }    

    ## convenience function for computing relative position of splitting node
    pos_frac <- function(node) {
      if(is.terminal(node)) 0.5 else {
        width_kids <- sapply(kids_node(node), width)
        nk <- length(width_kids)
        rval <- if(nk %% 2 == 0) sum(width_kids[1:(nk/2)]) else
	  mean(cumsum(width_kids)[nk/2 + c(-0.5, 0.5)])
	rval/sum(width_kids)
      }
    }

    ## extract information
    split <- split_node(node)
    kids <- kids_node(node)
    width_kids <- sapply(kids, width)
    nk <- length(width_kids)

    ### position of inner node
    x0 <- xlim[1] + pos_frac(node) * diff(xlim)
    y0 <- max(ylim)

    ### relative positions of kids
    xfrac <- sapply(kids, pos_frac)
    x1lim <- xlim[1] + cumsum(c(0, width_kids))/sum(width_kids) * diff(xlim)
    x1 <- x1lim[1:nk] + xfrac * diff(x1lim)
    if (!drop_terminal) {
        y1 <- rep(y0 - 1, nk)
    } else {
        y1 <- ifelse(sapply(kids, is.terminal), tnex - 0.5, y0 - 1)
    }

    ### draw edges
    for(i in 1:nk) grid.lines(x = unit(c(x0, x1[i]), "native"), y = unit(c(y0, y1[i]), "native"))

    ### create viewport for inner node
    in_vp <- viewport(x = unit(x0, "native"),
                      y = unit(y0, "native"),
                      width = unit(1, "native"),
                      height = unit(1, "native") - unit(1, "lines"), 
                      name = paste("Node", id_node(node), sep = ""))
    pushViewport(in_vp)
    if(debug) grid.rect(gp = gpar(lty = "dotted"))
    inner_panel(node)
    upViewport()

    ### position of labels
    y1max <- max(y1)
    ypos <- y0 - (y0 - y1max) * 0.5
    xpos <- x0 - (x0 - x1) * 0.5 * (y0 - y1max)/(y0 - y1)

    ### setup labels
    for(i in 1:nk) {
      sp_vp <- viewport(x = unit(xpos[i], "native"),
                        y = unit(ypos, "native"),
                        width = unit(diff(x1lim)[i], "native"),
                        height = unit(1, "lines"), 
                        name =  paste("edge", id_node(node), "-", i, sep = ""))
      pushViewport(sp_vp)
      if(debug) grid.rect(gp = gpar(lty = "dotted", col = 2))
      edge_panel(node, i)
      upViewport()
    }

    ## call workhorse for kids
    for(i in 1:nk) .plot_node(kids[[i]],
      c(x1lim[i], x1lim[i+1]), c(y1[i], 1), nx, ny, 
      terminal_panel, inner_panel, edge_panel,
      tnex = tnex, drop_terminal = drop_terminal, debug = debug)
}

#
#' @rdname  Plot.MEglmTree
#' @export
NodeBoxplot <- function(obj, fitted = TRUE,  
                         col = "black",
		         fill = "lightgray",
		         width = 0.5,
		         yscale = NULL,
		         ylines = 3,
			     cex = 0.5,
		         id = TRUE,
                 mainlab = NULL, 
			 gp = gpar(), ...)
{
node <- NULL 

  if(fitted) Target = "fitted"
  else  Target = "Target"

   ## obtain dependent variable
  y <- obj$fitted.probs
  obj <- obj$tree.fit 
     
  mf <- model.frame(obj)   
  mf[, "nodes"] <- predict(obj, mf, type = "node") 
  mf[, "fitted"] <- y  
  stopifnot(is.numeric(y))
    if (is.null(yscale)) 
        yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
         
### panel function for boxplots in nodes
    rval <- function(node) {

## extract data
	nid <- id_node(node)
#	dat <- data_party(obj, nid)
#	dat <- subset(mf, nodes == nid)
	dat <- mf[mf$nodes == nid, ]	
	yn <- dat[[Target]] 
	wn <- NULL ##dat[["(weights)"]]
	if(is.null(wn)) 
    wn <- rep(1, length(yn))    
## parameter setup
	x <- boxplot(rep.int(yn, wn), plot = FALSE)
        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1), 
                                         c("lines", "null", "lines")),  
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_boxplot", nid, sep = ""),
			   gp = gp)

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
        if (is.null(mainlab)) {	
	  mainlab <- if(id) {
	    function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
  	  } else {
	    function(id, nobs) sprintf("n = %s", nobs)
	  }
        }
	if (is.function(mainlab)) {
          mainlab <- mainlab(names(obj)[nid], sum(wn))
	}
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col = 2, layout.pos.row = 2,
                         xscale = c(0, 1), yscale = yscale,
			 name = paste("node_boxplot", nid, "plot", 
                         sep = ""))

        pushViewport(plot)
	
	xl <- 0.5 - width/4
	xr <- 0.5 + width/4

        ## box & whiskers
        grid.lines(unit(c(xl, xr), "npc"), 
                   unit(x$stats[1], "native"), gp = gpar(col = col))
        grid.lines(unit(0.5, "npc"), 
                   unit(x$stats[1:2], "native"), gp = gpar(col = col, lty = 2))
        grid.rect(unit(0.5, "npc"), unit(x$stats[2], "native"), 
                  width = unit(width, "npc"), height = unit(diff(x$stats[c(2, 4)]), "native"),
                  just = c("center", "bottom"), 
                  gp = gpar(col = col, fill = fill))
        grid.lines(unit(c(0.5 - width/2, 0.5+width/2), "npc"), 
                   unit(x$stats[3], "native"), gp = gpar(col = col, lwd = 2))
        grid.lines(unit(0.5, "npc"), unit(x$stats[4:5], "native"), 
                   gp = gpar(col = col, lty = 2))
        grid.lines(unit(c(xl, xr), "npc"), unit(x$stats[5], "native"), 
                   gp = gpar(col = col))

        ## outlier
        n <- length(x$out)
        if (n > 0) {
            index <- 1:n ## which(x$out > yscale[1] & x$out < yscale[2])
            if (length(index) > 0)
                grid.points(unit(rep.int(0.5, length(index)), "npc"), 
                            unit(x$out[index], "native"),
                            size = unit(cex, "char"), gp = gpar(col = col))
        }
	
        grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
        upViewport(2)
    }
    
    return(rval)
}
class(NodeBoxplot) <- "grapcon_generator"
#
#' @rdname  Plot.MEglmTree
#' @export
NodeBarplot <- function(obj, fitted = TRUE, 
                 col = "black",
      		     fill = NULL,
			     beside = NULL,
		         ymax = NULL,
		         ylines = NULL,
		         widths = 1,
		         gap = NULL,
			     reverse = NULL,
		         id = TRUE,
                 mainlab = NULL,
			     gp = gpar(), ...)
{   
	if(fitted) Target = "fitted"
	else  Target = "Target"
 
    ## extract response
   ## obtain dependent variable
   y <- factor(obj$fitted.class)  
   stopifnot(is.factor(y) || isTRUE(all.equal(round(y), y)) || is.data.frame(y))

  levels(y) <- c("No", "Yes")   
  obj <- obj$tree.fit      
  mf <- model.frame(obj)   
  mf[, "nodes"] <- predict(obj, mf, type = "node") 
  mf[, "fitted"] <- y  
    
    ## FIXME: This could be avoided by
    ##   predict_party(obj, nodeids(obj, terminal = TRUE), type = "prob")
    ## but only for terminal nodes   
 
     probs_and_n <- function(x) {
    ## get id of this current node 
		nid <- id_node(x)
#	        y1 <- subset(mf, nodes == nid)[[Target]]
		y1 <- mf[mf$nodes == nid, ][[Target]]
#      
      if(!is.factor(y1)) {
        if(is.data.frame(y1)) {
	  y1 <- t(as.matrix(y1))
	} else {
          y1 <- factor(y1, levels = min(y):max(y))
	}
      }
      w <- NULL
      if(is.null(w)) w <- rep.int(1L, length(y1))
      sumw <- if(is.factor(y1)) tapply(w, y1, sum) else drop(y1 %*% w)
      sumw[is.na(sumw)] <- 0
      prob <- c(sumw/sum(w), sum(w))
      names(prob) <- c(if(is.factor(y1)) levels(y1) else rownames(y1), "nobs")
      prob
    }
    probs <- do.call("rbind", nodeapply(obj, nodeids(obj), probs_and_n, by_node = TRUE))
    nobs <- probs[, "nobs"]
    probs <- probs[, -ncol(probs), drop = FALSE]
    
    if(is.factor(y)) {
        ylevels <- levels(y)
	if(is.null(beside)) beside <- if(length(ylevels) < 3L) FALSE else TRUE
        if(is.null(ymax)) ymax <- if(beside) 1.1 else 1
	if(is.null(gap)) gap <- if(beside) 0.1 else 0
    } else {
        if(is.null(beside)) beside <- TRUE
        if(is.null(ymax)) ymax <- if(beside) max(probs) * 1.1 else max(probs)
        ylevels <- colnames(probs)
        if(length(ylevels) < 2) ylevels <- ""
	if(is.null(gap)) gap <- if(beside) 0.1 else 0
    }
    if(is.null(reverse)) reverse <- !beside
    if(is.null(fill)) fill <- gray.colors(length(ylevels))
    if(is.null(ylines)) ylines <- if(beside) c(3, 2) else c(1.5, 2.5)

    ### panel function for barplots in nodes
    rval <- function(node) {
    
        ## id
	nid <- id_node(node)
    
        ## parameter setup
        pred <- probs[nid,]
	if(reverse) {
	  pred <- rev(pred)
	  ylevels <- rev(ylevels)
	}
        np <- length(pred)
	nc <- if(beside) np else 1

	fill <- rep(fill, length.out = np)	
        widths <- rep(widths, length.out = nc)
	col <- rep(col, length.out = nc)
	ylines <- rep(ylines, length.out = 2)

	gap <- gap * sum(widths)
        yscale <- c(0, ymax)
        xscale <- c(0, sum(widths) + (nc+1)*gap)

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines[1], 1, ylines[2]), c("lines", "null", "lines")),
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_barplot", nid, sep = ""),
			   gp = gp)

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
        if (is.null(mainlab)) {	
	  mainlab <- if(id) {
	    function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
  	  } else {
	    function(id, nobs) sprintf("n = %s", nobs)
	  }
        }
	if (is.function(mainlab)) {
          mainlab <- mainlab(names(obj)[nid], nobs[nid])
	}
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                         xscale=xscale, yscale=yscale,
			 name = paste("node_barplot", node$nodeID, "plot", 
                         sep = ""))

        pushViewport(plot)
	
	if(beside) {
  	  xcenter <- cumsum(widths+gap) - widths/2
	  for (i in 1:np) {
            grid.rect(x = xcenter[i], y = 0, height = pred[i], 
                      width = widths[i],
	              just = c("center", "bottom"), default.units = "native",
	              gp = gpar(col = col[i], fill = fill[i]))
	  }
          if(length(xcenter) > 1) grid.xaxis(at = xcenter, label = FALSE)
	  grid.text(ylevels, x = xcenter, y = unit(-1, "lines"), 
                    just = c("center", "top"),
	            default.units = "native", check.overlap = TRUE)
          grid.yaxis()
	} else {
  	  ycenter <- cumsum(pred) - pred

	  for (i in 1:np) {
            grid.rect(x = xscale[2]/2, y = ycenter[i], height = min(pred[i], ymax - ycenter[i]), 
                      width = widths[1],
	              just = c("center", "bottom"), default.units = "native",
	              gp = gpar(col = col[i], fill = fill[i]))
	  }
          if(np > 1) {
	    grid.text(ylevels[1], x = unit(-1, "lines"), y = 0,
                      just = c("left", "center"), rot = 90,
	              default.units = "native", check.overlap = TRUE)
	    grid.text(ylevels[np], x = unit(-1, "lines"), y = ymax,
                      just = c("right", "center"), rot = 90,
	              default.units = "native", check.overlap = TRUE)
	  }
          if(np > 2) {
	    grid.text(ylevels[-c(1,np)], x = unit(-1, "lines"), y = ycenter[-c(1,np)],
                      just = "center", rot = 90,
	              default.units = "native", check.overlap = TRUE)
	  }
          grid.yaxis(main = FALSE)	
	}
	
        grid.rect(gp = gpar(fill = "transparent"))
        upViewport(2)
    }
    
    return(rval)
}
class(NodeBarplot) <- "grapcon_generator"


