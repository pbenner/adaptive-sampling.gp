# Copyright (C) 2013 Philipp Benner
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

library("ggplot2")
library("gridExtra")
library("scales")

plot.gp.1d <- function(gp, s=NULL)
{
  col <- rgb(8/255, 81/255, 156/255, alpha=0.625)
  
  p <- plot(gp$x, gp$mu, 'n', xlab="x", ylab="p", ylim=c(0,1))

  var <- diag(gp$sigma)
  if (is.null(gp$range)) {
    z1  <- gp$mu + 2*sqrt(var)
    z2  <- gp$mu - 2*sqrt(var)
  }
  else {
    z1  <- bound(gp$mu + 2*sqrt(var))
    z2  <- bound(gp$mu - 2*sqrt(var))
  }
  
  polygon(c(gp$x, rev(gp$x)), c(z1, rev(z2)),
     col = col, border = NA)

  lines(gp$x, gp$mu, 'l', lwd=3)

  if (!is.null(s)) {
    # if this is a scalar, then it specifies the number of samples
    if (!is.matrix(s)) {
      s <- samples(gp, s)
    }
    # otherwise it contains already samples
    for (i in 1:dim(s)[1]) {
      lines(gp$x, bound(s[i,]), 'l', lwd=0.5)
    }
  }
  return (p)
}

plot.gp.2d <- function(gp, counts=NULL, f=NULL, main="", plot.variance=TRUE)
{
  p1 <- ggplot(data = data.frame(x = gp$x[,1], y = gp$x[,2], z = gp$mu),
               aes_string(x = "x", y = "y", z = "z")) +
        geom_tile(aes_string(fill="z"), limits=c(min(gp$mu), max(gp$mu))) +
        stat_contour() +
        scale_fill_gradient2(limits=gp$range, low=muted("green"), mid="white", high=muted("red")) +
        ggtitle("Expected value")

  if (plot.variance) {
    p2 <- ggplot(data = data.frame(x = gp$x[,1], y = gp$x[,2], z = diag(gp$sigma)),
                 aes_string(x = "x", y = "y", z = "z")) +
                 geom_tile(aes_string(fill="z")) +
                 stat_contour() +
                 scale_fill_gradient(limits=c(0, 0.1), low="white", high=muted("green")) +
                 ggtitle("Variance")
  }
  else {
    p2 <- NULL
  }
  
  if (!is.null(f)) {
    p3 <- ggplot(data = data.frame(x = gp$x[,1], y = gp$x[,2], z = f(gp$x)),
                 aes_string(x = "x", y = "y", z = "z")) +
                 geom_tile(aes_string(fill="z")) +
                 stat_contour() +
                 scale_fill_gradient2(limits=gp$range, low=muted("green"), mid="white", high=muted("red"), midpoint=0.5) +
                 ggtitle("Ground truth")

    p4 <- ggplot(data = data.frame(x = counts[,1], y = counts[,2]),
                 aes_string(x = "x", y = "y")) +
                 stat_bin2d() +
                 scale_fill_gradient(low="white", high=muted("red")) +
                 ggtitle("Counts")
    
    g <- grid.arrange(p1, p2, p3, p4, ncol=2, main=textGrob(main, vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)))
  }
  else {
    if (plot.variance) {
      g <- grid.arrange(p1, p2, ncol=2, main=textGrob(main, vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)))
    }
    else {
      g <- grid.arrange(p1, ncol=1, main=textGrob(main, vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)))
    }
  }
  return (g)
}

#' Plot a Gaussian process
#' 
#' @param x Gaussian process
#' @param y unused
#' @param ... arguments to be passed to methods
#' @method plot gp
#' @S3method plot gp

plot.gp <- function(x, y=NULL, ...)
{
  gp <- x
  
  if (dim(gp) == 1) {
    p <- plot.gp.1d(gp, ...)
  }
  else if (dim(gp) == 2) {
    p <- plot.gp.2d(gp, ...)
  }
  else {
    stop("Gaussian process has invalid dimension.")
  }
  return (p)
}

#' Plot an experiment
#' 
#' @param x an object of class experiment
#' @param y positions where to plot the posterior
#' @param ... arguments to be passed to methods
#' @method plot experiment
#' @S3method plot experiment

plot.experiment <- function(x, y, ...)
{
  gp <- posterior(x, y)
  p  <- plot(gp, counts=get.counts(x), ...)
  return (p)
}
