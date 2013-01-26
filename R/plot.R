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

plot.gp.1d <- function(gp, s=NULL, xlab=NULL, ylab=NULL, ...)
{
  result <- predictive(gp)
  
  col <- rgb(10/255, 40/255, 10/255, alpha=0.5)

  if (is.null(gp$range)) {
    z1  <- result$mean + 2*sqrt(result$variance)
    z2  <- result$mean - 2*sqrt(result$variance)
  }
  else {
    z1  <- bound(result$mean + 2*sqrt(result$variance))
    z2  <- bound(result$mean - 2*sqrt(result$variance))
  }

  # plot returns nothing, so don't bother to catch the
  # return value
  plot(cbind(gp$x, gp$x, gp$x), cbind(result$mean, z1, z2),
       'n', xlab=xlab, ylab=ylab, ylim=gp$range)
  
  polygon(c(gp$x, rev(gp$x)), c(z1, rev(z2)),
     col = col, border = NA)

  lines(gp$x, result$mean, 'l', lwd=3)

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
}

plot.gp.2d <- function(gp, counts=NULL, f=NULL, main="", plot.variance=TRUE,
                       low=muted("green"), mid="white", high=muted("red"),
                       midpoint=NULL, ...)
{
  # initialize all plot objects
  p1 <- NULL
  p2 <- NULL
  p3 <- NULL
  p4 <- NULL
  # compute predictive expectation and variance
  result <- predictive(gp)
  # a midpoint is computed from gp$range, but only if it is
  # not not given as an argument to this function
  if (is.null(midpoint) && !is.null(gp$range)) {
    midpoint <- sum(gp$range)/2
  }

  # first, plot the expectation
  p1 <- ggplot(data = data.frame(x = gp$x[,1], y = gp$x[,2], z = result$mean),
               aes_string(x = "x", y = "y", z = "z")) +
        geom_tile(aes_string(fill="z"), limits=c(min(result$mean), max(result$mean))) +
        stat_contour() +
        scale_fill_gradient2(limits=gp$range, low=low, mid=mid, high=high, midpoint=midpoint) +
        ggtitle("Expected value")

  # plot varience only if plot.variance is TRUE
  if (plot.variance) {
    p2 <- ggplot(data = data.frame(x = gp$x[,1], y = gp$x[,2], z = result$variance),
                 aes_string(x = "x", y = "y", z = "z")) +
                 geom_tile(aes_string(fill="z")) +
                 stat_contour() +
                 scale_fill_gradient(low=mid, high=high) +
                 ggtitle("Variance")
  }
  # check whether a ground truth is available
  if (!is.null(f)) {
    p3 <- ggplot(data = data.frame(x = gp$x[,1], y = gp$x[,2], z = f(gp$x)),
                 aes_string(x = "x", y = "y", z = "z")) +
                 geom_tile(aes_string(fill="z")) +
                 stat_contour() +
                 scale_fill_gradient2(limits=gp$range, low=low, mid=mid, high=high, midpoint=midpoint) +
                 ggtitle("Ground truth")
  }
  # check if counts are available
  if (!is.null(counts)) {
    p4 <- ggplot(data = data.frame(x = counts[,1], y = counts[,2]),
                 aes_string(x = "x", y = "y")) +
                 stat_bin2d() +
                 scale_fill_gradient(low=mid, high=high) +
                 ggtitle("Counts")
  }
  # use different grids depending on what plots are available
  if (!is.null(p2) && !is.null(p3) && !is.null(p4)) {
    grid.arrange(p1, p2, p3, p4, ncol=2, main=textGrob(main, vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)))
  }
  else if (!is.null(p2)) {
    grid.arrange(p1, p2, ncol=2, main=textGrob(main, vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)))
  }
  else {
    grid.arrange(p1, ncol=1, main=textGrob(main, vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)))
  }
  # return a list of all plot objects
  return (invisible(list(p1=p1, p2=p2, p3=p3, p4=p4)))
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
    plot.gp.1d(gp, ...)
  }
  else if (dim(gp) == 2) {
    plot.gp.2d(gp, ...)
  }
  else {
    stop("Gaussian process has invalid dimension.")
  }
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
  plot(gp, counts=get.counts(x), ...)
}
