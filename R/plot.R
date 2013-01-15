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
library(scales)

plot.gp.1d <- function(gp, s=NULL)
{
  col <- rgb(8/255, 81/255, 156/255, alpha=0.625)
  
  plot(gp$x, gp$mu, 'n', xlab="x", ylab="p", ylim=c(0,1))

  var <- diag(gp$sigma)
  z1  <- bound(gp$mu + 2*sqrt(var))
  z2  <- bound(gp$mu - 2*sqrt(var))
  
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
}

plot.gp.2d <- function(gp, f=NULL, main="")
{
  p1 <- ggplot(data = data.frame(x = gp$x[,1], y = gp$x[,2], z = gp$mu),
               aes(x = x, y = y, z = z)) +
        stat_contour() +
        geom_tile(aes(fill=z)) +
        scale_fill_gradient2(limits=c(0, 1), low=muted("green"), mid="white", high=muted("red"), midpoint=0.5) +
        ggtitle("Expected value")
  
  p2 <- ggplot(data = data.frame(x = gp$x[,1], y = gp$x[,2], z = diag(gp$sigma)),
               aes(x = x, y = y, z = z)) +
        geom_tile(aes(fill=z)) +
        stat_contour() +
        scale_fill_gradient(limits=c(0, 0.1), low="white", high=muted("green")) +
        ggtitle("Variance")
  
  if (!is.null(f)) {
    p3 <- ggplot(data = data.frame(x = gp$x[,1], y = gp$x[,2], z = f(gp$x)),
                 aes(x = x, y = y, z = z)) +
                 stat_contour() +
                 geom_tile(aes(fill=z)) +
                 scale_fill_gradient2(limits=c(0, 1), low=muted("green"), mid="white", high=muted("red"), midpoint=0.5) +
                 ggtitle("Ground truth")

    grid.arrange(p1, p2, p3, ncol=2, main=textGrob(main, vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)))
  }
  else {
    grid.arrange(p1, p2, ncol=2, main=textGrob(main, vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)))
  }
}

plot.gp <- function(gp, ...)
{
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
