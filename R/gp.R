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

library("mvtnorm")
library("Matrix")
library("ggplot2")

new.gp <- function(x, mu.prior, kernelf, sigma=NULL)
{
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  mu <- rep(mu.prior[1], dim(x)[1])

  if (is.null(sigma)) {
    sigma <- as.matrix(nearPD(kernelf(x, x))$mat)
  }
  
  gp <- list(x        = x,              # where to evaluate the gp
             kernelf  = kernelf,        # kernel function
             mu.prior = mu.prior,       # prior mean
             mu       = mu,             # mean
             sigma    = sigma)
  class(gp) <- "gp"

  return (gp)
}

partially.apply <- function(f, ...) {
  capture <- list(...)
  function(...) {
    do.call(f, c(capture, list(...)))
  }
}

posterior <- function(gp, ...)
{
  UseMethod("posterior")
}

posterior.gp <- function(gp, xp, yp, noise=NULL)
{
  # check arguments
  if (!is.matrix(xp)) {
    xp <- as.matrix(xp)
  }
  if (!is.matrix(yp)) {
    yp <- as.matrix(yp)
  }
  if (!is.matrix(noise)) {
    noise <- as.matrix(noise)
  }
  # compute posterior...
  k0 <- gp$kernelf(xp, xp)     # K(X , X )
  k1 <- gp$kernelf(xp, gp$x)   # K(X , X*)
  k2 <- t(k1)                  # K(X*, X )
  k3 <- gp$kernelf(gp$x, gp$x) # K(X*, X*)

  mu <- gp$mu.prior

  # add noise to measurements?
  if (is.null(noise)) {
    A <- k0
  }
  else if (dim(noise)[1] == 1) {
    A <- k0 + noise
  }
  else {
    A <- k0 + diag(as.vector(noise))
  }
  L        <- chol(A)
  Linv     <- solve(L)

  gp$mu    <- drop(mu + (k2 %*% Linv) %*% (t(Linv) %*% (yp - mu)))
  gp$sigma <- k3 - (k2 %*% Linv) %*% (t(Linv) %*% k1)

  # the resulting matrix is usually not positive definite due to
  # numerical errors; use nearPD to compute the nearest positive
  # definite matrix
  gp$sigma <- as.matrix(nearPD(gp$sigma)$mat)

  return (gp)
}

samples <- function(gp, ...)
{
  UseMethod("samples")
}

samples.gp <- function(gp, n=1)
{
  x <- rmvnorm(n=n, mean=gp$mu, sigma=gp$sigma, method="chol")

  return (x)
}

bound <- function(line, range=c(0,1))
{
  # bound function values in the inverval given by range
  tmp0 <- sapply(line, function(x) min(x, range[2]))
  tmp1 <- sapply(tmp0, function(x) max(x, range[1]))
  return (tmp1)
}

plot.gp <- function(gp, s=NULL)
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

# Example
################################################################################

gp <- new.gp(1:100/20, 0.5, kernel.exponential(1, 1))

xp <- c(1, 2, 3)
yp <- c(0.7, 0.7, 0.7)
# measurement noise
ep <- c(0.01, 0.01, 0.01)

gp <- posterior(gp, xp, yp, ep)

plot(gp)

plot(gp, samples(gp, 100))

# 2D
################################################################################

data     <- as.matrix(expand.grid(x = 1:20/4, y = 1:20/4))
gp       <- new.gp(data, 0.5, kernel.exponential(2, 1))

xp <- matrix(0,2, nrow=2)
xp[1,] <- c(2,2)
xp[1,] <- c(4,4)

yp <- c(0.2, 0.8)
# measurement noise
ep <- c(0.001, 0.001)

gp <- posterior(gp, xp, yp, ep)


p <- ggplot(data = data.frame(x = data[,1], y = data[,2], z = gp$mu),
            aes(x = x, y = y, z = z))
p + geom_tile(aes(fill=z)) + stat_contour()

p <- ggplot(data = data.frame(x = data[,1], y = data[,2], z = diag(gp$sigma)),
            aes(x = x, y = y, z = z))
p + geom_tile(aes(fill=z)) + stat_contour()

