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

#' @name adaptive.sampling.gp
#' @docType package
#' @title Implements a predictive approach to non-parametric inference for adaptive sampling
#' @author Philipp Benner \email{Philipp.Benner@@mis.mpg.de}, Tobias Elze \email{Tobias.Elze@@schepens.harvard.edu}
#' @references
#'  Poppe, S, Benner, P, Elze, T.
#'  A predictive approach to nonparametric inference for adaptive
#'  sequential sampling of psychophysical experiments.
#'  Journal of Mathematical Psychology 56 (2012) 179-195
#' @import mvtnorm
#' @import ggplot2
#' @import scales
#' @importFrom Matrix nearPD
#' @importFrom gridExtra grid.arrange
#' @importFrom nnet which.is.max
#' @useDynLib adaptive.sampling.gp
NULL

library("Matrix")

#' Generate a new Gaussian process
#' 
#' @param x positions of random variables
#' @param mu.prior prior mean
#' @param kernelf kernel function
#' @param range restrict the GP to a certain interval
#' @param sigma covariance matrix (optional)
#' @export

new.gp <- function(x, mu.prior, kernelf, range=NULL, sigma=NULL)
{
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  mu <- rep(mu.prior[1], dim(x)[1])
  
  gp <- list(x        = x,              # where to evaluate the gp
             kernelf  = kernelf,        # kernel function
             mu.prior = mu.prior,       # prior mean
             range    = range,          # restricted range
             mu       = mu,             # mean
             sigma    = sigma)
  class(gp) <- "gp"

  return (gp)
}

#' Show the dimension of the GP
#' 
#' @param x Gaussian process
#' @method dim gp
#' @S3method dim gp

dim.gp <- function(x) {
  dim(x$x)[2]
}

#' Compute posterior of a Gaussian process or an experiment
#' 
#' @param model Gaussian process or experiment
#' @param ... arguments to be passed to methods
#' @export

posterior <- function(model, ...)
{
  UseMethod("posterior")
}

#' Compute posterior of a Gaussian process
#' 
#' @param model Gaussian process
#' @param xp positions of measurements
#' @param yp measured values
#' @param noise uncertainty of measurements (optional)
#' @param ... unused
#' @method posterior gp
#' @S3method posterior gp

posterior.gp <- function(model, xp, yp, noise=NULL, ...)
{
  gp <- model
  # check arguments
  if (!is.matrix(xp)) {
    xp <- as.matrix(xp)
  }
  if (!is.matrix(yp)) {
    yp <- as.matrix(yp)
  }
  if (!is.null(noise) && !is.matrix(noise)) {
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
    A <- k0 + diag(rep(noise, dim(xp)[1]))
  }
  else {
    A <- k0 + diag(as.vector(noise))
  }
  L        <- chol(A)
  Linv     <- solve(L)

  gp$mu    <- drop(mu + (k2 %*% Linv) %*% (t(Linv) %*% (yp - mu)))
  gp$sigma <- k3 - (k2 %*% Linv) %*% (t(Linv) %*% k1)

  if (!is.null(gp$range)) {
    gp$mu  <- bound(gp$mu, gp$range)
  }

  # the resulting matrix is usually not positive definite due to
  # numerical errors; use nearPD to compute the nearest positive
  # definite matrix
  gp$sigma <- as.matrix(nearPD(gp$sigma)$mat)

  return (gp)
}
