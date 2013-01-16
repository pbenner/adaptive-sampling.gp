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

library("mvtnorm")
library("Matrix")

#' Generate a new Gaussian process
#' 
#' @param x positions of random variables
#' @param mu.prior prior mean
#' @param kernelf kernel function
#' @param sigma covariance matrix (optional)
#' @exportClass gp
#' @export

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

#' Show the dimension of the GP
#' 
#' @param gp Gaussian process
#' @method dim gp

dim.gp <- function(gp) {
  dim(gp$x)[2]
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
#' @param gp Gaussian process
#' @param xp positions of measurements
#' @param yp measured values
#' @param noise uncertainty of measurements (optional)
#' @method posterior gp

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

#' Draw samples from a process model
#' 
#' @param model Gaussian process or experment
#' @param ... arguments to be passed to methods
#' @export

samples <- function(model, ...)
{
  UseMethod("samples")
}

#' Draw samples from a process model
#' 
#' @param gp Gaussian process or experment
#' @param n number of samples
#' @method samples gp

samples.gp <- function(gp, n=1)
{
  x <- rmvnorm(n=n, mean=gp$mu, sigma=gp$sigma, method="chol")

  return (x)
}
