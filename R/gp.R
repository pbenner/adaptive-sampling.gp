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
#' @param link link (response) function
#' @export

new.gp <- function(x, mu.prior, kernelf, range=NULL, link=NULL)
{
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  
  gp <- list(x        = x,              # where to evaluate the gp
             kernelf  = kernelf,        # kernel function
             mu.prior = mu.prior,       # prior mean
             range    = range,          # restricted range
             link     = link,           # link (response) function
             # these are computed by posterior()
             mu       = NULL,           # posterior mean
             sigma    = NULL)           # posterior variance
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

#' Compute predictive posterior expectation and variance
#' 
#' @param model GP object
#' @param ... unused
#' @method predictive gp
#' @S3method predictive gp

predictive.gp <- function(model, ...)
{
  gp <- model
  
  return (predictive(gp$link, gp$mu, diag(gp$sigma)))
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
  # xp and yp can be NULL
  if (is.null(xp) || is.null(yp)) {
    gp$mu    <- as.matrix(rep(gp$mu.prior, dim(gp$x)[1]))
    gp$sigma <- gp$kernelf(gp$x, gp$x)
    gp$sigma <- as.matrix(nearPD(gp$sigma)$mat)
    return (gp)
  }
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

  # check if a link is available
  if (!is.null(gp$link)) {
    # if so, then call a specialized method
    # approximate the posterior of the gaussian process
    result   <- approximate.posterior(xp, yp, k0, gp$link)
    # which is a gaussian with mean mu and covariance sigma
    gp$mu    <- k2 %*% result$d
    v        <- solve(result$L) %*% (sqrt(result$W) %*% k1)
    gp$sigma <- k3 - t(v) %*% v
  }
  else {
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

    gp$mu    <- drop(gp$mu.prior + (k2 %*% Linv) %*% (t(Linv) %*% (yp - gp$mu.prior)))
    gp$sigma <- k3 - (k2 %*% Linv) %*% (t(Linv) %*% k1)

    if (!is.null(gp$range)) {
      gp$mu  <- bound(gp$mu, gp$range)
    }
  }
  # the resulting matrix is usually not positive definite due to
  # numerical errors; use nearPD to compute the nearest positive
  # definite matrix
  gp$sigma <- as.matrix(nearPD(gp$sigma)$mat)

  return (gp)
}

approximate.posterior.derivative <- function(f, yp, K, link, N)
{
  # d: d/dx log p(y|f)
  d <- as.matrix(rep(0, N))
  # W: negative Hessian of log p(y|f)
  W <- diag(N)
  for (i in 1:N) {
    # current f value at x[[i]]
    fx  <- f[[i]]
    # counts at position x[[i]]
    c1  <- yp[[i,1]]
    c2  <- yp[[i,2]]
    # value of the response derivative evaluated at fx
    Nfx <- link$response.derivative(fx)
    # response evaluated at fx
    Pfx <- link$response(fx)
    d[[i]]   <- c1*Nfx/(0+Pfx) - c2*Nfx/(1-Pfx)
    W[[i,i]] <- c2*(-fx*Nfx/(1-Pfx) + Nfx^2/(1-Pfx)^2) -
                c1*(-fx*Nfx/(0+Pfx) - Nfx^2/(0+Pfx)^2)
  }
  return (list(d = d, W = W))
}

approximate.posterior.step <- function(f, yp, K, link, N)
{
  derivative <- approximate.posterior.derivative(f, yp, K, link, N)
  # d: d/dx log p(y|f)
  d <- derivative$d
  # W: negative Hessian of log p(y|f)
  W <- derivative$W
  B <- diag(N) + sqrt(W) %*% K %*% sqrt(W)
  L <- t(chol(B))
  b <- W %*% f + d
  a <- b - sqrt(W) %*% solve(t(L)) %*% (solve(L) %*% (sqrt(W) %*% K %*% b))
  f <- K %*% a

  return (f)
}

approximate.posterior <- function(xp, yp, K, link, epsilon=0.00001)
{
  # number of positions where measurements are available
  N <- dim(xp)[[1]]
  # f, fold
  f     <- as.matrix(rep(0, N))
  f.old <- as.matrix(rep(0, N))
  repeat {
    # run Newton steps until convergence
    f <- approximate.posterior.step(f, yp, K, link, N)
    if (norm(f - f.old) < epsilon) {
      break
    }
    f.old <- f
  }
  # evaluate the derivative at the current position
  derivative <- approximate.posterior.derivative(f, yp, K, link, N)
  d <- derivative$d
  W <- derivative$W
  B <- diag(N) + sqrt(W) %*% K %*% sqrt(W)
  L <- t(chol(B))
  # and compute the result, what is needed later...
  result <- list(# for the expectation we need the derivative
                 d = d,
                 # and for the variance L and W
                 L = L,
                 W = W)

  return (result)
}
