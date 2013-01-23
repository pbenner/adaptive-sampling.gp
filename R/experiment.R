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

key2value <- function(key) drop(sapply(strsplit(key, ","), FUN=as.numeric))
value2key <- function(value) toString(value)

#' Generate a new experiment
#' 
#' @param alpha Dirichlet pseudocounts
#' @param kernel.type partially applied kernel function
#' @export

new.experiment <- function(kernelf=kernel.exponential(1, 1))
{
  experiment        <- list(data    = new.env(),# experimental data
                            kernelf = kernelf)  # kernel function
  class(experiment) <- "experiment"

  return (experiment)
}

#' Add measurements
#' 
#' @param experiment an object of class experiment
#' @param ... arguments to be passed to methods
#' @export

add.measurement <- function(experiment, ...)
{
  UseMethod("add.measurement")
}

#' Add measurements
#' 
#' @param experiment an object of class experiment
#' @param x position of the measurement
#' @param counts measured binomial data
#' @param ... unused
#' @method add.measurement experiment
#' @S3method add.measurement experiment

add.measurement.experiment <- function(experiment, x, counts, ...)
{
  key <- value2key(x)

  if (is.null(experiment$data[[key]])) {
    experiment$data[[key]] <- counts
  }
  else {
    experiment$data[[key]] <- experiment$data[[key]] + counts
  }
  if (all(experiment$data[[key]] == 0)) {
    rm(list = key, envir=experiment$data)
  }
}

get.counts <- function(experiment)
{
  counts <- NULL
  for (key in ls(envir=experiment$data)) {
    value <- key2value(key)
    for (i in 1:sum(experiment$data[[key]])) {
      counts <- rbind(counts, value)
    }
  }
  return (counts)
}

#' Compute posterior of an experiment
#' 
#' @param model an object of class experiment
#' @param x positions where to evaluate the posterior
#' @param ... unused
#' @method posterior experiment
#' @S3method posterior experiment

posterior.experiment <- function(model, x, ...)
{
  experiment <- model
  # construct the gaussian process
  gp         <- new.gp(x, 0.0, experiment$kernelf)
  link       <- new.link()

  if (length(experiment$data) > 0) {
    # we have measurements...
    xp <- matrix(0, dim(gp), nrow=length(experiment$data)) # position
    yp <- matrix(0, 2,       nrow=length(experiment$data)) # counts

    for (i in 1:length(experiment$data)) {
      key     <- ls(envir=experiment$data)[i]
      xt      <- key2value(key)
      counts  <- experiment$data[[key]]

      xp[i,] <- xt
      yp[i,] <- counts
    }
    # evaluate the kernel
    k0 <- experiment$kernelf(xp, xp)  # K(X , X )
    k1 <- experiment$kernelf(xp, x )  # K(X , X*)
    k2 <- t(k1)                       # K(X*, X )
    k3 <- experiment$kernelf(x,  x )  # K(X*, X*)
    # approximate the posterior of the gaussian process
    result   <- approximate.posterior(xp, yp, k0, link)
    # which is a gaussian with mean mu and covariance sigma
    gp$mu    <- k2 %*% result$d
    v        <- solve(result$L) %*% (sqrt(result$W) %*% k1)
    gp$sigma <- k3 - t(v) %*% v

    # link...
    mu <- gp$mu
    s  <- as.matrix(diag(gp$sigma))
    z  <- mu/sqrt(1+s)
    gp$mu <- pnorm(mu + s * dnorm(z) / (pnorm(z) * sqrt(1 + s)))
  }

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
