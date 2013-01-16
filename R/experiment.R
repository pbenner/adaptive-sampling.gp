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
#' @exportClass experiment
#' @export

new.experiment <- function(alpha,
                           # the kernel.type is a partially applied
                           # kernel function, which still has a free
                           # variance parameter
                           kernel.type=partially.apply(kernel.exponential, 1.0))
{
  experiment        <- list(alpha       = alpha,       # Dirichlet pseudo counts
                            data        = new.env(),   # experimental data
                            kernel.type = kernel.type) # kernel function
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
#' @method add.measurement experiment

add.measurement.experiment <- function(experiment, x, counts)
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
#' @param experiment an object of class experiment
#' @param x positions where to evaluate the posterior
#' @method posterior experiment

posterior.experiment <- function(experiment, x)
{
  # get prior pseudocounts
  alpha             <- experiment$alpha
  # determine prior expectation and standard deviation
  prior.moments     <- dirichlet.moments(alpha)
  # the prior expectation sets the mean of the GP
  prior.expectation <- prior.moments$expectation[1]
  # and the variance is used in the kernel function
  prior.variance    <- prior.moments$variance[1]
  # generate a kernel which is limited by the prior variance
  kernelf           <- experiment$kernel.type(prior.variance)
  # construct the gaussian process
  gp                <- new.gp(x, prior.expectation, kernelf)

  if (length(experiment$data) > 0) {
    # we have measurements...
    xp <- matrix(0, dim(gp), nrow=length(experiment$data)) # position
    yp <- matrix(0, 1,       nrow=length(experiment$data)) # mean
    ep <- matrix(0, 1,       nrow=length(experiment$data)) # variance

    for (i in 1:length(experiment$data)) {
      key     <- ls(envir=experiment$data)[i]
      xt      <- key2value(key)
      counts  <- experiment$data[[key]]

      moments <- dirichlet.moments(counts + alpha)

      xp[i,] <- xt
      yp[i,] <- moments$expectation[1]
      ep[i,] <- moments$variance[1]
    }
    # compute the posterior of the gaussian process
    gp <- posterior(gp, xp, yp, ep)
  }

  return (gp)
}
