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
#' @param prior.mean prior mean of the Gaussian process
#' @param kernelf kernel function (prior covariance matrix)
#' @param type type of the experiment (e.g. bernoulli or gaussian)
#' @export

new.experiment <- function(prior.mean=0.5,
                           kernelf=kernel.exponential(1, 0.25),
                           type="bernoulli")
{
  experiment        <- list(data       = new.env(), # experimental data
                            kernelf    = kernelf,   # kernel function
                            type       = type,      # type of the experiment
                            prior.mean = prior.mean)
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

add.measurement.experiment <- function(experiment, x, data, ...)
{
  # select an appropriate method according to the type of the experiment
  if (experiment$type == "bernoulli") {
    return (add.measurement.experiment.bernoulli(experiment, x, data, ...))
  }
  if (experiment$type == "gaussian") {
    return (add.measurement.experiment.gaussian(experiment, x, data, ...))
  }
  stop("Unknown experiment type.")
}

add.measurement.experiment.bernoulli <- function(experiment, x, data, ...)
{
  key <- value2key(x)

  if (is.null(experiment$data[[key]])) {
    experiment$data[[key]] <- data
  }
  else {
    # for bernoulli experiments, simply sum up the counts
    experiment$data[[key]] <- experiment$data[[key]] + data
  }
  if (all(experiment$data[[key]] == 0)) {
    rm(list = key, envir=experiment$data)
  }
}

add.measurement.experiment.gaussian <- function(experiment, x, data, ...)
{
  key <- value2key(x)

  # convert data to matrix if necessary
  if (is.numeric(data)) {
    data <- t(data)
  }
  if (is.null(experiment$data[[key]])) {
    experiment$data[[key]] <- data
  }
  else {
    # with gaussian experiments, save all recordings in a matrix
    experiment$data[[key]] <- rbind(experiment$data[[key]], data)
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
  # select an appropriate method according to the type of the experiment
  if (experiment$type == "bernoulli") {
    return (posterior.experiment.bernoulli(model, x, ...))
  }
  if (experiment$type == "gaussian") {
    return (posterior.experiment.gaussian(model, x, ...))
  }
  stop("Unknown experiment type.")
}

posterior.experiment.bernoulli <- function(model, x, ...)
{
  experiment <- model
  # construct the gaussian process with a link function
  link       <- new.link()
  mu         <- link$link(experiment$prior.mean)
  gp         <- new.gp(x, mu, experiment$kernelf, range=c(0,1), link=link)

  xp         <- NULL
  yp         <- NULL

  # if we have measurements, add them to xp and yp
  if (length(experiment$data) > 0) {
    xp <- matrix(0, dim(gp), nrow=length(experiment$data)) # position
    yp <- matrix(0, 2,       nrow=length(experiment$data)) # counts

    for (i in 1:length(experiment$data)) {
      key     <- ls(envir=experiment$data)[i]
      xt      <- key2value(key)
      counts  <- experiment$data[[key]]

      xp[i,] <- xt
      yp[i,] <- counts
    }
  }
  # compute the posterior, whether or not measurements
  # exist
  gp <- posterior(gp, xp, yp)
  gp
}

posterior.experiment.gaussian <- function(model, x, ...)
{
  experiment <- model
  # do not use a link function for gaussian experiments
  gp         <- new.gp(x, experiment$prior.mean, experiment$kernelf)

  xp         <- NULL # position
  yp         <- NULL # mean
  ep         <- NULL # variance

  # if we have measurements, add them to xp, yp, and ep
  if (length(experiment$data) > 0) {
    # loop through entries
    for (i in 1:length(experiment$data)) {
      key <- ls(envir=experiment$data)[i]
      # for each entry, there might be several measurements
      for (j in 1:nrow(experiment$data[[key]])) {
        xp <- rbind(xp, key2value(key))
        yp <- rbind(yp, experiment$data[[key]][j, 1])
        ep <- rbind(ep, experiment$data[[key]][j, 2])
      }
    }
  }
  # compute the posterior, whether or not measurements
  # exist
  gp <- posterior(gp, xp, yp, noise=ep)
  gp
}
