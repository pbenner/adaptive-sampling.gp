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
  # construct the gaussian process with a link function
  gp         <- new.gp(x, 0.0, experiment$kernelf, range=c(0,1), link=new.link())

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
    posterior(gp, xp, yp)
  }

  return (gp)
}
