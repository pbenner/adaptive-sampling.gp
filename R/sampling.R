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
library("nnet") # which.is.max

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
#' @param model Gaussian process or experment
#' @param n number of samples
#' @param ... unused
#' @method samples gp
#' @S3method samples gp

samples.gp <- function(model, n=1, ...)
{
  gp <- model
  x  <- rmvnorm(n=n, mean=gp$mu, sigma=gp$sigma, method="chol")

  return (x)
}

#' Compute the utility of an experiment
#' 
#' @param experiment an object of class experiment
#' @param ... arguments to be passed to methods
#' @export

utility <- function(experiment, ...)
{
  UseMethod("utility")
}

#' Compute the utility of an experiment
#' 
#' @param experiment an object of class experiment
#' @param x positions where to evaluate the experiment
#' @param ... unused
#' @method utility experiment
#' @S3method utility experiment

utility.experiment <- function(model, x, ...)
{
  experiment <- model
  # select an appropriate method according to the type of the experiment
  if (experiment$type == "bernoulli") {
    return (utility.experiment.bernoulli(model, x, ...))
  }
  if (experiment$type == "gaussian") {
    return (utility.experiment.gaussian(model, x, ...))
  }
  stop("Unknown experiment type.")
}

utility.experiment.bernoulli <- function(experiment, x, ...)
{
  if (is.vector(x)) {
    x <- as.matrix(x)
  }
  
  L   <- dim(x)[1] # number of possible stimuli
  gp0 <- posterior(experiment, x, enforce.pd=TRUE)
  ut  <- as.matrix(rep(0.0, L))
  p   <- predictive(gp0)$mean
  
  for (i in 1:L) {
    add.measurement(experiment, x[i,], c( 1, 0))
    gp1 <- posterior(experiment, x, enforce.pd=TRUE)
    add.measurement(experiment, x[i,], c(-1, 1))
    gp2 <- posterior(experiment, x, enforce.pd=TRUE)
    add.measurement(experiment, x[i,], c( 0,-1))

    ut[i] <- ut[i] +    p[i] *kl.divergence(gp0, gp1)
    ut[i] <- ut[i] + (1-p[i])*kl.divergence(gp0, gp2)
  }
  return (ut)
}

# TODO: Is this version equivalent to the above one?
utility.experiment.bernoulli <- function(experiment, x, ...)
{
  if (is.vector(x)) {
    x <- as.matrix(x)
  }
  
  L   <- dim(x)[1] # number of possible stimuli
  ut  <- as.matrix(rep(0.0, L))
  
  for (i in 1:L) {
    gp0 <- posterior(experiment, x[i,], enforce.pd=TRUE)
    p   <- predictive(gp0)$mean
    add.measurement(experiment, x[i,], c( 1, 0))
    gp1 <- posterior(experiment, x[i,], enforce.pd=TRUE)
    add.measurement(experiment, x[i,], c(-1, 1))
    gp2 <- posterior(experiment, x[i,], enforce.pd=TRUE)
    add.measurement(experiment, x[i,], c( 0,-1))

    ut[i] <- ut[i] +    p *kl.divergence(gp0, gp1)
    ut[i] <- ut[i] + (1-p)*kl.divergence(gp0, gp2)
  }
  return (ut)
}

utility.experiment.gaussian <- function(experiment, x, ...)
{
  if (is.vector(x)) {
    x <- as.matrix(x)
  }
  
  L   <- dim(x)[1] # number of possible stimuli
  ut  <- as.matrix(rep(0.0, L))

  for (i in 1:L) {
    
  }
  return (ut)
}

#' Simulate an experiment by taking samples from a ground truth
#' 
#' @param experiment an object of class experiment
#' @param x positions where to evaluate the experiment
#' @param gt the ground truth
#' @param sdf for a Gaussian experiment sdf is a function that determines how much
#' we trust a measurement at a given position
#' @param N number of samples
#' @param verbose print sampling step
#' @export

sample.with.gt <- function(experiment, x, gt, sdf=function(x, mean) 1, N=1, verbose=FALSE)
{
  if (is.vector(x)) {
    x <- as.matrix(x)
  }
  for (i in 1:N) {
    if (verbose) {
      print(sprintf("Sampling step... %d", i))
    }
    # compute utility for every position
    ut <- utility(experiment, x)
    # select best position; if multiple global maxima exist then
    # choose one of them at random
    k  <- which.is.max(ut)
    # draw a new sample from the ground truth
    if (experiment$type == "bernoulli") {
      add.measurement(experiment, x[k,], gt(x[k,]))
    }
    if (experiment$type == "gaussian") {
      # draw a mean from the ground truth
      mean <- gt(x[k,])
      # how much do we trust this measurement?
      sd   <- sdf(x[k,])
      # add this to the experimental data
      add.measurement(experiment, x[k,], c(mean, sd))
    }
  }
}

#' Generate a new ground truth from a vector or function of binomial
#' parameters
#' 
#' @param f function or vector defining the binomial parameters
#' @param ... arguments to be passed to methods
#' @export

new.gt <- function(f, ...)
{
  UseMethod("new.gt")
}

#' Generate a new ground truth from a vector of binomial parameters
#' 
#' @param f positions
#' @param y values of the ground truth
#' @param ... unused
#' @method new.gt numeric
#' @S3method new.gt numeric

new.gt.numeric <- function(f, y, ...)
{
  new.gt.matrix(as.matrix(f), as.matrix(y), ...)
}

#' Generate a new ground truth from a matrix of binomial parameters
#' 
#' @param f positions
#' @param y values of the ground truth
#' @param type either "bernoulli" or "gaussian"
#' @param sd standard deviation if type is gaussian
#' @param ... unused
#' @method new.gt matrix
#' @S3method new.gt matrix

new.gt.matrix <- function(f, y, type="bernoulli", sd=NULL, ...)
{
  xp <- f

  if (type == "bernoulli") {
    gt <- function(xt) {
      p <- apply(xp, 1, function(x) all(x == as.matrix(xt)))
      if (runif(1, 0, 1) <= y[p]) {
        return (c(1,0))
      }
      else {
        return (c(0,1))
      }
    }
  }
  else {
    gt <- function(xt) {
      p <- apply(xp, 1, function(x) all(x == as.matrix(xt)))
      if (length(sd) == 1) {
        rnorm(1, mean=y[p], sd=sd)
      }
      else {
        rnorm(1, mean=y[p], sd=sd[p])
      }
    }
  }
  return (gt)
}

#' Generate a new ground truth from a function
#' 
#' @param f function that specifies the ground truth
#' @param ... unused
#' @method new.gt function
#' @S3method new.gt function

new.gt.function <- function(f, ...)
{
  gt <- function(xt) {
    if (runif(1, 0, 1) <= f(xt)) {
      return (c(1,0))
    }
    else {
      return (c(0,1))
    }
  }
  return (gt)
}
