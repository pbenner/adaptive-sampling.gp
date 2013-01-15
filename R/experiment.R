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

add.measurement <- function(experiment, ...)
{
  UseMethod("add.measurement")
}

add.measurement.experiment <- function(experiment, x, counts)
{
  key <- toString(x)

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

dirichlet.moments <- function(alpha)
{
  expectation <- rep(0, length(alpha))
  variance    <- rep(0, length(alpha))

  alpha0      <- sum(alpha)
  
  for (i in 1:length(alpha)) {
    expectation[i] <- alpha[i]/alpha0
    variance[i]    <- alpha[i]*(alpha0 - alpha[i])/(alpha0^2*(alpha0 + 1))
  }
  
  result <- list(expectation = expectation, variance = variance)

  return (result)
}

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
    xp <- c() # position
    yp <- c() # mean
    ep <- c() # variance

    for (key in ls(envir=experiment$data)) {
      xt      <- as.numeric(key)
      counts  <- experiment$data[[key]]

      moments <- dirichlet.moments(counts + alpha)

      xp <- append(xp, xt)
      yp <- append(yp, moments$expectation[1])
      ep <- append(ep, moments$variance[1])
    }
    # compute the posterior of the gaussian process
    gp <- posterior(gp, xp, yp, ep)
  }

  return (gp)
}

plot.experiment <- function(experiment, x)
{
  gp <- posterior(e, x)
  plot(gp)
}

# Example
################################################################################

e <- new.experiment(c(2.0,2.0))
add.measurement(e, 1, c(100,3))
add.measurement(e, 2, c( 90,1))
add.measurement(e, 3, c(100,4))

gp <- posterior(e, 1:100/20)
plot(gp)

e <- new.experiment(c(2.0,2.0))
add.measurement(e, 1, c( 1,4))
add.measurement(e, 2, c( 1,3))
add.measurement(e, 3, c( 1,6))

gp <- posterior(e, 1:100/20)
plot(gp, samples(gp, 10))
