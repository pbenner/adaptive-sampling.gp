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

new.link <- function(type = "probit")
{
  if (type == "probit") {
    # qnorm: quantile function  (qnorm = pnorm^-1)
    # pnorm: cumulative density (pnorm = Int dnorm)
    # dnorm: density function   (dnorm = d/dx pnorm)
    result   <- list(link                = qnorm, # link function
                     response            = pnorm, # inverse link function (response)
                     response.derivative = dnorm) # derivative of the response function
    class(result) <- c("link.probit", "link")
  }
  else {
    stop("Unknown type.")
  }
  return (result)
}

#' Compute predictive posterior expectation and variance
#' 
#' @param model Gaussian process or link function
#' @param ... arguments to be passed to methods
#' @export

predictive <- function(model, ...)
{
  UseMethod("predictive")
}

#' Compute predictive posterior expectation and variance
#' 
#' @param model NULL object
#' @param p.mean posterior mean
#' @param p.variance posterior variance
#' @param ... unused
#' @method predictive NULL
#' @S3method predictive NULL

predictive.NULL <- function(model, p.mean, p.variance, ...)
{
  return (list(mean = p.mean, variance = p.variance))
}

#' Compute predictive posterior expectation and variance
#' 
#' @param model probit link object
#' @param p.mean posterior mean
#' @param p.variance posterior variance
#' @param ... unused
#' @method predictive link.probit
#' @S3method predictive link.probit

predictive.link.probit <- function(model, p.mean, p.variance, ...)
{
  N        <- dim(p.mean)[1]
  mean     <- as.matrix(rep(0, N))
  variance <- as.matrix(rep(0, N))

  for (i in 1:N) {
    # prepare the covariance matrix
    sigma <- matrix(c(1+p.variance[i], p.variance[i], p.variance[i], 1+p.variance[i]), 2, 2)
    # predictive mean
    mean[i]     <- pnorm(p.mean[i], mean=0, sd=sqrt(1+p.variance[i]))
    # and predictive variance
    variance[i] <- as.numeric(pmvnorm(upper=c(p.mean[i], p.mean[i]), mean=c(0, 0), sigma=sigma)) - mean[i]^2
  }
  return (list(mean = mean, variance = variance))
}
