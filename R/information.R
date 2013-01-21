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

#' Compute the Kullback-Leibler divergence between two models
#' 
#' @param model for instance a Gaussian process
#' @param ... arguments to be passed to methods
#' @export

kl.divergence <- function(model, ...)
{
  UseMethod("kl.divergence")
}

#' Compute the Kullback-Leibler divergence between two Gaussian
#' processes
#' 
#' @param model first Gaussian process
#' @param gp1 second Gaussian process
#' @param ... unused
#' @method kl.divergence gp
#' @S3method kl.divergence gp

kl.divergence.gp <- function(model, gp1, ...)
{
  gp0    <- model
  mu0    <- gp0$mu
  mu1    <- gp1$mu
  sigma0 <- gp0$sigma
  sigma1 <- gp1$sigma

  N      <- length(mu0)

  L0inv  <- solve(chol(sigma0))
  L1inv  <- solve(chol(sigma1))
  
  tmp1   <- log(det((sigma1 %*% L0inv) %*% t(L0inv)))
  tmp2   <- sum(diag(L1inv %*% (t(L1inv) %*% sigma0)))
  tmp3   <- (mu0 - mu1) %*% (L1inv %*% (t(L1inv) %*% (mu0 - mu1)))
  
  return (1/2*drop(tmp1 + tmp2 + tmp3 - N))
}
