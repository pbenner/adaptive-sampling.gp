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

library(nnet) # which.is.max

utility <- function(experiment, ...)
{
  UseMethod("utility")
}

utility.experiment <- function(experiment, x)
{
  L   <- length(x) # number of possible stimuli
  gp0 <- posterior(experiment, x)
  ut  <- rep(0.0, L)
  
  for (i in 1:L) {
    add.measurement(experiment, x[i], c( 1, 0))
    gp1 <- posterior(experiment, x)
    add.measurement(experiment, x[i], c(-1, 1))
    gp2 <- posterior(experiment, x)
    add.measurement(experiment, x[i], c( 0,-1))

    p     <- bound(gp0$mu[i], c(0, 1))
    ut[i] <- ut[i] +    p *kl.divergence(gp0, gp1)
    ut[i] <- ut[i] + (1-p)*kl.divergence(gp0, gp2)
  }
  return (ut)
}

sample.with.gt <- function(experiment, x, gt, N=1)
{
  for (i in 1:N) {
    # compute utility for every position
    ut <- utility(experiment, x)
    # select best position; if multiple global maxima exist then
    # choose one of them at random
    k  <- which.is.max(ut)
    # draw a new sample from the ground truth
    counts          <- c(0, 0)
    counts[gt(x[k])] <- 1
    add.measurement(experiment, x[k], counts)
  }
}

new.gt <- function(x, y)
{
  gt <- function(xt) {
    if (runif(1, 0, 1) <= y[x == xt]) {
      return (1)
    }
    else {
      return (2)
    }
  }
  return (gt)
}
