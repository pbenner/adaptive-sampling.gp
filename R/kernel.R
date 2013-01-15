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

dyn.load("../src/gp.so")

kernel.exponential.c.1d <- function(x, y, l, var)
{
  storage.mode(x)   <- "double"
  storage.mode(y)   <- "double"
  storage.mode(l)   <- "double"
  storage.mode(var) <- "double"

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  .Call("exponential_kernel_1d", x, y, l, var)
}

kernel.exponential.c.2d <- function(x, y, l, var)
{
  storage.mode(x)   <- "double"
  storage.mode(y)   <- "double"
  storage.mode(l)   <- "double"
  storage.mode(var) <- "double"

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  .Call("exponential_kernel_2d", x, y, l, var)
}

kernel.exponential <- function(l, var)
{
  f <- function(x, y) {
    n <- dim(x)[1]
    m <- dim(y)[1]
    result      <- matrix(0.0, n*m)
    dim(result) <- c(n, m)
    for (i in 1:n) {
      for (j in 1:m) {
        result[i,j] <- var*exp(-1.0/(2.0*l^2)*(x[i]-y[j])^2)
      }
    }
    return (result)
  }
  return (f)
}

kernel.exponential <- function(l, var)
{
  f <- function(x, y=NULL) {
    if (is.null(y)) {
      y <- x
    }
    if (is.null(dim(x)) || dim(x)[2] == 1) {
      kernel.exponential.c.1d(x, y, l, var)
    }
    else if (dim(x)[2] == 2) {
      kernel.exponential.c.2d(x, y, l, var)
    }
    else {
      stop("x has wrong dimension")
    }
  }
  return (f)
}
