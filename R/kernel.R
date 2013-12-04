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

kernel.exponential.c.1d <- function(x, y, l, var)
{
  storage.mode(x)   <- "double"
  storage.mode(l)   <- "double"
  storage.mode(var) <- "double"
  if (!is.null(y)) {
    storage.mode(y)   <- "double"
  }

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.null(y) && !is.matrix(y)) {
    y <- as.matrix(y)
  }

  if (is.null(y)) {
    .Call("exponential_kernel_1d", x, x, l, var, PACKAGE="adaptive.sampling.gp")
  }
  else {
    .Call("exponential_kernel_1d", x, y, l, var, PACKAGE="adaptive.sampling.gp")
  }
}

kernel.exponential.c.2d <- function(x, y, l, var)
{
  storage.mode(x)   <- "double"
  storage.mode(l)   <- "double"
  storage.mode(var) <- "double"
  if (!is.null(y)) {
    storage.mode(y)   <- "double"
  }

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.null(y) && !is.matrix(y)) {
    y <- as.matrix(y)
  }

  if (is.null(y)) {
    .Call("exponential_kernel_2d", x, x, l, var, PACKAGE="adaptive.sampling.gp")
  }
  else {
    .Call("exponential_kernel_2d", x, y, l, var, PACKAGE="adaptive.sampling.gp")
  }
}

kernel.exponential.c.3d <- function(x, y, l, var)
{
  storage.mode(x)   <- "double"
  storage.mode(l)   <- "double"
  storage.mode(var) <- "double"
  if (!is.null(y)) {
    storage.mode(y)   <- "double"
  }

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.null(y) && !is.matrix(y)) {
    y <- as.matrix(y)
  }

  if (is.null(y)) {
    .Call("exponential_kernel_3d", x, x, l, var, PACKAGE="adaptive.sampling.gp")
  }
  else {
    .Call("exponential_kernel_3d", x, y, l, var, PACKAGE="adaptive.sampling.gp")
  }
}

kernel.exponential.spherical.c <- function(phi, theta, l, var)
{
  storage.mode(phi)   <- "double"
  storage.mode(l)   <- "double"
  storage.mode(var) <- "double"
  if (!is.null(theta)) {
    storage.mode(theta)   <- "double"
  }

  if (!is.matrix(phi)) {
    phi <- as.matrix(phi)
  }
  if (!is.null(theta) && !is.matrix(theta)) {
    theta <- as.matrix(theta)
  }

  if (is.null(theta)) {
    .Call("exponential_kernel_spherical", phi, phi, l, var, PACKAGE="adaptive.sampling.gp")
  }
  else {
    .Call("exponential_kernel_spherical", phi, theta, l, var, PACKAGE="adaptive.sampling.gp")
  }
}

kernel.exponential.c.1d.sparse <- function(x, y, l, var, n)
{
  storage.mode(x)   <- "double"
  storage.mode(l)   <- "double"
  storage.mode(var) <- "double"
  storage.mode(n)   <- "double"
  if (!is.null(y)) {
    storage.mode(y)   <- "double"
  }

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.null(y) && !is.matrix(y)) {
    y <- as.matrix(y)
  }

  if (is.null(y)) {
    nx <- 0
    ny <- min(dim(x)[1]-1, n)
    result <- .Call("exponential_kernel_1d_sparse", x, x, l, var, nx, ny, PACKAGE="adaptive.sampling.gp")
    bandSparse(dim(x)[1], dim(x)[1], 0:ny, result, symmetric=TRUE)
  }
  else {
    nx <- min(dim(x)[1]-1, n)
    ny <- min(dim(y)[1]-1, n)
    result <- .Call("exponential_kernel_1d_sparse", x, y, l, var, nx, ny, PACKAGE="adaptive.sampling.gp")
    bandSparse(dim(x)[1], dim(y)[1], -nx:ny, result, symmetric=FALSE)
  }
}

kernel.exponential.c.2d.sparse <- function(x, y, l, var, n)
{
  storage.mode(x)   <- "double"
  storage.mode(l)   <- "double"
  storage.mode(var) <- "double"
  storage.mode(n)   <- "double"
  if (!is.null(y)) {
    storage.mode(y)   <- "double"
  }

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.null(y) && !is.matrix(y)) {
    y <- as.matrix(y)
  }

  if (is.null(y)) {
    nx <- 0
    ny <- min(dim(x)[1]-1, n)
    result <- .Call("exponential_kernel_2d_sparse", x, x, l, var, nx, ny, PACKAGE="adaptive.sampling.gp")
    bandSparse(dim(x)[1], dim(x)[1], 0:ny, result, symmetric=TRUE)
  }
  else {
    nx <- min(dim(x)[1]-1, n)
    ny <- min(dim(y)[1]-1, n)
    result <- .Call("exponential_kernel_2d_sparse", x, y, l, var, nx, ny, PACKAGE="adaptive.sampling.gp")
    bandSparse(dim(x)[1], dim(y)[1], -nx:ny, result, symmetric=FALSE)
  }
}

#' Generate a squared exponential kernel
#' 
#' @param l length scale
#' @param var noise variance
#' @param n approximate the kernel with a symmetric sparse band matrix of n diagonals
#' @export

kernel.exponential <- function(l, var, n=NULL)
{
  f <- function(x, y=NULL) {
    # select an appropriate kernel function
    if (is.null(dim(x)) || dim(x)[2] == 1) {
      if (is.null(n)) {
        kernel.exponential.c.1d(x, y, l, var)
      }
      else {
        kernel.exponential.c.1d.sparse(x, y, l, var, n)
      }
    }
    else if (dim(x)[2] == 2) {
      if (is.null(n)) {
        kernel.exponential.c.2d(x, y, l, var)
      }
      else {
        kernel.exponential.c.2d.sparse(x, y, l, var, n)
      }
    }
    else if (dim(x)[2] == 3) {
      kernel.exponential.c.3d(x, y, l, var)
    }
    else {
      stop("x has wrong dimension")
    }
  }
  return (f)
}
