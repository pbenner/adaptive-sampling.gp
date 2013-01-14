
require("mvtnorm")
require("Matrix")

new.gp <- function(x, mu.prior, kernelf, sigma=NULL)
{
  mu <- rep(mu.prior[1], length(x))

  if (is.null(sigma)) {
    sigma <- as.matrix(nearPD(kernelf(x, x))$mat)
  }
  
  gp <- list(x        = x,              # where to evaluate the gp
             kernelf  = kernelf,        # kernel function
             mu.prior = mu.prior,       # prior mean
             mu       = mu,             # mean
             sigma    = sigma)
  class(gp) <- "gp"

  return (gp)
}

partially.apply <- function(f, ...) {
  capture <- list(...)
  function(...) {
    do.call(f, c(capture, list(...)))
  }
}

kernel.exponential <- function(l, var)
{
  f <- function(x, y) {
    n <- length(x)
    m <- length(y)
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

posterior <- function(gp, ...)
{
  UseMethod("posterior")
}

posterior.gp <- function(gp, xp, yp, noise=NULL)
{
  if (length(xp) > 0) {
    k0 <- gp$kernelf(xp, xp)     # K(X , X )
    k1 <- gp$kernelf(xp, gp$x)   # K(X , X*)
    k2 <- t(k1)                  # K(X*, X )
    k3 <- gp$kernelf(gp$x, gp$x) # K(X*, X*)

    mu <- gp$mu.prior

    # add noise to measurements?
    if (is.null(noise)) {
      A <- k0
    }
    else if (length(noise) == 1) {
      A <- k0 + noise
    }
    else {
      A <- k0 + diag(noise)
    }
    L        <- chol(A)
    Linv     <- solve(L)

    gp$mu    <- drop(mu + (k2 %*% Linv) %*% (t(Linv) %*% (yp - mu)))
    gp$sigma <- k3 - (k2 %*% Linv) %*% (t(Linv) %*% k1)

    # the resulting matrix is usually not positive definite due to
    # numerical errors; use nearPD to compute the nearest positive
    # definite matrix
    gp$sigma <- as.matrix(nearPD(gp$sigma)$mat)
  }
  return (gp)
}

samples <- function(gp, ...)
{
  UseMethod("samples")
}

samples.gp <- function(gp, n=1)
{
  x <- rmvnorm(n=n, mean=gp$mu, sigma=gp$sigma, method="chol")

  return (x)
}

bound <- function(line, range=c(0,1))
{
  # bound function values in the inverval given by range
  tmp0 <- sapply(line, function(x) min(x, range[2]))
  tmp1 <- sapply(tmp0, function(x) max(x, range[1]))
  return (tmp1)
}

plot.gp <- function(gp, s=NULL)
{
  col <- rgb(8/255, 81/255, 156/255, alpha=0.625)
  
  plot(gp$x, gp$mu, 'n', xlab="x", ylab="p", ylim=c(0,1))

  var <- diag(gp$sigma)
  z1  <- bound(gp$mu + 2*sqrt(var))
  z2  <- bound(gp$mu - 2*sqrt(var))
  
  polygon(c(gp$x, rev(gp$x)), c(z1, rev(z2)),
     col = col, border = NA)

  lines(gp$x, gp$mu, 'l', lwd=3)

  if (!is.null(s)) {
    # if this is a scalar, then it specifies the number of samples
    if (!is.matrix(s)) {
      s <- samples(gp, s)
    }
    # otherwise it contains already samples
    for (i in 1:dim(s)[1]) {
      lines(gp$x, bound(s[i,]), 'l', lwd=0.5)
    }
  }
}

# Example
################################################################################

gp <- new.gp(1:100/20, 0.0, kernel.exponential(10, 1))

xp <- c(1, 2, 3)
yp <- c(0.7, 0.7, 0.7)
# measurement noise
ep <- c(0.01, 0.01, 0.01)

gp <- posterior(gp, xp, yp, ep)

plot(gp)

plot(gp, samples(gp, 100))

# Experiment
################################################################################

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

# Information measures
################################################################################

kl.divergence <- function(gp, ...)
{
  UseMethod("kl.divergence")
}

kl.divergence.gp <- function(gp0, gp1)
{
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

# Example
################################################################################

e <- new.experiment(c(2.0,2.0))
#add.measurement(e, 1, c( 1,4))
#add.measurement(e, 2, c( 1,3))
gp0 <- posterior(e, 1:100/20)

add.measurement(e, 2, c( 0,1))
gp1 <- posterior(e, 1:100/20)

add.measurement(e, 2.5, c( 0,1))
gp2 <- posterior(e, 1:100/20)

kl.divergence(gp0, gp1)
kl.divergence(gp0, gp2)

# Sampling
################################################################################

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

# Example
################################################################################

e  <- new.experiment(c(2.0,2.0))
x  <- 1:10/2
gt <- new.gt(x,
             c(0.977895, 0.959606, 0.927331, 0.872769, 0.786948,
               0.666468, 0.523201, 0.388603, 0.307012, 0.327954))

sample.with.gt(e, x, gt)
plot(e, 1:100/20)
