
require("mvtnorm")

new.gp <- function(x, expectation, kernelf)
{
  # the prior expectation can be either a constant
  # or an array of the same length as x
  if (length(expectation) == length(x)) {
    mu <- expectation
  }
  else {
    mu <- rep(expectation[1], length(x))
  }
  
  gp <- list(x       = x,              # where to evaluate the gp
             kernelf = kernelf,        # kernel function
             mu      = mu,             # (prior) mean
             sigma   = kernelf(x, x))  # (prior) covariance
  class(gp) <- "gp"

  return (gp)
}

kernel.exponential <- function(variance, l)
{
  f <- function(x, y) {
    n <- length(x)
    m <- length(y)
    result      <- matrix(0.0, n*m)
    dim(result) <- c(n, m)
    for (i in 1:n) {
      for (j in 1:m) {
        result[i,j] <- variance*exp(-1.0/(2.0*l^2)*(x[i]-y[j])^2)
      }
    }
    return (result)
  }
  
  return (f)
}

kernel.bernoulli <- function(stddev, prime.x, prime.stddev, l)
{
  # stddev       : default standard deviation
  # prime.x      : a list of positions for which there are measurements
  # prime.stddev : standard deviation for each of the positions

  # if a given position xi exists in prime.x then output
  # its standard deviation, otherwise return the default deviation
  f.stddev <- function(xi) {
    if (xi %in% prime.x) {
      return (prime.stddev[prime.x == xi][1])
    }
    else {
      return (stddev)
    }
  }
  
  f <- function(x, y) {
    n <- length(x)
    m <- length(y)
    result      <- matrix(0.0, n*m)
    dim(result) <- c(n, m)
    for (i in 1:n) {
      for (j in 1:m) {
        result[i,j] <- f.stddev(x[i])*f.stddev(y[j])*exp(-1.0/(2.0*l^2)*(x[i]-y[j])^2)
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
  k0 <- gp$kernelf(xp, xp)     # K(X , X )
  k1 <- gp$kernelf(xp, gp$x)   # K(X , X*)
  k2 <- t(k1)                  # K(X*, X )
  k3 <- gp$kernelf(gp$x, gp$x) # K(X*, X*)

  mu <- 0.5

  # add noise to measurements?
  if (is.null(noise)) {
    A <- k0
  }
  else {
    A <- k0 + diag(noise)
  }
  L        <- chol(A)
  Linv1    <- solve(L)
  Linv2    <- solve(t(L))

  gp$mu    <- mu + (k2 %*% Linv1) %*% (Linv2 %*% (yp - mu))
  gp$sigma <- k3 - (k2 %*% Linv1) %*% (Linv2 %*% k1)

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

plot.gp <- function(gp, samples=NULL)
{
  col <- rgb(8/255, 81/255, 156/255, alpha=0.625)
  
  plot(gp$x, gp$mu, 'n', xlab="x", ylab="p", ylim=c(0,1))

  var <- diag(gp$sigma)
  z1  <- gp$mu + sqrt(var)
  z2  <- gp$mu - sqrt(var)
  
  polygon(c(gp$x, rev(gp$x)), c(z1, rev(z2)),
     col = col, border = NA)

  lines(gp$x, gp$mu, 'l', lwd=3)

  if (!is.null(samples)) {
    for (i in 1:dim(samples)[1]) {
      lines(gp$x, samples[i,], 'l', lwd=0.5)
    }
  }
}

# Example
################################################################################

gp <- new.gp(1:100/20, 0.0, kernel.exponential(1, 10))

gp <- new.gp(1:10/2, 0.0, kernel.bernoulli(0.1, c(1, 2), c(0.1, 0.1), 1))

xp <- c(1, 2, 3)
yp <- c(0.7, 0.7, 0.7)
# measurement noise
#ep <- c(0.01, 0.01, 0.01)
ep <- c(0.0, 0.0, 0.0)

gp <- posterior(gp, xp, yp, ep)

plot(gp)

plot(gp, samples(gp, 100))

# Experiment
################################################################################

new.experiment <- function(alpha)
{
  experiment        <- list(alpha = alpha,     # Dirichlet pseudo counts
                            data  = new.env()) # experimental data
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
  print (xp)
  print (yp)
  print (ep)
  # construct the gaussian process
  kernelf <- kernel.bernoulli(sqrt(prior.variance), xp, sqrt(ep), 1)
#  kernelf <- kernel.exponential(1, 5)
  gp      <- new.gp(x, prior.expectation, kernelf)
  gp      <- posterior(gp, xp, yp, ep)

  return (gp)
}

# Example
################################################################################

e <- new.experiment(c(1,1))
add.measurement(e, 1, c(10,3))
add.measurement(e, 2, c( 9,4))
add.measurement(e, 3, c(10,3))

e <- new.experiment(c(1,1))
add.measurement(e, 1, c( 1,2))
add.measurement(e, 2, c( 1,2))
add.measurement(e, 3, c( 1,2))

gp <- posterior(e, 1:100/20)
plot(gp)

gp <- posterior(e, 1:10/2)
plot(gp)
