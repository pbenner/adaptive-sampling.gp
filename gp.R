
require("mvtnorm")

new.gp <- function(x, kernelf)
{
  gp <- list(x       = x,                   # where to evaluate the gp
             kernelf = kernelf,             # kernel function
             mu      = rep(0.5, length(x)), # (prior) mean
             sigma   = kernelf(x, x))       # (prior) covariance
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

posterior <- function(gp, ...)
{
  UseMethod("posterior")
}

posterior.gp <- function(gp, xp, yp, ep)
{
  k0 <- gp$kernelf(xp, xp)     # K(X , X )
  k1 <- gp$kernelf(xp, gp$x)   # K(X , X*)
  k2 <- t(k1)                  # K(X*, X )
  k3 <- gp$kernelf(gp$x, gp$x) # K(X*, X*)

  mean     <- 1/2

  A        <- k0 + diag(ep)
  L        <- chol(A)
  Linv     <- solve(L) %*% solve(t(L))

  gp$mu    <- mean + k2 %*% Linv %*% (yp - mean)
  gp$sigma <- k3   - k2 %*% Linv %*% k1
  
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
  z1  <- gp$mu + 2*sqrt(var)
  z2  <- gp$mu - 2*sqrt(var)
  
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

gp <- new.gp(0:100/20, kernel.exponential(1, 10))

xp <- c(1, 2, 3)
yp <- c(0.7, 0.7, 0.7)
ep <- c(0.01, 0.01, 0.01)

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

posterior.experiment <- function(experiment, x)
{
  gp <- new.gp(x, kernel.exponential(1, 10))
  xp <- c() # position
  yp <- c() # mean
  ep <- c() # variance

  for (key in ls(envir=experiment$data)) {
    x      <- as.numeric(key)
    alpha  <- experiment$alpha
    counts <- experiment$data[[key]]

    a0     <- sum(alpha + counts)
    a1     <- alpha[1] + counts[1]
    y      <- a1/a0
    e      <- a1*(a0 - a1)/(a0^2*(a0 + 1))

    xp <- append(xp, x)
    yp <- append(yp, y)
    ep <- append(ep, e)
  }
  gp <- posterior(gp, xp, yp, ep)

  return (gp)
}

# Example
################################################################################

e <- new.experiment(c(1,2))
add.measurement(e, 1, c(10,3))
add.measurement(e, 2, c( 9,4))
add.measurement(e, 3, c(10,3))

gp <- posterior(e, 1:100/20)
plot(gp)
