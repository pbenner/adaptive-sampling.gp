
kernel.exponential <- function(variance, l)
{
  f <- function(x, y) {
    n <- length(x)
    m <- length(y)
    result      <- matrix(0.0, n*m)
    dim(result) <- c(n, m)
    for (i in 1:n) {
      for (j in 1:m) {
        result[i,j] <- variance*exp(-1.0/(2.0*l^2)*norm(as.matrix(x[i]-y[j]), "F")^2)
      }
    }
    return (result)
  }
  
  return (f)
}

posterior.gp <- function(x, kernelf, xp, yp, ep)
{
  A     <- kernelf(xp, xp) + diag(ep)
  Ainv  <- solve(A)

  mean  <- 1/2
  mu    <- mean + kernelf(x, xp) %*% Ainv %*% (yp - mean)
  sigma <- kernelf(x, x ) - kernelf(x, xp) %*% Ainv %*% kernelf(xp, x)

  gp    <- list(x = x, mu = mu, sigma = sigma)
  class(gp) <- "gp"

  return (gp)
}

plot.gp <- function(gp)
{
  col <- rgb(8/255, 81/255, 156/255, alpha=0.625)
  
  plot(gp$x, gp$mu, 'n', xlab="x", ylab="p", ylim=c(0,1))

  var <- diag(gp$sigma)
  z1  <- gp$mu + 2*sqrt(var)
  z2  <- gp$mu - 2*sqrt(var)
  
  polygon(c(gp$x, rev(gp$x)), c(z1, rev(z2)),
     col = col, border = NA)

  lines(gp$x, gp$mu, 'l', lwd=3)
}

# Example
################################################################################

x <- 0:100/20
kernelf <- kernel.exponential(1, 1)
xp <- c(1, 2, 3)
yp <- c(0.7, 0.7, 0.7)
ep <- c(0.01, 0.01, 0.01)

gp <- posterior.gp(x, kernelf, xp, yp, ep)
plot(gp)
