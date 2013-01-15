
# 1-dimensional sampling
################################################################################
e  <- new.experiment(c(2.0,2.0))
x  <- 1:10/2
gt <- new.gt(x,
             c(0.977895, 0.959606, 0.927331, 0.872769, 0.786948,
               0.666468, 0.523201, 0.388603, 0.307012, 0.327954))

sample.with.gt(e, x, gt)
plot(e, 1:100/20)

# 2-dimensional sampling
################################################################################

# function that defines the gt
f <- function(z)
{
  if (is.vector(z)) {
    z <- t(as.matrix(z))
  }
  
  x <- z[,1]
  y <- z[,2]

  mux <- 4
  muy <- 3
  
  result <- exp(-((x-mux)^2 + (y-muy)^2 + (x-mux)*(y-muy)))

  return (result)
}

# generate an experiment
e           <- new.experiment(c(1.0,1.0))
# positions where to evaluate the Gaussian process
x           <- as.matrix(expand.grid(x = 0:20/4, y = 0:20/4))
# positions where samples can be drawn
x.sampling  <- as.matrix(expand.grid(x = 1:4, y = 1:4))

sample.with.gt(e, x.sampling, new.gt.f(f))
