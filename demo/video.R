# 1-dimensional sampling
################################################################################
x        <- 0:100/20
theta    <- (1-tanh(-50:50/20))/2.5+0.1

e        <- new.experiment(c(1.0,1.0),
                           partially.apply(kernel.exponential, 0.8))
sample.x <- 0:10/2
gt       <- new.gt(x, theta)

for (i in 1:200) {
  sample.with.gt(e, sample.x, gt)
  png(filename=sprintf("plot_%03d.png", i),
      width = 1200, height = 800)
  plot(e, x)
  title(main = sprintf("%d samples", i))
  lines(x, theta, 'l', lwd=3, col='red')
  dev.off()
}

# mencoder mf://plot_*.png -mf type=png:fps=4 -ovc lavc -lavcopts
# vcodec=mpeg4 -oac copy -o plot.avi

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

for (i in 1:250) {
  sample.with.gt(e, x.sampling, new.gt.f(f))
  png(filename=sprintf("plot_%03d.png", i),
      width = 1200, height = 800)
  plot(e, x, f)
#  title(main = sprintf("%d samples", i))
  dev.off()
}

# mencoder mf://plot_*.png -mf type=png:fps=4 -ovc lavc -lavcopts
# vcodec=mpeg4 -oac copy -o plot.avi
