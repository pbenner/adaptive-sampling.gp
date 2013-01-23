
library(adaptive.sampling.gp)

# Experiment 1
################################################################################
e <- new.experiment()
add.measurement(e, 1, c(100,3))
add.measurement(e, 2, c( 90,1))
add.measurement(e, 3, c(100,4))

gp <- posterior(e, 1:100/20)
plot(gp)

# Experiment 2
################################################################################
e <- new.experiment(kernelf=kernel.exponential(1.0,0.1))
add.measurement(e, 1, c(1,40))
add.measurement(e, 2, c( 1,3))
add.measurement(e, 3, c( 1,6))

gp <- posterior(e, 1:100/20)
plot(gp, samples(gp, 10))

# Experiment 3
################################################################################

e <- new.experiment()
add.measurement(e, c(1,1), c( 1,4))
add.measurement(e, c(2,2), c( 1,3))
add.measurement(e, c(3,3), c( 1,6))

x  <- as.matrix(expand.grid(x = 1:20/4, y = 1:20/4))
gp <- posterior(e, x)
plot(gp)
