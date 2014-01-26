
require(adaptive.sampling.gp)

# data
################################################################################

# x a matrix of coordinates where each row
# contains [theta phi]
# theta: polar angle     ( latitude) [0,  pi]
# phi  : azimuthal angle (longitude) [0, 2pi]
x <- c()
# the first two points should result in the maximal distance
x <- rbind(x, c(1/2*pi, 0.0*2*pi))
x <- rbind(x, c(1/2*pi, 1/2*2*pi))
# a point almost exactly in-between
x <- rbind(x, c(1/4*pi, 1/(4+0.1)*2*pi))

# sperical kernel
################################################################################

kernel <- kernel.exponential(1, 0.5, spherical=T)
kernel(x)
