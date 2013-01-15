gp <- new.gp(1:100/20, 0.5, kernel.exponential(1, 1))

xp <- c(1, 2, 3)
yp <- c(0.7, 0.7, 0.7)
# measurement noise
ep <- c(0.01, 0.01, 0.01)

gp <- posterior(gp, xp, yp, ep)

plot(gp)

plot(gp, samples(gp, 100))
