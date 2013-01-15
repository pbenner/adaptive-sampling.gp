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
