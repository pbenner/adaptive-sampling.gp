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
