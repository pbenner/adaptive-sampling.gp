e <- new.experiment(kernelf=kernel.exponential(1.0,0.5), prior.mean=0.5)

add.measurement(e, 1, c( 5, 4))
add.measurement(e, 2, c( 1, 3))
add.measurement(e, 3, c( 2, 6))

gp1 <- posterior(e, 1:100/20)
gp2 <- gp1
gp2$link  <- NULL
gp2$range <- NULL
p1 <- plot(gp1, main="Probability space", xlabel="x", ylabel="p")
p2 <- plot(gp2, main="Physical space", xlabel="x")
pdf(file="illustration.pdf")
grid.arrange(p1, p2, ncol=1)
dev.off()
